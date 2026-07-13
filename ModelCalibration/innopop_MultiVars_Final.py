#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:52:50 2025

This scripts uses the model developed by Marty Anderies in xppaut and moved into python in xppaut_arcaheology2 and
fits it to real world data.
The fit is done via simulated annealing.

The logic of the fitting is the following:
    A group of key parameters is fixed, as they are not-changeable by societies (i.e. intrinsic birth-death rates are physiological)
    The four key parameters society can alter are A , be, s and phi. 
    A is the technology infrastructure
    s is storage or the ability of society go generate surplus
    phi is the ability of society to convert suruplus into increased net population growth (either by increasing birth rate or reducing death rate)
    be is .... 
    
    We want to fit time-variant parameters as we do think that societies can alter those A, s, phi and be parameters. However, the variation must be within certain constraints, 
    unless there is a huge technnological innovation (still thinking about how to go for that...)

@author: jbaggio
"""

# Usual Suspects
import numpy as np
import scipy as sp
import pandas as pd

# Optimization algorithm and solver for the ODE
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.optimize import dual_annealing
from scipy.optimize import differential_evolution
from scipy.optimize import minimize
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import dendrogram

#time series decomposition via kernel or via topological decomposition (2d using MKDE and percap)
import ruptures as rpt
from scipy.signal import savgol_filter
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import RobustScaler


#for graphing results
import matplotlib.patches as patches

#for clustering time-series to assess A and be shape affecting declines
from fastdtw import fastdtw
from scipy.signal import find_peaks, periodogram
from numpy.linalg import norm as linorm



'''

Functions to simulate the model, define the objectve function and fit the model tot he data

'''

def assemble_params(fixed_params, var_par, x):
    """
    fixed_params : dict[str, float]
    var_par      : dict[str, (low, high)]
    x            : optimizer vector
    """
    params = dict(fixed_params)
    params.update(dict(zip(var_par.keys(), x)))
    return params


#scaling for a series in case values are too low for initial conditions to escape a 1, 0, 0 trap
scaling = 100

# Define the ODE system with fixed parameters per segment
def pow_clip(x, a):
    '''
    this function clips power at 0
    '''
    return x**a if x > 0 else 0

def gr1(x, t):
    return np.tanh(20 * (x - t))

def innovation_society(t, y, params):
    """
    params contains ALL parameter values, regardless of origin.
    """

    R, K, N = y

    # Just read values — NO role logic
    A      = params['A']
    be     = params['be']
    s      = params['s']
    phi    = params['phi']
    An     = params['An']
    d      = params['d']
    al     = params['al']
    tr     = params['tr']
    depon  = params['depon']
    Z      = params['Z']
    dn     = params['dn']
    ksted  = params['ksted']
    nsted  = params['nsted']
    rsted  = params['rsted']

    K_per_N = K / N if N > 0 else 0.0
    Y = N * (An + A * pow_clip(K_per_N, al)) * R

    dRdt = (depon * gr1(R, tr) * R * (Z - R) + (1 - depon) * R * (Z - R) - be * Y) * rsted
    dKdt = (s * Y - d * K) * ksted
    dNdt = (-dn * N + phi * Y) * nsted

    return [dRdt, dKdt, dNdt]




def detect_change_points(series, method="pelt_auto", max_breaks=10, penalties=None):
    """
    Detect change points using kernel RBF and PELT or BIC.
    """
    sizing = np.ceil(len(series) * 0.1) 
    algo = rpt.KernelCPD(kernel="rbf", min_size = sizing).fit(series)
    n = len(series)

    # BIC over k
    
    if method == "bic":
        bic_scores = []
        for n_bkps in range(1, max_breaks + 1):
            bkps = algo.predict(n_bkps=n_bkps)
            cost = algo.cost.sum_of_costs(bkps)
            bic = cost + n_bkps * np.log(n)
            bic_scores.append((n_bkps, bic, bkps))
        best_n_bkps, best_bic, best_bkps = min(bic_scores, key=lambda x: x[1])
        return np.unique(best_bkps)

    # PELT with theoretical penalty
    elif method == "penalty":
        penalty_value = 2 * np.log(max(n,2))  # common choice
        bkps = algo.predict(pen=penalty_value)
        return np.unique(bkps)

    # NEW: Automatic penalty sweep with PELT
    elif method == "pelt_auto":
        if penalties is None:
            base = np.log(max(n, 2))            
            penalties = [0.5 * base, 1.0 * base, 2.0 * base, 4.0 * base, 8.0 * base]

        results = []
        for pen in penalties:
            bkps = algo.predict(pen=pen)
            cost = algo.cost.sum_of_costs(bkps)
            bic = cost + len(bkps) * np.log(n)
            results.append((pen, bic, bkps))

        # Choose penalty giving lowest BIC
        best_pen, best_bic, best_bkps = min(results, key=lambda x: x[1])
        return np.unique(best_bkps)

    else:
        raise ValueError("Method must be 'bic', 'penalty', or 'pelt_auto'")

                
def geometric_features( x, y, t, eps=1e-3, smooth=True, smooth_window=5, smooth_polyorder=3, winsorize=True, winsor_limits=(0.01, 0.01)):
    
    """
    Compute geometric features for a 2D trajectory (x(t), y(t)):
    - speed
    - curvature
    - radial speed dr
    - angular speed dtheta
    - radius

    Parameters
    ----------
    x, y : Coordinates of the trajectory.
    t : Time stamps (not necessarily uniform, but must be increasing).
    eps             :Floor for curvature denominator.
    smooth          :If True, apply Savitzky-Golay smoothing to x and y before derivatives.
    smooth_window   :Window length for smoothing (must be odd and <= len(x)).
    smooth_polyord  :Polynomial order for Savitzky-Golay filter.
    winsorize       : If True, winsorize features per column (clip extreme quantiles).
    winsor_limits   :(lower, upper) quantile limits for winsorization, e.g. (0.01, 0.01).
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)

    # Ensure ordering by time (just in case)
    order = np.argsort(t)
    x, y, t = x[order], y[order], t[order]

    # Optional smoothing before derivatives
    if smooth and len(x) >= smooth_window and smooth_window % 2 == 1:
        x_s = savgol_filter(x, window_length=smooth_window,
                            polyorder=smooth_polyorder, mode="interp")
        y_s = savgol_filter(y, window_length=smooth_window,
                            polyorder=smooth_polyorder, mode="interp")
    else:
        x_s, y_s = x, y

    # First and second derivatives
    dx = np.gradient(x_s, t, edge_order=2)
    dy = np.gradient(y_s, t, edge_order=2)
    ddx = np.gradient(dx, t, edge_order=2)
    ddy = np.gradient(dy, t, edge_order=2)

    # Speed
    speed = np.sqrt(dx**2 + dy**2)

    # Curvature
    denom = (dx**2 + dy**2)**1.5
    denom = np.maximum(denom, eps)
    curvature = (dx * ddy - dy * ddx) / denom

    # Polar coordinates
    radius = np.sqrt(x_s**2 + y_s**2)
    dr = np.gradient(radius, t, edge_order=2)

    theta = np.unwrap(np.arctan2(y_s, x_s))
    dtheta = np.gradient(theta, t, edge_order=2)

    feats = np.column_stack([speed, curvature, dr, dtheta, radius])

    # Replace NaNs and infinities
    feats = np.nan_to_num(feats, nan=0.0, posinf=0.0, neginf=0.0)

    # Optional winsorization per feature (column-wise)
    if winsorize:
        lower_q = winsor_limits[0]
        upper_q = 1.0 - winsor_limits[1]
        for j in range(feats.shape[1]):
            q_low, q_high = np.quantile(feats[:, j], [lower_q, upper_q])
            feats[:, j] = np.clip(feats[:, j], q_low, q_high)

    return feats


def detect_2D_change_points(x, y, t, method="pelt_auto", penalties=None, geom = False, rbf_gamma='median', max_pairs_for_gamma = 200000000):
    """
    Detect regime shifts in a 2D trajectory using geometric features + KernelCPD.

    Parameters
    ----------
    x, y : Phase-plane trajectory.
    t : Time (for ordering and derivatives).
    method : {"pelt_auto", "penalty"}
        "pelt_auto" sweeps over penalties and selects by a BIC-like criterion.
        "penalty" uses a single default penalty.
    penalties : List of penalty values to sweep over for "pelt_auto".
    rbf_gamma : Gamma parameter for RBF kernel. If None, use kernel defaults.
    """
    
    if geom == True:
        # geometric features (with optional smoothing/winsorizing) detecting big regime changes
        features = geometric_features(x, y, t)
    else:
        # directly on data to detect less dramatic changes
        features = np.column_stack([x,y])

    # robustly scale features for kernel-based change detection
    features = RobustScaler().fit_transform(features)
    n = features.shape[0]
    
    if rbf_gamma == "median":
        # To avoid O(n^2) for very long series, subsample
        if n > max_pairs_for_gamma:
            idx = np.random.choice(n, size=max_pairs_for_gamma, replace=False)
            feats_sample = features[idx]
        else:
            feats_sample = features

        # pairwise distances on the sample
        dists = pairwise_distances(feats_sample, metric="euclidean")
        # take only upper triangle (excluding diagonal)
        iu = np.triu_indices_from(dists, k=1)
        dist_vals = dists[iu]
        # avoid zeros
        dist_vals = dist_vals[dist_vals > 0]

        if len(dist_vals) == 0:
            # fallback: use gamma = 1.0
            gamma_val = 1.0
        else:
            sigma = np.median(dist_vals)
            gamma_val = 1.0 / (2.0 * sigma**2)
    elif isinstance(rbf_gamma, (float, int)):
        gamma_val = float(rbf_gamma)
    elif rbf_gamma is None:
        gamma_val = None  # let ruptures / sklearn decide
    else:
        raise ValueError("rbf_gamma must be 'median', float, int or None.")

    # Prepare KernelCPD with optional gamma
    kernel_params = {}
    if rbf_gamma is not None:
        # passed as kwargs to sklearn.metrics.pairwise_kernels
        kernel_params["gamma"] = gamma_val

    sizing = np.ceil(len(x) * 0.1) 
    algo = rpt.KernelCPD(kernel="rbf", min_size = sizing, params = kernel_params).fit(features)

    if method == "pelt_auto":
        if penalties is None:
            # reasonable default grid scaled with log(n)
            base = np.log(max(n, 2))
            penalties = [0.5 * base, 1.0 * base, 2.0 * base, 4.0 * base, 8.0 * base]

        candidates = []
        for pen in penalties:
            bkps = algo.predict(pen=pen)
            # cost of this segmentation
            cost = algo.cost.sum_of_costs(bkps)
            # number of segments K = len(bkps)
            K = len(bkps)
            bic = cost + (K - 1) * np.log(n)
            candidates.append((pen, bic, bkps))

        best_pen, best_bic, best_bkps = min(candidates, key=lambda x: x[1])
        return np.unique(best_bkps)

    elif method == "penalty":
        # simple default: scaled with log(n)
        pen_val = 5.0 * np.log(max(n, 2))
        bkps = algo.predict(pen=pen_val)
        return np.unique(bkps)

    else:
        raise ValueError("method must be 'pelt_auto' or 'penalty'")

# Simulate the model segment-wise
def simulate_model(data_series, time_points, fixed_params, change_points, restart=True, var_par=None, max_iters=1000, algorithm='annealing', objective='rmse'):
    '''
    Note that the default algorithm is evolution as it allows for complex patterns fitting. However, at times it may does not converges even with a 
    high number of iteration (defined in fit_all_columns. Still, near optimal solutions can be valid especially if the loss, here RMSE is low enough.
    
    Parameters
    ----------
    data_series : the real world data
    time_points : time points determining data segmentation and where parameter values may change
    fixed_params : fixed parameters
    change_points : change points in the data detected
    param_bounds : lower and upper bound for varying parameters
    restart : whether we start the model at t=0 at each time-segement or continue. Does not make a difference as we start with initial conditions given by the previous endpoint.
    algorithm : the actual algorithm used, default annealing, two more options in case one want to experiment are differential evolution and l-bfgs.
    objective: what objective is minimized (either rmse = RMSE + derivative RMSE,  dtw = rmse + dtw  - dynamic time worping), and shape = rmse + corr + ts features)
    

    Returns
    -------
    outcome per time-step Y
    full_params outputs the time-varying per time-segment 
    full_loss outputs the loss function per time-segment
    '''
    
    segments = [0] + list(change_points) + [len(data_series)]
    full_results = []

    y0 = [1, 0, data_series[0]]

    if var_par is None:
        raise ValueError("var_par must be provided")

    # extract variable names + bounds ONCE
    var_names = list(var_par.keys())
    bounds = list(var_par.values())
    n_vars = len(var_names)

    for i in range(len(segments) - 1):
        start_idx = segments[i]
        end_idx = segments[i + 1]

        if restart:
            segment_length = end_idx - start_idx
            t_eval = np.arange(segment_length)
            t_span = (0, segment_length - 1)
            data_segment = data_series[start_idx:end_idx]
        else:
            t_eval = time_points[start_idx:end_idx]
            data_segment = data_series[start_idx:end_idx]
            t_span = (t_eval[0], t_eval[-1])

        if len(t_eval) < 2:
            continue

        # ---------- OPTIMIZATION ----------
        if algorithm == 'annealing':
            result = dual_annealing(
                segment_objective,
                bounds,
                args=(y0, t_span, t_eval, data_segment,
                      fixed_params, var_par, objective),
                maxiter=max_iters
            )

        elif algorithm == 'evolution':
            result = differential_evolution(
                segment_objective,
                bounds,
                args=(y0, t_span, t_eval, data_segment,
                      fixed_params, var_par, objective),
                strategy='best2bin',
                init='latinhypercube',
                popsize=25,
                maxiter=max_iters,
                tol=0.01,
                mutation=(0.4, 0.8),
                recombination=0.8,
                polish=True
            )

        else:  # L-BFGS-B
            lows = np.array([b[0] for b in bounds])
            highs = np.array([b[1] for b in bounds])
            x0 = np.random.uniform(lows, highs)

            result = minimize(
                segment_objective,
                x0,
                args=(y0, t_span, t_eval, data_segment,
                      fixed_params, var_par, objective),
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': max_iters}
            )

        if not result.success:
            print(f"Optimization warning for segment {start_idx}-{end_idx}: {result.message}")

        # build full param dict here
        params = assemble_params(fixed_params, var_par, result.x)

        # ODE receives dict 
        sol = solve_ivp(
            lambda t, y: innovation_society(t, y, params),
            t_span, y0, t_eval=t_eval, method='RK45'
        )

        if not sol.success or np.any(np.isnan(sol.y)) or np.any(np.isinf(sol.y)):
            print(f"Simulation failed for segment {start_idx}-{end_idx}. Skipping.")
            continue

        R, K, N = sol.y
        K_per_N = np.where(N > 0, K / N, 0.0)

        # NO positional indexing
        Y = N * (
            params['An']
            + params['A'] * np.maximum(K_per_N, 0)**params['al']
        ) * R

        loss = np.sqrt(np.mean((Y - data_segment) ** 2))

        for j, _ in enumerate(t_eval):
            full_results.append({
                'time': time_points[start_idx + j],
                'YPredicted': Y[j],
                'YReal': data_segment[j],
                'KPredicted': K[j],
                'NPredicted': N[j],
                'RPredicted': R[j],
                'loss': loss,
                **params  # ✅ all parameters, fixed or variable
            })

        y0 = [R[-1], K[-1], N[-1]]

    return full_results


def extract_shape_features(series, t_eval):
    """Extract shape-based features from a time series."""
    peaks, _ = find_peaks(series)
    valleys, _ = find_peaks(-series)
    amplitude = np.max(series) - np.min(series)
    mean_val = np.mean(series)
    std_val = np.std(series)
    freqs, power = periodogram(series, fs=1/(t_eval[1] - t_eval[0]))
    dominant_freq = freqs[np.argmax(power)] if len(freqs) > 0 else 0
    return np.array([len(peaks), len(valleys), amplitude, mean_val, std_val, dominant_freq])

#note: Before even running it, this objective is: 
    #the most informationally rich,
    #the most aligned with nonlinear feedback systems,
    #the only one accounting for oscillation count, dominant frequency, shape complexity,
    #the only one sensitive to non-monotonic dynamics.
    #hence hybrid/shape objective should be preferred

def segment_objective(x, y0, t_span, t_eval, data_segment, fixed_params, var_par, objective,
                      w_rmse=0.0, w_deriv_rmse=0.0,
                      w_dtw=0.0, w_corr=0.0, w_shape=0.0):
    """
    Unified hybrid loss objective for optimization.

    Each term is multiplied by a weight. If a weight is 0, that term is *not computed* at all (for efficiency), and contributes 0 to the loss.
    Also note that weights are theoretically justified and then we use dimensionless values in the loss. Values
    of the weights are given by 
    - RMSE must lead — it preserves basin membership
    - dRMSE helps but must be weaker (to avoid collapse)
    - Spearman pushes for global monotonicity or cyclic ordering
    - Shape features push for matching oscillation count and frequency
    - DTW enforces measures characteristics similar to derivative rmse but is more powerful, using both is akin to double counting.
    
    Note also that all parts of the loss are rendered dimensionless and scaled appropriately (theoretically)

    Parameters
    ----------
    x : result.x from simulated annealing
    fixed_params : array-like, parameters that are fixed
    var_par : array like parameters of the model to optimize.
    y0 : array-like, Initial conditions for ODE solver.
    t_span : tuple, Time span for integration (t0, tf).
    t_eval : array-like , Time points for evaluation.
    data_segment : array-like, Observed data segment.
    fixed_params : dict, Fixed model parameters; must contain keys 'An' and 'al'.
    objective: whether i want to have the loss function as rmse, dtw or shape, this changes the weights
    
    w_rmse : default 0.0       ,  Weight for RMSE term.
    w_deriv_rmse : default 0.0 , Weight for derivative RMSE term.
    w_dtw : default 0.0        , Weight for DTW term.
    w_corr : default 0.0       , Weight for (1 - Spearman correlation) term.
    w_shape : default 0.0      , Weight for feature-based (shape) distance term.

    Returns
    ------- 
    loss :   loss value to minimize. 
    """
   
    # assemble full parameter dict 
    params = assemble_params(fixed_params, var_par, x)

    # scaling
    sigma_y = np.std(data_segment)
    dy_true = np.gradient(data_segment, t_eval)
    sigma_dy = np.std(dy_true)
    len_ds = len(data_segment)

    sigma_y  = sigma_y  if sigma_y  > 1e-8 else 1.0
    sigma_dy = sigma_dy if sigma_dy > 1e-8 else 1.0

    # objective weights
    if objective == 'shape':
        w_rmse, w_deriv_rmse, w_corr, w_shape, w_dtw = 1.0, 0.4, 0.3, 0.25, 0.0
    elif objective == 'dtw':
        w_rmse, w_dtw = 1.0, 1.0
    else:  # rmse / hybrid
        w_rmse, w_deriv_rmse = 1.0, 0.3

    # solve ODE
    sol = solve_ivp(
        lambda t, y: innovation_society(t, y, params), t_span, y0, t_eval=t_eval, method='RK45' )

    if not sol.success or np.any(~np.isfinite(sol.y)):
        return np.inf

    R, K, N = sol.y
    K_per_N = np.where(N > 0, K / N, 0.0)

    #  model output
    Y = N * (params['An'] + params['A'] * np.maximum(K_per_N, 0)**params['al']) * R

    if Y.shape != data_segment.shape:
        return np.inf

    loss = 0.0

    # RMSE 
    if w_rmse:
        rmse = np.sqrt(np.mean((Y - data_segment) ** 2))
        loss += w_rmse * (rmse / sigma_y)

    # derivative RMSE 
    if w_deriv_rmse:
        dy_pred = np.gradient(Y, t_eval)
        rmse_d = np.sqrt(np.mean((dy_true - dy_pred) ** 2))
        loss += w_deriv_rmse * (rmse_d / sigma_dy)

    # DTW 
    if w_dtw:
        radius = min(len(Y), 10)
        dtw_score, _ = fastdtw(data_segment, Y, radius=radius, dist=2)
        loss += w_dtw * (dtw_score / max(len_ds, 1))

    # Spearman
    if w_corr:
        corr, _ = spearmanr(data_segment, Y)
        corr = 0.0 if np.isnan(corr) else corr
        loss += w_corr * ((1.0 - corr) / 2.0)

    # shape
    if w_shape:
        ft = extract_shape_features(data_segment, t_eval)
        fp = extract_shape_features(Y, t_eval)
        loss += w_shape * (linorm(ft - fp) / (linorm(ft) + 1e-8) / 6.0)

    return loss



#Fit data and record time-varying parameters
def fit_all_kde(df, fixed_params, var_par, max_iters=1000, algorithm='annealing', objective='hybrid', null_model=False, detection='2d', geom=False):

    '''
    Parameters
    ----------
    df : dataframe with real data and all N KDE models performed
    fixed_params : fixed parameters
    var_par : variable parameters with lower and upper bounds 
    max_iters : max iterations, The default is 1000.
    algorithm : type of optimization (annealing, lbfgsb or differential evolution. The default is 'annealing'.
    null_model: flag to fit the whole time-series or the segmented time series. Segmented is the default.

    Returns
    -------
    detailed_results : list of dictionaries with time-varying parameters and predictions
    '''
    detailed_results = []

    for region in np.unique(df['region_id']):
        print(f"Processing column: {region}")
        df_filt = df[df['region_id'] == region]
        #mean value of all the iterations and scale up by 1000
        orig_series = df_filt['MKDE'].dropna()
        series = orig_series.copy()
        time_points = series.index.to_numpy()
        values = series.values * scaling
        
        #assess change points based on full kde exploration
        coldrop = ['region_id', 'MKDE', 'StKDE', 'calBP', 'PerCap', 'PeriodID']
        df_change = df_filt.drop(columns=coldrop, errors='ignore')

                
        if null_model == True:
            changes = []
        else:
            #the following code does the following: looks at geometric fitures of KDE values and perCAP and then select breakpoints based on those features
            if detection == '2d':
                list_changes = []
                for col in df_change:
                    cols_needed = ['time', 'PerCap', col]
                    sub = df_filt[cols_needed].dropna()
                    time_values = sub['time'].values
                    percap_values = sub['PerCap'].values

                    kde_values = sub[col].values
                    
                    cp_kde = detect_2D_change_points(kde_values, percap_values, time_values, method='pelt_auto', geom = geom, rbf_gamma = 'median')
                    cp_kde = cp_kde[cp_kde < len(series)]
                    list_changes.append(cp_kde)
                
            else:
            #the following code does the following: calculates change point for each KDE, then select the median lenght
                list_changes = []
                for col in df_change:
                    kde_series = df_change[col].dropna()
                    kde_series = kde_series
                    kde_values = kde_series.values
                    cp_kde = detect_change_points(kde_values, method='pelt_auto')
                    cp_kde = cp_kde[cp_kde < len(series)]
                    list_changes.append(cp_kde)
            
            #filter out the median of change points based on only the number of KDE with the same number of break points (mode).
            #to avoid time-points not consistent with the index we use ceiling
            leng = [len(cps) for cps in list_changes]
            mode_leng, n_obs = sp.stats.mode(leng)
            filt_list = [cps for cps in list_changes if len(cps) == mode_leng]
            filt_mat = np.array(filt_list)
            changes = np.ceil(np.median(filt_mat, axis=0)).astype(int).tolist()
            print(f"Selected breakpoints: {len(np.unique(changes))}")

        results = simulate_model(
            data_series=values,
            time_points=time_points,
            fixed_params=fixed_params,
            change_points=changes,
            var_par=var_par,
            max_iters=max_iters,
            algorithm=algorithm,
            objective=objective
        )   

        # Ensure results is a list of dictionaries
        if isinstance(results, list):
            for res in results:
                if isinstance(res, dict):
                    res['region_id'] = region
                    detailed_results.append(res)
                else:
                    print(f"Warning: Unexpected result type for column {col}, skipping entry.")
        else:
            print(f"Warning: simulate_model did not return a list for column {col}, skipping.")

    return detailed_results


'''

Graphing functions

'''

# Function to compute period ranges for each region (assumes)
def compute_period_ranges(df_label):
    '''
    Parameters
    ----------
    df_label : dataframe that has time and period

    Returns
    -------
    period_ranges_by_region : start and end CalBP for the different periods
    '''
    
    period_ranges_by_region = {}

    for region in df_label.columns:
        region_series = df_label[region].dropna()
        region_series = region_series.sort_index(ascending=False)  # Ensure descending CalBP

        ranges = []
        current_period = None
        start_time = None
        prev_time = None

        for time, period in region_series.items():
            if period != current_period:
                if current_period is not None and prev_time is not None:
                    ranges.append((start_time, prev_time, current_period))
                current_period = period
                start_time = time
            prev_time = time

        if current_period is not None and prev_time is not None:
            ranges.append((start_time, prev_time, current_period))

        period_ranges_by_region[region] = ranges

    return period_ranges_by_region

# Function to highlight periods on graphs
def highlight_periods(ax, period_ranges):
    
    '''
    Parameters
    ----------
    ax : figure axs
    period_ranges : period_ranges calculated with compute_period_ranges
   
    Returns
    -------
    Highlights the different plot regions based on period.
    '''
    
    ymin, ymax = ax.get_ylim()
    ymid = (ymin + ymax) / 2
    colors = ['#e0e0e0', '#c0c0c0']

    for i, (start, end, name) in enumerate(period_ranges):
        left, right = min(start, end), max(start, end)
        ax.axvspan(left, right, facecolor=colors[i % 2], alpha=0.3)
        ax.text((left + right) / 2, ymid, name, ha='center', va='center', rotation=90, fontsize=14)
        
        
# Function to highlight change points
def highlight_change_points(ax, change_times, color='red', linestyle='--',
                            linewidth=1.5, alpha=0.8, zorder=5,
                            label=None):
    """
    Draw vertical lines at change points.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    change_times : array-like
        1D iterable of x-coordinates (time / CalBP) where change points occur.
    color : str
        Line color.
    linestyle : str
        Matplotlib line style (e.g. '--', '-.', ':').
    linewidth : float
        Line width.
    alpha : float
        Line transparency.
    zorder : int or float
        Z-order; higher values plotted on top.
    label : str or None
        Optional label for legend (will be used only for the first line).
    """
    change_times = np.asarray(change_times)
    if change_times.size == 0:
        return

    # Use label only on the first line so legend doesn't duplicate
    first = True
    for t in change_times:
        if first and label is not None:
            ax.axvline(x=t, color=color, linestyle=linestyle,
                       linewidth=linewidth, alpha=alpha,
                       zorder=zorder, label=label)
            first = False
        else:
            ax.axvline(x=t, color=color, linestyle=linestyle,
                       linewidth=linewidth, alpha=alpha,
                       zorder=zorder)

#draw triangle cells (diagonals)
def draw_diagonal_cell(ax, i, j, color_stability, color_trajectory):
    # Upper-left triangle (stability)
    tri1 = patches.Polygon(
        [(j, i+1), (j, i), (j+1, i+1)],
        closed=True, facecolor=color_stability, edgecolor='none'
    )
    ax.add_patch(tri1)

    # Bottom-right triangle (trajectory)
    tri2 = patches.Polygon(
        [(j+1, i), (j, i), (j+1, i+1)],
        closed=True, facecolor=color_trajectory, edgecolor='none'
    )
    ax.add_patch(tri2)

#draw split cells
def draw_split_cell(ax, row, col, stab_color, traj_color):
    # left half: stability
    rect1 = patches.Rectangle((col, row), 0.5, 1,
                              facecolor=stab_color, edgecolor='none')
    ax.add_patch(rect1)

    # right half: trajectory
    rect2 = patches.Rectangle((col + 0.5, row), 0.5, 1,
                              facecolor=traj_color, edgecolor='none')
    ax.add_patch(rect2)



'''

Stablity/Trajectory analysis functions

'''

#Numerical Jacobian
def numerical_jacobian(f, x_star, t=0.0, args=(), h=1e-6):
    """
    Compute numerical Jacobian of f at x_star using central differences.
    f(t, x, *args) -> array-like of length n.
    """
    x_star = np.asarray(x_star, dtype=float)
    n = x_star.size
    J = np.zeros((n, n))

    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        x_plus  = x_star + h * e
        x_minus = x_star - h * e

        f_plus  = np.asarray(f(t, x_plus, *args), dtype=float)
        f_minus = np.asarray(f(t, x_minus, *args), dtype=float)

        J[:, j] = (f_plus - f_minus) / (2.0 * h)

    return J


#Stability classification from eigenvalues
def classify_stability(J, tol=1e-10, mode="full", return_overshoot=False):
    """
    Classify local stability based on eigenvalues of the Jacobian.

    Parameters
    ----------
    J : ndarray, shape (n, n)
        Jacobian matrix evaluated at a fixed point.
    tol : float
        Tolerance for treating real parts as zero.

    Returns
    -------
    label : str
        One of 'stable_node', 'stable_spiral', 'unstable_node',
        'unstable_spiral', 'saddle', 'center', or 'indeterminate'.
    OR 
    label: str complex_unstable, complex_stable or real.
    
    eigvals : ndarray
        The eigenvalues of J.
    
    overshoot : ndarray, overshoot severity
    			
    	
    """
    eigvals = np.linalg.eigvals(J)

    real_parts = np.real(eigvals)
    imag_parts = np.imag(eigvals)

    has_imag = np.any(np.abs(imag_parts) > tol)

    # Bifurcation aware classification
    if mode == "bifurcation":
        
        # detect complex eigenvalues
        complex_mask = np.abs(imag_parts) > tol
        has_complex = np.any(complex_mask)

        #  classification 
        max_real = np.max(real_parts)
        
        if has_complex and max_real > tol:
            label = "complex_unstable"
        elif has_complex:
            label = "complex_stable"
        else:
            label = "real"

        #  overshoot computation 
        if return_overshoot:
            overshoots = []
            for r, im in zip(real_parts, imag_parts):
                val = im / (abs(r) + 0.001)
                val = val + 0.001
                overshoots.append(np.log10(val))

            # take maximum overshoot
            overshoot = np.nanmax(overshoots) if overshoots else np.nan

            return label, eigvals, overshoot
        else:
            return label, eigvals

    # Full classification

    neg = real_parts < -tol
    pos = real_parts > tol
    zero = np.abs(real_parts) <= tol
    has_imag = np.any(np.abs(imag_parts) > tol)

    if np.all(neg):
        label = 'stable_spiral' if has_imag else 'stable_node'
    elif np.all(pos):
        label = 'unstable_spiral' if has_imag else 'unstable_node'
    elif np.any(pos) and np.any(neg):
        label = 'saddle'
    elif np.all(zero) and has_imag:
        label = 'center'
    else:
        label = 'indeterminate'

    if return_overshoot:
        return label, eigvals, np.nan
    else:
        return label, eigvals

#Find a fixed point numerically
def find_fixed_point_for_regime(params, x0, t=0.0, tol=1e-9):
    """
    Find fixed point x* for a given parameter set.
    """
    #def f_wrapped(t, x, params, fixed_params):
     #   return innovation_society(t, x, params, fixed_params)

    def g(x):
        return np.asarray(innovation_society(t, x, params), dtype=float)

    sol = root(g, np.asarray(x0, dtype=float), tol=tol)
    
    if not sol.success:
        raise RuntimeError(f"Root finding failed: {sol.message}")
    return sol.x, sol



#Simulate trajectories
def simulate_trajectory(params, y0, t_span=(0, 200), t_eval=None, max_step=0.5, **kwargs):
    
    """
    Simulate dy/dt = innovation_society(t, y, params, fixed_params)
    """
    def rhs(t, y):
        return innovation_society(t, y, params)

    sol = solve_ivp(rhs, t_span=t_span, y0=y0, t_eval=t_eval, max_step=max_step, **kwargs)
    return sol


#Empirical classification of trajectories (distance-based)
def detect_fp_approach_extended(states,
                                tail_window=100,
                                eps=1e-6,
                                collapse_threshold=1e-6,
                                diverge_threshold=1e6):
    """
    Classify trajectory based on its long-term behavior.
    
    Uses only the TAIL of the trajectory to avoid false collapse detection
    when ICs start very close to zero for N or K (K always starts actually at 0).
    """

    states = np.asarray(states)
    T = states.shape[0]

    # Adjust tail window if the trajectory is short
    if T < tail_window:
        tail_window = T

    # Extract tail portion
    tail = states[-tail_window:]

    # --- Corrected collapse detection ---
    # Collapse only if the entire tail is near zero
    if np.all(tail < collapse_threshold):
        return "collapse"

    # Divergence detected from the tail
    if np.any(tail > diverge_threshold):
        return "diverge"

    # --- FP approach classification logic ---
    fp = tail.mean(axis=0)

    # Distances from FP
    d = np.linalg.norm(states - fp, axis=1)
    d = d[d > eps]  # Ignore tiny numerical noise

    if d.size < 10:
        return "direct_fp"

    dd = np.diff(d)
    signs = np.sign(dd)
    signs_nz = signs[signs != 0]

    if signs_nz.size < 2:
        return "direct_fp"

    sign_changes = np.sum(np.diff(signs_nz) != 0)

    if sign_changes >= 3:
        return "spiral_fp"
    elif np.any(dd > 0):
        return "hook_fp"
    else:
        return "direct_fp"
    


'''

Functions for Clustering

'''

#calculate distance for clustering stability-trajectory 

def gower_categorical_df(df_cat: pd.DataFrame) -> np.ndarray:
    """
    Compute a Gower-like distance matrix for a DataFrame of categorical values,
    allowing for missing values (NaN).
    
    Distance between two rows = 
        (# of differing columns, counting missing vs observed as differing)
        / (# of columns where at least one value is observed)
    """
    X = df_cat.values  # (n, p) array of objects / strings / NaN
    n, p = X.shape
    D = np.zeros((n, n), dtype=float)

    for i in range(n):
        xi = X[i, :]
        for j in range(i + 1, n):
            xj = X[j, :]

            # columns where at least one of xi,xj is not NaN
            mask = ~pd.isna(xi) | ~pd.isna(xj)

            if not mask.any():
                d = 0.0  # if both rows are entirely missing on all shared columns
            else:
                diffs = np.zeros(p, dtype=float)

                for k in range(p):
                    if not mask[k]:
                        continue

                    vi = xi[k]
                    vj = xj[k]

                    # Both missing --> shouldn't happen because mask excludes this case
                    if pd.isna(vi) and pd.isna(vj):
                        continue

                    # One missing, one not --> count as mismatch
                    if pd.isna(vi) or pd.isna(vj):
                        diffs[k] = 1.0
                    else:
                        diffs[k] = 0.0 if vi == vj else 1.0

                d = diffs[mask].sum() / mask.sum()

            D[i, j] = D[j, i] = d

    return D



def build_region_sequences(df_segments_all, param_names):
    """
    Build parameter sequences per (fit_type, region_id).
    """
    seq_ids = []
    seqs = []

    grouped = df_segments_all.groupby(["fit_type", "region_id"], sort=False)
    for (fit_type, region), g in grouped:
        g_sorted = g.sort_values("segment_index")

        seq_ids.append((fit_type, region))
        seqs.append({
            p: g_sorted[p].values.astype(float)
            for p in param_names
        })

    return seq_ids, seqs


def region_dtw_distance(seq_i, seq_j, param_names, radius=10, weights=None):
    """
    Robust multivariate DTW region-distance using fastdtw.
    """
    if weights is None:
        weights = {p: 1.0 for p in param_names}

    dist_sq = 0.0

    for p in param_names:
        x = np.asarray(seq_i[p], dtype=float)
        y = np.asarray(seq_j[p], dtype=float)

        if len(x) == 0 or len(y) == 0:
            continue

        d, _ = fastdtw(x, y, radius=radius, dist=2)
        dist_sq += weights.get(p, 1.0) * d**2

    return np.sqrt(dist_sq)


def build_region_distance_matrix(seq_ids, seqs, param_names, weights=None):
    """
    Build NxN DTW distance matrix for region sequences.
    """
    n = len(seqs)
    D = np.zeros((n, n), float)

    for i in range(n):
        for j in range(i + 1, n):
            d = region_dtw_distance(seqs[i], seqs[j], param_names=param_names, weights=weights)
            
            D[i, j] = D[j, i] = d

    return D


def get_branch_colors(Z, labels):
    """
    Returns a list of colors for leaves in the exact order of `labels`.
    """
    dend = dendrogram(Z, labels=labels, no_plot=True)
    leaf_order = dend["ivl"]
    leaf_colors = dend["leaves_color_list"]

    # Map: label → branch color
    color_map = dict(zip(leaf_order, leaf_colors))

    # Return colors in the original label order
    return [color_map[label] for label in labels]


def build_region_sequence_embeddings(df_seg_ft, params_cols, segment_order_col=None):
    """
    For a given fit_type subset (df_seg_ft), build region-level embeddings
    of sequences of parameter vectors, with padding to equal length.

    Each region r with segments in order:
        (A1, be1, s1, phi1), ..., (Ak, bek, sk, phik)
    is turned into a vector:
        [A1, be1, s1, phi1, A2, be2, s2, phi2, ..., Ak, bek, sk, phik, ... (padded)]

    Padding strategy:
        - Find max_segments across regions for this fit_type
        - If a region has fewer segments, repeat its last segment params
          until reaching max_segments.

    Parameters
    ----------
    df_seg_ft : DataFrame
        Subset of df_segments_all for a specific fit_type.
        Must contain columns: 'region_id', params_cols, and optionally segment_order_col.
    params_cols : list of str
        Names of parameter columns, e.g. ["A", "be", "s", "phi"].
    segment_order_col : str or None
        Column to sort segments within a region. If None, original index order is used.

    Returns
    -------
    X_reg_seq : ndarray, shape (n_regions, max_segments * len(params_cols))
        Region-level embedding matrix.
    region_ids : list
        List of region_ids corresponding to rows of X_reg_seq.
    max_segments : int
        Maximum number of segments per region for this fit_type.
    """

    # Number of segments per region
    seg_counts = df_seg_ft.groupby("region_id").size()
    max_segments = seg_counts.max()

    region_ids = []
    embeddings = []

    for region_id, df_reg in df_seg_ft.groupby("region_id"):
        if segment_order_col is not None and segment_order_col in df_reg.columns:
            df_reg_sorted = df_reg.sort_values(segment_order_col)
        else:
            # Fall back on the existing index order (assumed temporal)
            df_reg_sorted = df_reg.sort_index()

        seq = df_reg_sorted[params_cols].values  # shape (n_seg, n_params)
        n_seg, n_params = seq.shape

        if n_seg < max_segments:
            # Pad by repeating the last segment's parameters
            pad_val = seq[-1]  # shape (n_params,)
            pad = np.repeat(pad_val[None, :], max_segments - n_seg, axis=0)
            seq_padded = np.vstack([seq, pad])
        elif n_seg > max_segments:
            # Truncate if somehow there are more segments than max_segments
            seq_padded = seq[:max_segments, :]
        else:
            seq_padded = seq

        embedding = seq_padded.flatten()  # length max_segments * n_params
        embeddings.append(embedding)
        region_ids.append(region_id)

    X_reg_seq = np.vstack(embeddings)
    return X_reg_seq, region_ids, int(max_segments)

'''

Utilitiy to extract segments

'''


def extract_all_segments(dfall, fit_types, regions, param_names):
    """
    Build a long DataFrame with one row per (fit_type, region, segment),
    containing the segment-level parameter values.
    """
    segment_rows = []

    for df_fit, fit_type in zip(dfall, fit_types):
        for reg in regions:
            df_reg = df_fit[df_fit['region_id'] == reg].copy()
            if df_reg.empty:
                continue

            df_reg = df_reg.reset_index(drop=True)

            # dominant segmentation based on all parameters
            idx_arrays = []
            for p in param_names:
                _, idx = np.unique(df_reg[p], return_index=True)
                idx_arrays.append(idx)

            max_len = max(len(idx) for idx in idx_arrays)
            longest_idx = [idx for idx in idx_arrays if len(idx) == max_len][0]
            longest_idx_sorted = np.sort(longest_idx)

            df_segments = df_reg.iloc[longest_idx_sorted].copy()
            if df_segments.empty:
                continue

            for seg_idx, (_, seg_row) in enumerate(df_segments.iterrows()):
                row = {
                    "fit_type": fit_type,
                    "region_id": reg,
                    "segment_index": seg_idx,
                }
                for p in param_names:
                    row[p] = float(seg_row[p])

                segment_rows.append(row)

    df_segments_all = pd.DataFrame(segment_rows)
    return df_segments_all

#check symmetry for matrices

def is_symmetric_exact(matrix):
    arr = np.array(matrix)
    # Check if square and compare to its transpose
    if arr.shape[0] != arr.shape[1]:
        return False
    return (arr.T == arr).all() #
    

def is_symmetric_float(matrix, rtol=1e-05, atol=1e-08):
    arr = np.array(matrix)
    if arr.shape[0] != arr.shape[1]:
        return False
    # Compares two arrays for equality within a tolerance
    return np.allclose(arr, arr.T, rtol=rtol, atol=atol)
    



