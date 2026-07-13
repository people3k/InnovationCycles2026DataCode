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
    
    Here we present the main script based on functions written in innopop.py

@author: jbaggio
"""

# Usual Suspects
import numpy as np
import pandas as pd
import os
import itertools
import scipy as sp

#for graphing results
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm


#for clustering 
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import pairwise_distances, silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

#for the main analysis functions:
import innopop_MultiVars_Final as inp


'''

Functions to simulate the model, define the objectve function and fit the model tot he data

'''


'''
Script to run the functions
'''



# Fixed parameters
fixed_params = {
    'An': 1, 
    'd': 0.5, 
    'al': 0.5, 
    'tr': 0, 
    'depon': 1, 
    'gentime': 10, 
    'Z': 1, 
    'dn': 0.2, 
    'ksted': 1, 
    'nsted': 1, 
    'rsted': 1, 
    's' : 0.25,
}

# Bounds of varying parameters
var_par = {   
    'A': (0.01, 10),    
    'be': (0.005, 5),
    'phi':(0.2, 0.4)}  

#parameter check no overlap between fixed and variables
overlap = set(fixed_params) & set(var_par)
if overlap:
    raise ValueError(f"Parameters cannot be both fixed and variable: {overlap}")
    

#load data
os.chdir ('/Users/freemanh/Library/CloudStorage/Dropbox/agg_macroecology/InnovationCycles2026DataCode/ModelCalibration/Data')

#load all kde
data = pd.read_csv('NARecessions.csv')

data['time']  = - data['calBP']
data = data.sort_values(['region_id', 'time']).reset_index(drop=True)

data_label = data.pivot_table(index='time', columns='region_id', values='PeriodID', aggfunc='first')   

os.chdir ('/Users/freemanh/Library/CloudStorage/Dropbox/agg_macroecology/InnovationCycles2026DataCode/ModelCalibration/FitResults_NAR')

#selected_regions = ['Levant2', 'North30']
#test_data = data[data['region_id'].isin(selected_regions)]
# data = test_data

#scaling for output of the original data
scaling = 100

# Fit the model to real data using only RMSE in the objective function 
print('Fitting time-varying parameters for RMSE')
fit_res_rmse = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 1500, algorithm = 'annealing', objective = 'rmse', detection = '1d', geom = False)
print('Fitting Null Models for RMSE')
fit_res_rmse_null = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 1500, algorithm = 'annealing', objective = 'rmse', null_model = True, detection = '1d', geom = False)

# Save dataframes of results in csv
df_fit_rmse = pd.DataFrame(fit_res_rmse)
df_null_rmse = pd.DataFrame(fit_res_rmse_null)
df_fit_rmse.to_csv('AllKDE_RMSE.csv')
df_null_rmse.to_csv('AllKDE_RMSE_Null.csv')

# Fit the model to real data using RMSE and Dynamic Time Warping in the objective function
print('Fitting time-varying parameters for DTW')
fit_res_dtw = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 200, algorithm = 'annealing', objective = 'dtw', detection = '1d', geom = False)
print('Fitting Null Models for DTW')
fit_res_dtw_null = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 200, algorithm = 'annealing', objective = 'dtw', null_model = True, detection = '1d', geom = False)

# Save dataframes of results in csv
df_fit_dtw = pd.DataFrame(fit_res_dtw)
df_null_dtw = pd.DataFrame(fit_res_dtw_null)
df_fit_dtw.to_csv('AllKDE_DTW.csv')
df_null_dtw.to_csv('AllKDE_DTW_Null.csv')

# Fit the model to real data using RMSE, Dynamic Time Warping and Time-Series shape in the objective function
print('Fitting time-varying parameters for SHAPE')
fit_res_hyb = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 1500, algorithm = 'annealing', objective = 'shape', detection = '1d', geom = False)
print('Fitting Null Models for SHAPE')
fit_res_hyb_null = inp.fit_all_kde(data, fixed_params, var_par, max_iters = 1500, algorithm = 'annealing', objective = 'shape', null_model = True, detection = '1d', geom = False)

# Save dataframes of results in csv
df_fit_hyb = pd.DataFrame(fit_res_hyb)
df_null_hyb = pd.DataFrame(fit_res_hyb_null)
df_fit_hyb.to_csv('AllKDE_Hyb_Features.csv')
df_null_hyb.to_csv('AllKDE_Hyb_Null.csv')


'''
Results 
one can open the df_fit dataframes if too long to compute for the visualization procedures


'''

os.chdir ('/Users/freemanh/Library/CloudStorage/Dropbox/agg_macroecology/InnovationCycles2026DataCode/ModelCalibration/FitResults_NAR')

df_fit_hyb = pd.read_csv('AllKDE_Hyb_Features.csv')
df_null_hyb = pd.read_csv('AllKDE_Hyb_Null.csv')
df_fit_dtw = pd.read_csv('AllKDE_DTW.csv')
df_null_dtw = pd.read_csv('AllKDE_DTW_Null.csv')
df_fit_rmse = pd.read_csv('AllKDE_RMSE.csv')
df_null_rmse = pd.read_csv ('AllKDE_RMSE_Null.csv')

#regions are the same for all fit dataframes so we calculate it once
regions = np.unique(df_fit_rmse['region_id'])

period_ranges = inp.compute_period_ranges(data_label)

param_names = list(var_par.keys())

'''

Stability , simulation and classification Analysis

'''

# store results
all_results = []

dfall = [df_fit_rmse, df_fit_dtw, df_fit_hyb]
fit_types = ["rmse", "dtw", "hybrid"]

for df_fit, fit_type in zip(dfall, fit_types):
    for reg in regions:
        df_reg = df_fit[df_fit['region_id'] == reg].copy().reset_index()
        if df_reg.empty:
            continue

        # dominant segmentation based on ALL variable parameters
        idx_arrays = []
        for p in param_names:
            _, idx = np.unique(df_reg[p], return_index=True)
            idx_arrays.append(idx)

        max_len = max(len(idx) for idx in idx_arrays)
        longest_idx = [idx for idx in idx_arrays if len(idx) == max_len][0]
        longest_idx_sorted = np.sort(longest_idx)

        # segment durations
        segment_starts = longest_idx_sorted
        segment_ends = np.append(longest_idx_sorted[1:], len(df_reg) - 1)
        times = df_reg.index
        t_lengths = np.abs(times[segment_starts] - times[segment_ends])

        # extract segment rows
        df_segments = df_reg.iloc[longest_idx_sorted].copy()
        if df_segments.empty:
            continue

        # variable parameter values per segment
        param_values = {p: df_segments[p].values for p in param_names}

        R_init = df_segments['RPredicted'].values
        K_init = df_segments['KPredicted'].values
        N_init = df_segments['NPredicted'].values

        # loop over segments 
        for ipar in range(len(df_segments)):

            # full parameter dictionary (fixed + variable)
            params = dict(fixed_params)
            for p in param_names:
                params[p] = param_values[p][ipar]

            # initial condition
            x0_guess = np.array(
                [R_init[ipar], K_init[ipar], N_init[ipar]],
                dtype=float
            )

            seg_duration = t_lengths[ipar]

            # base record
            result = {
                "fit_type": fit_type,
                "region_id": reg,
                "segment_index": ipar,
                **{p: params[p] for p in param_names},
                "R_init": R_init[ipar],
                "K_init": K_init[ipar],
                "N_init": N_init[ipar],
            }

            #  positive equilibrium check (model‑specific logic preserved) 
            A = params.get('A', 0.0)
            s = params.get('s', fixed_params.get('s', 0.0))
            be = params.get('be', 1.0)
            phi = params.get('phi', fixed_params.get('phi', 0.0))

            pos_condition = 5 * phi * (1 + A * np.sqrt(2 * s / (5 * phi)))

            if pos_condition > 1:
                result["pos_eq"] = "yes"
                r_eq = 1 / pos_condition
                y_eq = (np.tanh(20 * r_eq) * r_eq * (1 - r_eq)) / be
                k_eq = 2 * s * y_eq
                n_eq = 5 * phi * y_eq

                result.update({
                    "R_eq": r_eq,
                    "Y_eq": y_eq,
                    "K_eq": k_eq,
                    "N_eq": n_eq
                })

                x0_guess = np.array([r_eq, k_eq, n_eq], dtype=float)
            else:
                result.update({
                    "pos_eq": "no",
                    "R_eq": 1,
                    "Y_eq": 0,
                    "K_eq": 0,
                    "N_eq": 0
                })

            # fixed poin
            try:
                x_star, fp_info = inp.find_fixed_point_for_regime(params, x0=x0_guess)
                fp_success = True
                fp_message = fp_info.message
            except RuntimeError as e:
                fp_success = False
                fp_message = str(e)
                x_star = np.array([1, 0, 0], dtype=float)

            result.update({
                "fp_success": fp_success,
                "fp_message": fp_message,
                "R_star": x_star[0],
                "K_star": x_star[1],
                "N_star": x_star[2]
            })

            # Jacobian & stability
            if fp_success:
                try:
                    J = inp.numerical_jacobian(lambda t, x, p: inp.innovation_society(t, x, p), x_star, args=(params,))
                    
                    stability_label, eigvals, overshoot = inp.classify_stability(J, mode="bifurcation", return_overshoot=True)
                
                except Exception:
                    stability_label = "jacobian_failure"
                    eigvals = np.full(3, np.nan, dtype=complex)
                    overshoot = np.full(3, np.nan, dtype=complex)
            else:
                stability_label = "no_fp"
                eigvals = np.full(3, np.nan, dtype=complex)
                overshoot = np.full(3, np.nan, dtype=complex)

            result.update({
                "eig1_real": np.real(eigvals[0]),
                "eig1_imag": np.imag(eigvals[0]),
                "eig2_real": np.real(eigvals[1]),
                "eig2_imag": np.imag(eigvals[1]),
                "eig3_real": np.real(eigvals[2]),
                "eig3_imag": np.imag(eigvals[2]),
                "stability_label": stability_label,
                "overshoot_severity": overshoot
            })

            # trajectory classification
            if fp_success:
                try:
                    t_eval = np.linspace(0, seg_duration, max(2, seg_duration * 2))
                    sol = inp.simulate_trajectory(
                        params,
                        y0=x0_guess,
                        t_span=(t_eval[0], t_eval[-1]),
                        t_eval=t_eval,
                        max_step=0.5
                    )
                    if sol.success:
                        traj_label = inp.detect_fp_approach_extended(sol.y.T)
                    else:
                        traj_label = "sim_failure"
                except Exception:
                    traj_label = "sim_failure"
            else:
                traj_label = "no_fp"

            result["traj_label"] = traj_label

            all_results.append(result)

#Build final DataFrame
df_stab_sim = pd.DataFrame(all_results)
df_stab_sim['stab_traj'] = df_stab_sim['stability_label'] + ' | ' + df_stab_sim['traj_label']
df_stab_sim.to_csv('A_Stability_Trajectory.csv')


'''
        
stability analysis robustness to initial conditions (checking +- 25% of initial conditions)

'''
robust_results = []

# multipliers for ±25% variation in (R, K, N) initial conditions
scale_vals = np.linspace(0.75, 1.25, 5)
ic_multipliers = list(itertools.product(scale_vals, repeat=3))

for df_fit, fit_type in zip(dfall, fit_types):
    for reg in regions:
        print(f'Fit Type = {fit_type} and region = {reg}')
        df_reg = df_fit[df_fit['region_id'] == reg].copy()
        if df_reg.empty:
            continue

        # dominant segmentation based on ALL variable parameters
        idx_arrays = []
        for p in param_names:
            _, idx = np.unique(df_reg[p], return_index=True)
            idx_arrays.append(idx)

        max_len = max(len(idx) for idx in idx_arrays)
        longest_idx = [idx for idx in idx_arrays if len(idx) == max_len][0]
        longest_idx_sorted = np.sort(longest_idx)

        # segment durations
        segment_starts = longest_idx_sorted
        segment_ends = np.append(longest_idx_sorted[1:], len(df_reg) - 1)
        times = df_reg.index
        t_lengths = np.abs(times[segment_starts] - times[segment_ends])

        # extract segment rows
        df_segments = df_reg.iloc[longest_idx_sorted].copy()
        if df_segments.empty:
            continue

        # parameter values per segment
        param_values = {p: df_segments[p].values for p in param_names}

        R_init = df_segments['RPredicted'].values
        K_init = df_segments['KPredicted'].values
        N_init = df_segments['NPredicted'].values

        # loop over segments
        for ipar in range(len(df_segments)):

            # full parameter dictionary
            params = dict(fixed_params)
            for p in param_names:
                params[p] = param_values[p][ipar]

            # baseline predicted state
            x0_baseline = np.array(
                [R_init[ipar], K_init[ipar], N_init[ipar]],
                dtype=float
            )

            seg_duration = t_lengths[ipar]

            # loop over perturbed ICs
            for ic_id, (r_mult, k_mult, n_mult) in enumerate(ic_multipliers):

                x0_pert = x0_baseline * np.array(
                    [r_mult, k_mult, n_mult], dtype=float
                )

                # base record
                result = {
                    "fit_type": fit_type,
                    "region_id": reg,
                    "segment_index": ipar,
                    **{p: params[p] for p in param_names},

                    "R_pred": R_init[ipar],
                    "K_pred": K_init[ipar],
                    "N_pred": N_init[ipar],

                    "ic_id": ic_id,
                    "R_mult": r_mult,
                    "K_mult": k_mult,
                    "N_mult": n_mult,
                    "R_init_ic": x0_pert[0],
                    "K_init_ic": x0_pert[1],
                    "N_init_ic": x0_pert[2],
                }

                # positive equilibrium check (unchanged logic)
                A = params.get('A', 0.0)
                s = params.get('s', fixed_params.get('s', 0.0))
                be = params.get('be', 1.0)
                phi = params.get('phi', fixed_params.get('phi', 0.0))

                pos_condition = 5 * phi * (1 + A * np.sqrt(2 * s / (5 * phi)))

                if pos_condition > 1:
                    r_eq = 1 / pos_condition
                    y_eq = (np.tanh(20 * r_eq) * r_eq * (1 - r_eq)) / be
                    k_eq = 2 * s * y_eq
                    n_eq = 5 * phi * y_eq

                    result.update({
                        "pos_eq": "yes",
                        "R_eq": r_eq,
                        "Y_eq": y_eq,
                        "K_eq": k_eq,
                        "N_eq": n_eq
                    })

                    x0_fp = np.array([r_eq, k_eq, n_eq], dtype=float)
                else:
                    result.update({
                        "pos_eq": "no",
                        "R_eq": 1,
                        "Y_eq": 0,
                        "K_eq": 0,
                        "N_eq": 0
                    })
                    x0_fp = x0_baseline

                # fixed point
                try:
                    x_star, fp_info = inp.find_fixed_point_for_regime(
                        params, x0=x0_fp
                    )
                    fp_success = True
                    fp_message = fp_info.message
                except RuntimeError as e:
                    fp_success = False
                    fp_message = str(e)
                    x_star = np.array([1, 0, 0], dtype=float)

                result.update({
                    "fp_success": fp_success,
                    "fp_message": fp_message,
                    "R_star": x_star[0],
                    "K_star": x_star[1],
                    "N_star": x_star[2]
                })

                # Jacobian & stability
                if fp_success:
                    try:
                        J = inp.numerical_jacobian(
                            lambda t, x, p: inp.innovation_society(t, x, p),
                            x_star,
                            args=(params,)
                        )
                        stability_label, eigvals = inp.classify_stability(J)
                    except Exception:
                        stability_label = "jacobian_failure"
                        eigvals = np.full(3, np.nan, dtype=complex)
                else:
                    stability_label = "no_fp"
                    eigvals = np.full(3, np.nan, dtype=complex)

                result.update({
                    "eig1_real": np.real(eigvals[0]),
                    "eig1_imag": np.imag(eigvals[0]),
                    "eig2_real": np.real(eigvals[1]),
                    "eig2_imag": np.imag(eigvals[1]),
                    "eig3_real": np.real(eigvals[2]),
                    "eig3_imag": np.imag(eigvals[2]),
                    "stability_label": stability_label
                })

                # trajectory classification
                if fp_success:
                    try:
                        t_eval = np.linspace(0, seg_duration, max(2, seg_duration * 2))
                        sol = inp.simulate_trajectory(
                            params,
                            y0=x0_pert,
                            t_span=(t_eval[0], t_eval[-1]),
                            t_eval=t_eval,
                            max_step=0.5
                        )
                        if sol.success:
                            traj_label = inp.detect_fp_approach_extended(sol.y.T)
                        else:
                            traj_label = "sim_failure"
                    except Exception:
                        traj_label = "sim_failure"
                else:
                    traj_label = "no_fp"

                result["traj_label"] = traj_label

                robust_results.append(result)

# build DataFrame & save to a NEW file (no overwrite)
df_robust_stab_sim = pd.DataFrame(robust_results)

# combined stability/trajectory label
df_robust_stab_sim["stab_traj"] = df_robust_stab_sim["stability_label"] + " | " + df_robust_stab_sim["traj_label"]

# Save to a new CSV file
df_robust_stab_sim.to_csv("Stablity_Robust_Region.csv", index=False)

# summarize per region as % of each stability/trajectory regime
df_valid_stab_sim = df_robust_stab_sim.copy()
# count occurrences of each stab_traj per (fit_type, region_id)
combo_counts = (
    df_valid_stab_sim
    .groupby(['fit_type', 'region_id', 'stab_traj'])
    .size()
    .rename('count')
    .reset_index()
)

# convert counts to percentages within each (fit_type, region_id)
combo_counts['pct'] = (
    combo_counts
    .groupby(['fit_type', 'region_id'])['count']
    .transform(lambda x: x / x.sum())
)

# wide-format table: rows = (fit_type, region_id), columns = stab_traj, values = percentage
wide_robust_stab_sim = combo_counts.pivot_table(
    index=['fit_type', 'region_id'],
    columns='stab_traj', 
    values='pct',
    fill_value=np.nan
)

# save summary 
wide_robust_stab_sim.to_csv("Stability_Robust_Percentages.csv")


#now do the same adding segment id
# Count occurrences within (fit_type, region_id, segment_index, stab_traj)
combo_counts_segment = (
    df_valid_stab_sim
    .groupby(['fit_type', 'region_id', 'segment_index', 'stab_traj'])
    .size()
    .rename('count')
    .reset_index()
)

# Convert counts to percentages *within each (fit_type, region_id, segment_id)*
combo_counts_segment['pct'] = (
    combo_counts_segment
    .groupby(['fit_type', 'region_id', 'segment_index'])['count']
    .transform(lambda x: x / x.sum())
)

# Wide-format table:
wide_robust_stab_sim_segment = combo_counts_segment.pivot_table(
    index=['fit_type', 'region_id'],
    columns=['segment_index', 'stab_traj'],
    values='pct',
    fill_value=np.nan
).sort_index(axis=1)  # optional: sort columns by segment then regime

# Optional: save
wide_robust_stab_sim_segment.to_csv("Stability_Robust_Segment_Percentages.csv")


'''

Figures 

'''


'''
Main Paper  Figures
'''

# Mapping from fit type to dataframes and output suffix
fit_configs = {
    "RMSE": {
        "fit_df": df_fit_rmse,
        "null_df": df_null_rmse,
        "suffix": "RMSE"
    },
    "DTW": {
        "fit_df": df_fit_dtw,
        "null_df": df_null_dtw,
        "suffix": "DTW"
    },
    "Hybrid": {
        "fit_df": df_fit_hyb,
        "null_df": df_null_hyb,
        "suffix": "Hybrid"
    }
}

# Store residual dataframes by fit type & region
kde_res = {
    "RMSE": {},
    "DTW": {},
    "Hybrid": {}
}

for fit_type, cfg in fit_configs.items():
    fit_df = cfg["fit_df"]
    null_df = cfg["null_df"]
    suffix = cfg["suffix"]

    print(f"Fit Type: {fit_type}")

    for reg in regions:
        print(reg)
        # Filter data for this region
        df_reg = fit_df[fit_df['region_id'] == reg].copy()
        df_reg_null = null_df[null_df['region_id'] == reg].copy()

        # Load label dataframe (calBP and period name)
        df_label = pd.DataFrame(data_label[reg].dropna())
        df_label2 = df_label.reset_index().rename(columns={'time': 'calBP', reg: 'period'})

        # Add calBP to null dataframe
        df_reg_null['calBP'] = df_label2['calBP'].values

        # In some cases the last percap is NaN (North30, South30, etc.) this is true only if using 2d to detect changes
        if len(df_label2) != len(df_reg):
            print('Missing Last PerCap: ' + str(reg))
            df_label2 = df_label2[:-1]
            df_reg['calBP'] = df_label2['calBP'].values
        else:
            df_reg['calBP'] = df_label2['calBP'].values

        # Merge period labels
        df_reg = df_reg.merge(df_label2, on='calBP', how='left')
        # Store for later KDE or other analysis
        kde_res[fit_type][reg] = df_reg
        
        # Assess change points
        
        # Get parameters and indices
        idx_arrays = []
        for p in param_names:
            _, idx = np.unique(df_reg[p], return_index=True)
            idx_arrays.append(idx)

        max_len = max(len(idx) for idx in idx_arrays)
        longest_idx = [idx for idx in idx_arrays if len(idx) == max_len][0]
        longest_idx_sorted = np.sort(longest_idx)
        
        # Convert indices to x-axis time values
        change_times = df_reg['calBP'].iloc[longest_idx_sorted].values

        # Plotting
        #calculate grid
        #  Plotting 
        n_params = len(param_names)
        n_rows = 2 + n_params
        
        fig = plt.figure(figsize=(14, 4 * n_rows))
        fig.suptitle(f'Region {reg} – Fit Type: {fit_type}', fontsize=16)

        gs = fig.add_gridspec(n_rows, 2, height_ratios=[1, 1] + [0.8] * n_params)

        # Y real vs predicted vs null
        ax1 = fig.add_subplot(gs[0, :])
        sns.lineplot(data=df_reg, x='calBP', y='YReal', label='Real', ax=ax1, errorbar=None)
        sns.lineplot(data=df_reg, x='calBP', y='YPredicted', label='Predicted', ax=ax1, errorbar=None)
        sns.lineplot(data=df_reg_null, x='calBP', y='YPredicted', label='Null', ax=ax1, errorbar=None)
        ax1.set_title('Predicted vs Real Output')
        ax1.set_xlabel('CalBP')
        ax1.set_ylabel('Output')
        inp.highlight_periods(ax1, period_ranges.get(reg, []))
        inp.highlight_change_points(ax1, change_times)
    
        # R, K, N trajectories
        ax2 = fig.add_subplot(gs[1, :])
        sns.lineplot(data=df_reg, x='calBP', y='RPredicted', label='R', ax=ax2, errorbar=None)
        sns.lineplot(data=df_reg, x='calBP', y='NPredicted', label='N', ax=ax2, errorbar=None)
        sns.lineplot(data=df_reg, x='calBP', y='KPredicted', label='K', ax=ax2, errorbar=None)
        ax2.set_title('State Variables (R, K, N)')
        ax2.set_xlabel('CalBP')
        ax2.set_ylabel('Value')
        inp.highlight_periods(ax2, period_ranges.get(reg, []))
        inp.highlight_change_points(ax2, change_times)

        # One subplot per variable parameter
        for i, p in enumerate(param_names):
            ax = fig.add_subplot(gs[2 + i, :])
            sns.lineplot(
                data=df_reg,
                x='calBP',
                y=p,
                label=p,
                ax=ax,
                errorbar=None
            )
            ax.set_title(f'Parameter: {p}')
            ax.set_xlabel('CalBP')
            ax.set_ylabel(p)
            inp.highlight_periods(ax, period_ranges.get(reg, []))
            inp.highlight_change_points(ax, change_times)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'{suffix}_{reg}.pdf')        
        plt.show()
 

#stability figures: with full or bifurcation model

stability_type = 'bifurcation' #alternative is full

df_gr_stab = df_stab_sim.groupby('fit_type')

#Stability per region and segment
if stability_type == 'full':
    stability_order = [
        'stable_node', 'stable_spiral', 'saddle',
        'no_fp', 'collapse', 'indeterminate', 'jacobian_failure'
        ]  
else:
    stability_order = [
        'real', 'complex_stable', 'complex_unstable'
        ]
    
traj_order = ['direct_fp', 'hook_fp', 'spiral_fp',  'collapse', 'no_fp',]

# Build an ordered list of all combo labels (only those that actually occur)
combo_order = []
for stab in stability_order:
    for traj in traj_order:
        mask = (df_stab_sim['stability_label'] == stab) & (df_stab_sim['traj_label'] == traj)
        if mask.any():
            combo_order.append(f'{stab} | {traj}')

# Also include any rare combos that don't fit the pattern
for cl in df_stab_sim['stab_traj'].unique():
    if cl not in combo_order:
        combo_order.append(cl)

# Map each combo_label to an integer
combo_to_int = {label: i for i, label in enumerate(combo_order)}
int_to_combo = {i: label for label, i in combo_to_int.items()}

df_stab_sim['combo_id'] = df_stab_sim['stab_traj'].map(combo_to_int)

#  Build a custom color palette
# define stability colors 

if stability_type == 'full':
    stab_colors = {
        'stable_node'     : (0.25, 0.60, 0.30),  # green
        'stable_spiral'   : (0.15, 0.55, 0.55),  # teal
        'saddle'          : (0.85, 0.55, 0.10),  # orange
        'collapse'        : (0.80, 0.20, 0.20),  # red
        #these are all indicators failure to characterize local Jacobian stability
        'no_fp'           : (0.60, 0.60, 0.60),  # gray
        'indeterminate'   : (0.60, 0.60, 0.60),  # gray
        'jacobian_failure': (0.60, 0.60, 0.60),  # gray
        }
else: 
    stab_colors = {
        'real'             : (0.20, 0.40, 0.80),   # blue
        'complex_stable'   : ( 0.85, 0.15, 0.15),  # red
        'complex_unstable' : (0.75, 0.20, 0.75),   # magnta
        
        }


    # define per-trajectory colors
traj_colors = {
    'direct_fp': (0.90, 0.90, 0.30),  # yellowish
    'hook_fp'  : (0.95, 0.55, 0.55),  # pinkish
    'spiral_fp': (0.40, 0.80, 1.00),  # sky blue
    'no_fp'    : (0.70, 0.70, 0.70),  # gray
    'collapse' : (0.80, 0.20, 0.20),  # red
}

def parse_combo(label):
    stab, traj = [x.strip() for x in label.split('|')]
    return stab, traj

pivot_mat = (
    df_stab_sim
    .pivot_table(
        index='region_id',
        columns='segment_index',
        values='stab_traj',
        aggfunc='first'
    )
    .sort_index(axis=0)   # ensure regions sorted
    .sort_index(axis=1)   # ensure segments sorted
)

drawsplit = True
for fit_type, df_stab in df_gr_stab:
    pivot_mat = df_stab.pivot_table(
        index='region_id',
        columns='segment_index',
        values='combo_id',
        aggfunc='first'
    )

    fig, ax = plt.subplots(figsize=(10, 15))

    # Turn off automatic image
    ax.set_xlim(0, pivot_mat.shape[1])
    ax.set_ylim(0, pivot_mat.shape[0])
    ax.invert_yaxis()  # So row 0 appears at top

    # Draw each cell manually
    for i, region in enumerate(pivot_mat.index):
        for j, seg in enumerate(pivot_mat.columns):
            
            combo_id = pivot_mat.loc[region, seg]
            if pd.isna(combo_id):
                continue

            combo_label = int_to_combo[int(combo_id)]
            stab, traj = parse_combo(combo_label)

            # Lookup colors
            stab_color = stab_colors.get(stab, (0.2,0.2,0.2))
            traj_color = traj_colors.get(traj, (0.85,0.85,0.85))
            
            if drawsplit == True:
                inp.draw_split_cell(ax, i, j, stab_color, traj_color)
            else:
                inp.draw_diagonal_cell(ax, i, j, stab_color, traj_color)
                
    #draw lines for visualization
    # Horizontal lines for every region row
    for y in range(pivot_mat.shape[0] + 1):
        ax.plot([0, pivot_mat.shape[1]], [y, y],
                color='black', linewidth=0.2)

    # Vertical lines ONLY at segment boundaries
    for x in range(pivot_mat.shape[1] + 1):
        ax.plot([x, x], [0, pivot_mat.shape[0]],
                color='black', linewidth=0.8)

    # Axes labels
    ax.set_xticks(np.arange(pivot_mat.shape[1]) + 0.5)
    ax.set_xticklabels(pivot_mat.columns, rotation=90)
    ax.set_yticks(np.arange(pivot_mat.shape[0]) + 0.5)
    ax.set_yticklabels(pivot_mat.index)

    ax.set_xlabel("Segment Index")
    ax.set_ylabel("Region ID")
    ax.set_title(f"Stability and Trajectory — {fit_type}")
    
    # Build Stability legend (top-left triangle)
    stab_patches = [
        plt.Line2D([0], [0], marker='s', color='none',
               markerfacecolor=color, markersize=12, 
               label=label.replace('_', ' '))
        for label, color in stab_colors.items()
        ]

    legend1 = ax.legend(
        handles=stab_patches,
        title="Stability Class",
        loc='upper left',
        bbox_to_anchor=(1.02, 1.00),
        frameon=True
        )
    ax.add_artist(legend1)


    # Build Trajectory legend (bottom-right triangle)
    traj_patches = [
        plt.Line2D([0], [0], marker='s', color='none',
               markerfacecolor=color, markersize=12, 
               label=label.replace('_', ' '))
        for label, color in traj_colors.items()
        ]

    ax.legend(
        handles=traj_patches,
        title="Trajectory Type",
        loc='upper left',
        bbox_to_anchor=(1.0, 0.45),
        frameon=True
    )
    
    plt.tight_layout()
    plt.savefig(f"ARegion_Stablity_{fit_type}.pdf")
    plt.show()


'''
Figure for overshoot shaded areas with Y_real, Y_null and Y_fitted
'''

region_order = [
    "Continental US",
    "Southwest US",
    "N. Colorado Plat.",
    "S. Colorado Plat.",
    "Sonoran Desert",
    "Chihuahua Desert"
]

for fit_type, cfg in fit_configs.items():
    fit_df = cfg["fit_df"]
    null_df = cfg["null_df"]
    suffix = cfg["suffix"]

    print(f"Fit Type: {fit_type}")

    # create ONE figure with 6 rows
    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(14, 22), sharex=False)
    axes = axes.flatten()

    for i, reg in enumerate(region_order):
        ax = axes[i]
        print(reg)

        # filter data
        df_reg = fit_df[fit_df['region_id'] == reg].copy()
        df_reg_null = null_df[null_df['region_id'] == reg].copy()

        if df_reg.empty:
            ax.set_visible(False)
            continue

        # calBP
        df_label = pd.DataFrame(data_label[reg].dropna())
        df_label2 = df_label.reset_index().rename(columns={'time': 'calBP', reg: 'period'})

        df_reg_null['calBP'] = df_label2['calBP'].values

        if len(df_label2) != len(df_reg):
            df_label2 = df_label2[:-1]
            df_reg['calBP'] = df_label2['calBP'].values
        else:
            df_reg['calBP'] = df_label2['calBP'].values

        df_reg = df_reg.merge(df_label2, on='calBP', how='left')

        # segmentation
        idx_arrays = []
        for p in param_names:
            _, idx = np.unique(df_reg[p], return_index=True)
            idx_arrays.append(idx)

        max_len = max(len(idx) for idx in idx_arrays)
        longest_idx = [idx for idx in idx_arrays if len(idx) == max_len][0]
        longest_idx_sorted = np.sort(longest_idx)

        change_times = df_reg['calBP'].iloc[longest_idx_sorted].values

        segment_starts = longest_idx_sorted
        segment_ends = np.append(longest_idx_sorted[1:], len(df_reg) - 1)

        for seg_id, (start, end) in enumerate(zip(segment_starts, segment_ends), start=1):
            df_reg.loc[start:end, 'segment_index'] = seg_id
            
        df_reg['segment_index'] = df_reg['segment_index'] - 1

        # merge overshoot
        df_reg = df_reg.merge(
            df_stab[['region_id', 'segment_index', 'overshoot_severity']],
            on=['region_id', 'segment_index'],
            how='left'
        )

        # RAW (no normalization)
        cmap = plt.cm.viridis
        vmin = df_stab['overshoot_severity'].min()
        vmax = df_stab['overshoot_severity'].max()

        # shade
        for seg_id in sorted(df_reg['segment_index'].dropna().unique()):
            seg_data = df_reg[df_reg['segment_index'] == seg_id]

            if seg_data.empty:
                continue

            severity = seg_data['overshoot_severity'].iloc[0]
            if pd.isna(severity):
                continue

            # raw scaling (global)
            color_val = (severity - vmin) / (vmax - vmin)

            ax.axvspan(
                seg_data['calBP'].min(),
                seg_data['calBP'].max(),
                color=cmap(color_val),
                alpha=0.30
            )
        
        # vertical dashed lines at segment starts
        for t in change_times:
            ax.axvline(
                x=t,
                color='black',
                linestyle='--',
                linewidth=1,
                alpha=0.7
                )   

        # lines
        sns.lineplot(data=df_reg, x='calBP', y='YReal', ax=ax, label='Real', errorbar=None)
        sns.lineplot(data=df_reg, x='calBP', y='YPredicted', ax=ax, label='Predicted', errorbar=None)
        sns.lineplot(data=df_reg_null, x='calBP', y='YPredicted', ax=ax, label='Null', errorbar=None)

        # labels
        ax.set_title(reg)
        ax.set_xlabel('')
        ax.set_ylabel('Output')

        #inp.highlight_periods(ax, period_ranges.get(reg, []))
        #inp.highlight_change_points(ax, change_times)

        # only show legend once
        if i == 0:
            ax.legend(loc='upper right')
        else:
            ax.get_legend().remove()

    # shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_clim(vmin, vmax)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
    cbar.set_label('Overshoot Severity')

    fig.suptitle(f'{fit_type} – Overshoot Across Regions', fontsize=16)

    plt.subplots_adjust(hspace=0.4, top=0.93, right=0.88)
    plt.savefig(f'Overshoot_{suffix}.pdf')
    plt.show()




      
'''
Supplementary Figures (more detailed)
'''

for group, df_stab in df_gr_stab:
    # Pivot: rows = traj_label, columns = stability_label
    gr_pivot = pd.pivot_table(
        data=df_stab,
        index='traj_label',          # rows
        columns='stability_label',   # columns
        values='region_id',          # anything non-null, we'll just count
        aggfunc='count',
        fill_value=0
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.heatmap(
        gr_pivot,
        ax=ax,
        cmap='coolwarm',
        annot=True,       # show counts or proportions
        fmt='g'           # 'g' for general format (works for int or float)
    )
    ax.set_title(f"Fit Type = {group}")
    ax.set_xlabel("Stability")
    ax.set_ylabel("Trajectory")
    plt.tight_layout()
    plt.savefig(f'ATraj_Stab_{group}.pdf')
    plt.show()
    
#robustness heatmap for stability|trajectory

# Ensure a nice ordering (optional)
wide_robust_stab_sim = wide_robust_stab_sim.sort_index()  # sorts by fit_type, then region_id

for fit_type in wide_robust_stab_sim.index.get_level_values('fit_type').unique():
    subset = wide_robust_stab_sim.xs(fit_type, level='fit_type')  # rows = region_id

    # Optionally, sort regimes (columns) by average prevalence
    col_order = subset.mean(axis=0).sort_values(ascending=False).index
    subset = subset[col_order]

    plt.figure(figsize=(max(8, 0.5 * subset.shape[1]), 0.35 * subset.shape[0] + 2))
    sns.heatmap(
        subset,
        cmap="coolwarm",
        annot=False,
        cbar_kws={'label': 'Proportion of IC×segment runs'},
        vmin=0.0,
        vmax=1.0
    )
    plt.title(f"Distribution of stability–trajectory regimes by region (%) — {fit_type}")
    plt.xlabel("Stability | Trajectory")
    plt.ylabel("Region ID")
    plt.tight_layout()
    plt.savefig(f'ARegion_FitRobustness_{fit_type}.pdf')
    plt.show()

#robustness heatmap for stability|trajecotry with segments

wide_robust_stab_sim_segment = wide_robust_stab_sim_segment.sort_index()  # sorts by fit_type, then region_id

for fit_type in wide_robust_stab_sim_segment.index.get_level_values('fit_type').unique():
    subset = wide_robust_stab_sim_segment.xs(fit_type, level='fit_type')  # rows = region_id

    # Optionally, sort regimes (columns) by average prevalence
    subset = subset.sort_index(axis=1, level=['segment_index', 'stab_traj'])

    plt.figure(figsize=(max(8, 0.5 * subset.shape[1]), 0.35 * subset.shape[0] + 2))
    sns.heatmap(
        subset,
        cmap="coolwarm",
        annot=False,
        cbar_kws={'label': 'Proportion of IC×segment runs'},
        vmin=0.0,
        vmax=1.0
    )
    plt.title(f"Distribution of stability–trajectory regimes by region (%) — {fit_type}")
    plt.xlabel("Stability | Trajectory")
    plt.ylabel("Region ID")
    plt.tight_layout()
    plt.savefig(f'ARegion_Segment_FitRobustness_{fit_type}.pdf')
    plt.show()
    

# Store residual dataframes by fit type & region (unchanged logic)
kde_res = {k: {} for k in fit_configs.keys()}

for reg in regions:
    print(f"\n=== Region: {reg} ===")

    # collect data for all fit types for this region
    reg_dfs = []

    for fit_type, cfg in fit_configs.items():
        fit_df = cfg["fit_df"]
        null_df = cfg["null_df"]
        suffix = cfg["suffix"]

        print(f"  Fit type: {fit_type}")

        # Filter data for this region
        df_reg = fit_df[fit_df['region_id'] == reg].copy()
        df_reg_null = null_df[null_df['region_id'] == reg].copy()

        # Label dataframe (calBP and period name)
        df_label = pd.DataFrame(data_label[reg].dropna())
        df_label2 = df_label.reset_index().rename(columns={'time': 'calBP', reg: 'period'})

        # Add calBP to null dataframe
        df_reg_null['calBP'] = df_label2['calBP'].values

        # Adjust for possible last NaN percap
        if len(df_label2) != len(df_reg):
            print('    Missing Last PerCap:', reg, 'for fit type:', fit_type)
            df_label2 = df_label2[:-1]
            df_reg['calBP'] = df_label2['calBP'].values
        else:
            df_reg['calBP'] = df_label2['calBP'].values

        # Merge period labels
        df_reg = df_reg.merge(df_label2, on='calBP', how='left')

        # Residuals & loss
        df_reg['residuals'] = df_reg['YReal'] - df_reg['YPredicted']

        # If you already have Y_Loss_rmse / Y_Loss_DTW / Y_Loss_Shape columns,
        # you can adapt this – here I create a generic Y_Loss column:
        if 'Y_Loss' not in df_reg.columns:
            # Example: squared error as loss – change if you have a different definition
            df_reg['Y_Loss'] = df_reg['residuals'] ** 2

        # Keep track of fit type for plotting
        df_reg['fit_type'] = fit_type

        # Store for later KDE or analysis
        kde_res[fit_type][reg] = df_reg

        # Append to list to build a combined dataframe for this region
        reg_dfs.append(df_reg)

    # Combine all fit types for this region
    reg_all = pd.concat(reg_dfs, ignore_index=True)

    # ---- Plotting: one figure per region with all fit types ----
    n_params = len(param_names)
    n_panels = 2 + n_params + 3   # loss, Y, params, K/N/R
    n_cols = 2
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(18, 4 * n_rows),
        sharex=True
        )
    fig.suptitle(f"Region {reg}: Comparison of Fit Types", fontsize=18)

    axes = axes.ravel()
    panel = 0

    # ---- 1) Loss by fit type ----
    ax = axes[panel]
    sns.lineplot(data=reg_all, x='calBP', y='loss', hue='fit_type',
             ax=ax, errorbar=None)
    ax.set_title('Loss by Fit Type')
    ax.set_ylabel('Loss')
    inp.highlight_periods(ax, period_ranges.get(reg, []))
    panel += 1

    # ---- 2) Y real vs predicted ----
    ax = axes[panel]
    real_df = (
        reg_all
        .drop_duplicates(subset=['calBP'])
        .sort_values('calBP')
    )
    sns.lineplot(data=real_df, x='calBP', y='YReal',
             color='black', label='YReal', ax=ax, errorbar=None)
    sns.lineplot(data=reg_all, x='calBP', y='YPredicted',
             hue='fit_type', ax=ax, errorbar=None)
    ax.set_title('Output: Real vs Predicted')
    ax.set_ylabel('Y')
    inp.highlight_periods(ax, period_ranges.get(reg, []))
    panel += 1

    # ---- 3) Variable parameters (dynamic) ----
    for p in param_names:
        ax = axes[panel]
        sns.lineplot(
            data=reg_all,
            x='calBP',
            y=p,
            hue='fit_type',
            ax=ax,
            errorbar=None
        )
        ax.set_title(f'Parameter: {p}')
        ax.set_ylabel(p)
        inp.highlight_periods(ax, period_ranges.get(reg, []))
        panel += 1

    # ---- 4) State trajectories ----
    for state, label in zip(
            ['KPredicted', 'NPredicted', 'RPredicted'],
            ['K – Infrastructure', 'N – Population', 'R – Resources']
            ):
        ax = axes[panel]
        sns.lineplot(
            data=reg_all,
            x='calBP',
            y=state,
            hue='fit_type',
            ax=ax,
            errorbar=None
            )   
        ax.set_title(label)
        ax.set_ylabel(state.replace('Predicted', ''))
        inp.highlight_periods(ax, period_ranges.get(reg, []))
        panel += 1

    # ---- Remove unused axes (if any) ----
    for ax in axes[panel:]:
        ax.set_visible(False)

    # ---- Legend handling ----
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', title='Fit Type')

    for ax in axes:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'Fit_Comparison_{reg}.pdf')
    plt.show()






'''
Cluster stability-trajectory for different regions.

'''

# Assuming df_stab_sim and df_gr_stab = df_stab_sim.groupby('fit_type') already exist

all_cluster_assignments = []

for fit_type, df_stab in df_gr_stab:
    # 1) Region × segment categorical matrix of 'stab_traj' strings
    pivot_combo = df_stab.pivot_table(
        index='region_id',
        columns='segment_index',
        values='stab_traj',   # "stability | trajectory"
        aggfunc='first'
    ).sort_index(axis=0).sort_index(axis=1)

    # compute Gower-like distance
    D = inp.gower_categorical_df(pivot_combo)

    # convert to condensed form for linkage
    D_condensed = squareform(D, checks=False)

    # hierarchical clustering
    Z = linkage(D_condensed, method='average')  # 'ward', 'complete', etc. also possible

    # dendrogram
    fig, ax = plt.subplots(figsize=(8, 4))
    dendrogram(
        Z, 
        labels=pivot_combo.index.astype(str).tolist(),
        leaf_rotation=90,
        ax=ax
    )
    ax.set_title(f"Region Clustering by stab_traj Pattern — {fit_type}")
    ax.set_ylabel("Distance")
    plt.tight_layout()
    plt.savefig(f'Dendograms_Stability-Trajectory_{fit_type}.pdf')
    plt.show()

    # cut the tree into K clusters (choose K based on dendrogram / domain knowledge)
    K = 0.5
    cluster_labels = fcluster(Z, K, criterion='distance')

    cluster_df = pd.DataFrame({
        'fit_type': fit_type,
        'region_id': pivot_combo.index,
        'cluster': cluster_labels
    })
    all_cluster_assignments.append(cluster_df)

# Combine clusters for all fit_types
clusters_all = pd.concat(all_cluster_assignments, ignore_index=True)

# Optional: merge back to df_stab_sim for further analysis / plotting
df_stab_sim = df_stab_sim.merge(
    clusters_all,
    on=['fit_type', 'region_id'],
    how='left'
)


params_cols = ['A', 'be', 's', 'phi']
df_gr_stab = df_stab_sim.groupby('fit_type')

for fit_type, df_stab in df_gr_stab:
    print(f"Plotting fit_type = {fit_type}")

    if df_stab.empty:
        continue

    
    # All unordered parameter pairs
    pairs = list(itertools.combinations(params_cols, 2))  # 6 pairs
    n_rows, n_cols = 2, 3

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 10))
    axes = axes.ravel()

    for ax, (x_param, y_param) in zip(axes, pairs):

        if df_stab.empty:
            ax.set_visible(False)
            continue

        # For each subplot, we’ll collect handles for legend
        legend_handles = []
        legend_labels = []
        last_scatter = None  # for colorbar
        sc = ax.scatter(
            df_stab[x_param],
            df_stab[y_param],
            c=df_stab["cluster"],
            cmap="tab10",
            s=90,
            alpha=0.85,
            edgecolor="k",
            marker= 'o' #marker_map[cl]
        )

        last_scatter = sc  # remember last one for colorbar
        # If no points plotted (e.g. all empty for some reason), hide axis
        if last_scatter is None:
            ax.set_visible(False)
            continue

        ax.set_xlabel(x_param)
        ax.set_ylabel(y_param)
        ax.set_title(
            f"{fit_type}: {y_param} vs {x_param}\n"
            "(color = cluster)"
        )

        # Colorbar for Y_eq
        cbar = fig.colorbar(last_scatter, ax=ax)
        cbar.set_label("Cluster Label")

        # Legend for clusters (marker shapes)
        ax.legend(
            legend_handles,
            legend_labels,
            fontsize=12,
            loc="best",
            title="Clusters"
        )

    plt.tight_layout()
    plt.show()


"""
Clustering on A, be, s and phi to assess different regimes.

Each cluster represents a type of regime (combination of A, be, s, phi).

"""


# Your existing helper to build segment table
df_segments_all = inp.extract_all_segments(dfall, fit_types, regions, param_names)

params_cols = param_names

segment_order_col = "segment_index"

# Kmeans n-cluster selection 

for fit_type, df_seg_ft in df_segments_all.groupby("fit_type"):

    # build region-level padded sequence embeddings for this fit_type
    X_reg_seq, region_ids, max_segments = inp.build_region_sequence_embeddings(
        df_seg_ft,
        params_cols=params_cols,
        segment_order_col=segment_order_col
    )

    n_regions = X_reg_seq.shape[0]

    # Standardize region-level embedding features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_reg_seq)

    # Explore K via elbow plot + silhouette scores
    max_K = min(12, n_regions - 1)
    K_range = list(range(2, max_K + 1))

    inertias = []
    sil_scores = []

    for k in K_range:
        kmeans_tmp = KMeans(
            n_clusters=k,
            random_state=0,
            n_init=10
        )
        labels_tmp = kmeans_tmp.fit_predict(X_scaled)
        inertias.append(kmeans_tmp.inertia_)

        if len(np.unique(labels_tmp)) > 1:
            sil_scores.append(silhouette_score(X_scaled, labels_tmp))
        else:
            sil_scores.append(np.nan)

    # Plot elbow and silhouette
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Elbow
    axes[0].plot(K_range, inertias, "-o")
    axes[0].set_xlabel("Number of clusters K")
    axes[0].set_ylabel("Inertia (within-cluster SSE)")
    axes[0].set_title(f"Elbow plot (regions, {fit_type})")
    axes[0].grid(True)

    # Silhouette
    axes[1].plot(K_range, sil_scores, "-o")
    axes[1].set_xlabel("Number of clusters K")
    axes[1].set_ylabel("Silhouette score")
    axes[1].set_title(f"Silhouette scores (regions, {fit_type})")
    axes[1].grid(True)

    plt.suptitle(f"Region-level KMeans model selection for fit_type = {fit_type}")
    plt.tight_layout()
    plt.show()


# Choose n_clusters based on visual inspection 

best_k_map = {
    "rmse": 6,
    "dtw": 5,
    "hybrid": 6
}

# store region-level cluster labels here (optional)
region_cluster_records = []

# use Kmeans found and do clustermaps

for fit_type, df_seg_ft in df_segments_all.groupby("fit_type"):

    best_k = best_k_map.get(fit_type, 4)

    # Build region-level padded sequence embeddings
    X_reg_seq, region_ids, max_segments = inp.build_region_sequence_embeddings(
        df_seg_ft,
        params_cols=param_names,
        segment_order_col=segment_order_col
    )

    n_regions = X_reg_seq.shape[0]

    # Standardize region-level embedding features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_reg_seq)

    # --- Region-level KMeans clustering ---
    kmeans_reg = KMeans(
        n_clusters=best_k,
        random_state=0,
        n_init=10
    )
    region_labels = kmeans_reg.fit_predict(X_scaled)

    #  Region × region distance matrix in standardized sequence space 
    dist_mat_reg = pairwise_distances(X_scaled, metric="euclidean")

    if inp.is_symmetric_exact(dist_mat_reg):
        dist_mat_reg = dist_mat_reg
    elif inp.is_symmetric_float(dist_mat_reg):
        # Enforce symmetry + zero diagonal (for safety against FP asymmetry)
        dist_mat_reg = 0.5 * (dist_mat_reg + dist_mat_reg.T)
        np.fill_diagonal(dist_mat_reg, 0.0)
    else:
        print('pairwise distance does not return a symmetric matrix even with floating point tolerance')


    # labels for regions as strings
    region_labels_str = [str(rid) for rid in region_ids]

    df_dist_reg = pd.DataFrame(
        dist_mat_reg,
        index=region_labels_str,
        columns=region_labels_str
    )

    # hierarchical linkage for ordering & dendrogram
    Y = squareform(dist_mat_reg, checks=False)
    Z = linkage(Y, method="average")

    # colors for clusters (row/col color bars)
    unique_clusters = np.unique(region_labels)
    n_clusters_ft = len(unique_clusters)

    palette = sns.color_palette("tab10", n_clusters_ft)
    cluster_color_map = {cl: palette[i] for i, cl in enumerate(unique_clusters)}

    row_colors = [cluster_color_map[cl] for cl in region_labels]
    col_colors = row_colors  # symmetric

    # clustermap: heatmap = region distances, colors = region clusters
    g = sns.clustermap(
        df_dist_reg,
        row_linkage=Z,
        col_linkage=Z,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap="RdYlBu_r",
        figsize=(12, 12)
    )

    g.savefig(f"KMeans_RegionSequence_Cluster_{fit_type}.pdf")
    plt.show()




'''
DTW Clustering based on A, be, s and phi

'''

# build sequences and DTW distance matrix
seq_ids, seqs = inp.build_region_sequences(df_segments_all, param_names)
dtw_mat = inp.build_region_distance_matrix(seq_ids, seqs, param_names)

# build labels and DataFrame
region_labels = [f"{reg}" for (ft, reg) in seq_ids]
df_dtw_mat = pd.DataFrame(dtw_mat, index=region_labels, columns=region_labels)

fit_types_unique = sorted(set(ft for (ft, reg) in seq_ids))



k = 3  # distance for hierarchical clustering

for ft in fit_types_unique:

    # indices of regions belonging to this fit type
    idx = [i for i, (ft_i, reg_i) in enumerate(seq_ids) if ft_i == ft]
    if len(idx) <= 1:
        continue

    # extract submatrix
    dtw_sub = dtw_mat[np.ix_(idx, idx)]
    region_labels_sub = [f"{seq_ids[i][1]}" for i in idx]

    # linkage on condensed distances
    Y_sub = squareform(dtw_sub)
    Z_sub = linkage(Y_sub, method='average')

    # cut dendrogram into clusters
    cluster_labels = fcluster(Z_sub, t=k, criterion='distance')
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    # map clusters to colors
    palette = sns.color_palette("tab20", n_clusters)
    cluster_color_map = {i: palette[i-1] for i in range(1, n_clusters+1)}
    row_colors = [cluster_color_map[c] for c in cluster_labels]
    col_colors = row_colors

    # dataframe for heatmap
    df_dtw_sub = pd.DataFrame(dtw_sub, index=region_labels_sub, columns=region_labels_sub)

    # find positive minimum for log scale
    positive_vals = dtw_sub[dtw_sub > 0]
    if positive_vals.size == 0:
        # fallback: skip log scaling if all zeros
        vmin = None
        norm = None
    else:
        vmin = positive_vals.min()
        vmax = dtw_sub.max()
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # clustermap with log-scaled colormap
    g = sns.clustermap(
        df_dtw_sub,
        row_linkage=Z_sub,
        col_linkage=Z_sub,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap="RdYlBu_r",
        norm=norm,          # <-- key line for log scaling
        figsize=(12, 12)
    )

    plt.suptitle(f"DTW Cluster for {ft}", y=1.02)
    g.savefig(f"DTW Cluster_{ft}.pdf")
    plt.show()
    

'''
potential stacked bar visualization for robustness and fit type for stability trajectories

import numpy as np
import matplotlib.pyplot as plt

# choose  fit_type
fit_type = "rmse"
subset = combo_pct_wide.xs(fit_type, level='fit_type')  # (n_regions, n_regimes)

regions_order = subset_top.index

fig, ax = plt.subplots(figsize=(10, 0.4 * len(regions_order) + 2))

bottom = np.zeros(len(regions_order))
colors = plt.cm.tab20(np.linspace(0, 1, len(top_regimes)))

for i, regime in enumerate(top_regimes):
    vals = subset_top[regime].values
    ax.barh(regions_order, vals, left=bottom, color=colors[i], label=regime)
    bottom += vals

ax.set_xlabel("Proportion of IC×segment runs")
ax.set_ylabel("Region ID")
ax.set_title(f"Top {top_k} stability–trajectory regimes by region — {fit_type}")
ax.legend(title="Stability | Trajectory", bbox_to_anchor=(1.05, 1.0), loc="upper left")
plt.tight_layout()
plt.show()


# Define a list of marker styles to cycle through
marker_styles = ['o', 's', '^', 'v', 'D', 'P', 'X', '*', 'h', '+', 'x']

for fit_type, df_stab in df_gr_stab:
    print(f"Plotting fit_type = {fit_type}")

    if df_stab.empty:
        continue

    # Get all clusters that actually appear for this fit_type
    clusters = (
        df_stab["cluster"]
        .dropna()
        .unique()
    )
    clusters = np.sort(clusters)  # for reproducibility

    if len(clusters) == 0:
        print(f"No cluster labels found for fit_type = {fit_type}, skipping.")
        continue

    # Map each cluster to a marker (cycle if more clusters than markers)
    marker_map = {
        cl: marker_styles[i % len(marker_styles)]
        for i, cl in enumerate(clusters)
    }

    # All unordered parameter pairs
    pairs = list(itertools.combinations(params_cols, 2))  # 6 pairs
    n_rows, n_cols = 2, 3

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 10))
    axes = axes.ravel()

    for ax, (x_param, y_param) in zip(axes, pairs):

        if df_stab.empty:
            ax.set_visible(False)
            continue

        # For each subplot, we’ll collect handles for legend
        legend_handles = []
        legend_labels = []
        last_scatter = None  # for colorbar

        for cl in clusters:
            sub = df_stab[df_stab["cluster"] == cl]
            if sub.empty:
                continue

            sc = ax.scatter(
                sub[x_param],
                sub[y_param],
                c=sub["Y_eq"],
                cmap="RdYlBu",
                s=90,
                alpha=0.85,
                edgecolor="k",
                marker= 'o' #marker_map[cl]
            )

            last_scatter = sc  # remember last one for colorbar

            # One legend entry per cluster (marker only, no color)
            # legend_handles.append(
            #     plt.Line2D(
            #         [0], [0],
            #         marker=marker_map[cl],
            #         linestyle="",
            #         color="k",
            #         markersize=7,
            #         markerfacecolor="none"
            #     )
            # )
            # legend_labels.append(f"cluster {cl}")

        # If no points plotted (e.g. all empty for some reason), hide axis
        if last_scatter is None:
            ax.set_visible(False)
            continue

        ax.set_xlabel(x_param)
        ax.set_ylabel(y_param)
        ax.set_title(
            f"{fit_type}: {y_param} vs {x_param}\n"
            "(color = Y_eq, marker = cluster)"
        )

        # Colorbar for Y_eq
        cbar = fig.colorbar(last_scatter, ax=ax)
        cbar.set_label("Y_eq")

        # Legend for clusters (marker shapes)
        ax.legend(
            legend_handles,
            legend_labels,
            fontsize=12,
            loc="best",
            title="Clusters"
        )

    plt.tight_layout()
    plt.show()

'''

    
    