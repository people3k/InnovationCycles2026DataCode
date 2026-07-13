#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Null Model Ensemble Fitting + Plotting

Outputs saved to ./NASWNull/
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import innopop_MultiVarPars as inp

# ---------------------------
# PARAMETERS
# ---------------------------

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
    's': 0.25
}

var_par = {
    'A': (0.01, 10),
    'be': (0.005, 5),
    'phi': (0.2, 0.4)
}

MAX_ITERS = 1000
SCALING = 100
OBJECTIVES = ['shape'] 
# 'dtw' 'rmse',

# ---------------------------
# PATHS
# ---------------------------

BASE_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(BASE_DIR, "NASWNull")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------
# LOAD DATA
# ---------------------------

data = pd.read_csv('NARecessions.csv')
data['time'] = -data['calBP']
data = data.sort_values(['region_id', 'time']).reset_index(drop=True)

# KDE columns
kde_cols = [c for c in data.columns if c.startswith('V')]
kde_cols = sorted(kde_cols, key=lambda x: int(x[1:]))
kde_cols = kde_cols[:500]
kde_cols = ['MKDE'] + kde_cols

# ---------------------------
# FIT NULL MODEL
# ---------------------------

all_results = []

for obj in OBJECTIVES:
    print(f"\n===== OBJECTIVE: {obj} =====")

    for region in data['region_id'].unique():
        print(f"\nRegion: {region}")

        df_reg = data[data['region_id'] == region]

        for col in kde_cols:
            print(f"  Fitting {col}")

            series = df_reg[col].dropna()

            # FIX: use actual time values (not index)
            time_points = df_reg.loc[series.index, 'time'].values
            values = series.values * SCALING

            results = inp.simulate_model(
                data_series=values,
                time_points=time_points,
                fixed_params=fixed_params,
                change_points=[],
                var_par=var_par,
                max_iters=MAX_ITERS,
                algorithm='annealing',
                objective=obj
            )

            for r in results:
                r['region_id'] = region
                r['kde_id'] = col
                r['objective'] = obj

            all_results.extend(results)

df_all = pd.DataFrame(all_results)
df_all.to_csv(os.path.join(OUTPUT_DIR, "NullModel_AllFits_AllObjectives.csv"), index=False)

# ---------------------------
# COMPUTE 100% ENVELOPE (MIN/MAX)
# ---------------------------

df_envelope = (
    df_all
    .groupby(['region_id', 'time', 'objective'])
    .agg(
        Y_mean=('YPredicted', 'mean'),
        Y_min=('YPredicted', 'min'),
        Y_max=('YPredicted', 'max')
    )
    .reset_index()
)

df_envelope.to_csv(os.path.join(OUTPUT_DIR, "NullModel_Envelopes_AllObjectives.csv"), index=False)

# ---------------------------
# PLOTTING FUNCTION
# ---------------------------

def plot_objective_comparison(region):

    df_real = data[data['region_id'] == region].copy()
    df_real = df_real.sort_values('time')

    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

    for i, obj in enumerate(OBJECTIVES):

        ax = axes[i]

        df_obj = df_envelope[
            (df_envelope['region_id'] == region) &
            (df_envelope['objective'] == obj)
        ].sort_values('time')

        ax.plot(df_real['time'], df_real['MKDE'] * SCALING, label='Real (MKDE)', linewidth=2)
        ax.plot(df_obj['time'], df_obj['Y_mean'], label='Mean Fit', linewidth=2)

        # FIX: full min-max envelope
        ax.fill_between(
            df_obj['time'],
            df_obj['Y_min'],
            df_obj['Y_max'],
            alpha=0.3,
            label='Full Envelope (Min–Max)'
        )

        ax.set_title(f"{obj.upper()} Objective")
        ax.set_ylabel("Output (Y)")
        ax.legend()

    axes[-1].set_xlabel("Time")

    plt.suptitle(f"Null Model Comparison — {region}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    filename = f"NullModel_Comparison_{region}.pdf"
    plt.savefig(os.path.join(OUTPUT_DIR, filename))
    plt.close()

# ---------------------------
# GENERATE PLOTS
# ---------------------------

for region in data['region_id'].unique():
    print(f"Plotting {region}")
    plot_objective_comparison(region)

print(f"\nDONE. Outputs saved to: {OUTPUT_DIR}")
