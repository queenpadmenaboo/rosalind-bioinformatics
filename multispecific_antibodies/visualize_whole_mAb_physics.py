import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# --- 1. SETUP & DATA LOADING ---
DATA_DIR = "Whole_mAb"
physics_path = os.path.join(DATA_DIR, "Whole_mAb_3D_Physics_Features_Filtered.csv")
master_path = os.path.join(DATA_DIR, "Whole_mAb_Data_Sync_Master.csv")

if not os.path.exists(physics_path) or not os.path.exists(master_path):
    print(f"Error: Required CSV files not found in {DATA_DIR}.")
    exit()

df_p = pd.read_csv(physics_path)
df_m = pd.read_csv(master_path)

# Normalize and Merge
df_p['Therapeutic'] = df_p['Therapeutic'].str.strip().str.lower()
df_m['Therapeutic'] = df_m['Therapeutic'].str.strip().str.lower()
df = pd.merge(df_p, df_m, on='Therapeutic')
df['Isotype'] = df['CH1 Isotype_x'].astype(str)

df = df.loc[:, ~df.columns.str.endswith('_y')]
df.columns = df.columns.str.replace('_x', '')

# --- FEATURE DISTRIBUTION SUMMARY ---
print("\n" + "="*60)
print("FEATURE STATISTICS FOR ML MODEL:")
print(f"Total antibodies: {len(df)}")
print(f"\nHydro_SASA: mean={df['Hydro_SASA'].mean():.4f}, std={df['Hydro_SASA'].std():.4f}")
print(f"  Range: [{df['Hydro_SASA'].min():.4f}, {df['Hydro_SASA'].max():.4f}]")
print(f"\nNet_Charge: mean={df['Net_Charge'].mean():.2f}, std={df['Net_Charge'].std():.2f}")
print(f"  Range: [{df['Net_Charge'].min():.2f}, {df['Net_Charge'].max():.2f}]")
print("="*60 + "\n")

# --- 2. HISTOGRAM: Hydro_SASA Distribution ---
fig1 = px.histogram(
    df, 
    x='Hydro_SASA', 
    nbins=30,
    title="Physical Stickiness Distribution (n=300)",
    labels={'Hydro_SASA': 'Hydro_SASA (Surface Stickiness)'},
    template="plotly_dark",
    marginal="box"
)
mean_val = df['Hydro_SASA'].mean()
fig1.add_vline(x=mean_val, line_dash="dash", line_color="red", 
               annotation_text=f"Mean: {mean_val:.3f}")
fig1.write_html(os.path.join(DATA_DIR,"Whole_mAbs_Hydro_Histogram.html"))
print("Histogram saved: Whole_mAbs_Hydro_Histogram.html")

# --- 3. DEVELOPABILITY MAP: Scatter Plot ---
fig2 = px.scatter(
    df, 
    x='Net_Charge', 
    y='Hydro_SASA',
    color='Isotype',
    hover_name='Therapeutic',
    hover_data={'Hydro_SASA': ':.4f', 'Net_Charge': ':.2f', 'Isotype': True},
    title=f"Developability Map: Charge vs. Stickiness (n={len(df)})",
    labels={'Hydro_SASA': 'Hydro_SASA', 'Net_Charge': 'Net Charge (pH 7.4)'},
    template="plotly_dark"
)
fig2.update_traces(marker=dict(size=12, opacity=0.8, line=dict(width=1, color='white')))
fig2.update_layout(dragmode='zoom', hovermode='closest')
fig2.write_html(os.path.join(DATA_DIR,"Whole_mAbs_Developability_Map.html"))
print("Scatter plot saved: Whole_mAbs_Developability_Map.html")

# --- 4. ISOTYPE COMPARISON: Box Plot ---
fig3 = px.box(
    df, 
    x='Isotype', 
    y='Hydro_SASA',
    color='Isotype',
    points="all",
    title="Whole mAb Stickiness Comparison by Isotype",
    labels={'Hydro_SASA': 'Hydro_SASA', 'Isotype': 'Isotype Scaffold'},
    template="plotly_dark"
)
fig3.update_traces(marker=dict(size=6, opacity=0.5))
fig3.write_html(os.path.join(DATA_DIR,"Whole_mAbs_Isotype_Boxplot.html"))
print("Box plot saved: Whole_mAbs_Isotype_Boxplot.html")

# --- 5. FILTERING: EXTRACT STABLE CANDIDATES ---
# Data-driven thresholds based on your actual distribution
hydro_threshold = df['Hydro_SASA'].quantile(0.25)  # Bottom 25% (lowest stickiness)
charge_threshold = 2.0  # Reasonable charge threshold

stable_candidates = df[(df['Hydro_SASA'] < hydro_threshold) & (abs(df['Net_Charge']) > charge_threshold)].copy()
stable_candidates = stable_candidates.sort_values(by='Hydro_SASA')

output_csv = os.path.join(DATA_DIR, "Filtered_Stable_Whole_mAbs.csv")
stable_candidates.to_csv(output_csv, index=False)

# --- 6. DATA SUMMARY ---
print("\n" + "="*60)
print("--- ANALYSIS & FILTERING COMPLETE ---")
print("="*60)
print(f"Total Merged: {len(df)}")
print(f"Candidates Passing Stability Filter: {len(stable_candidates)}")
print(f"  (Hydro_SASA < {hydro_threshold:.4f} AND |Net_Charge| > {charge_threshold})")
print(f"Global Avg Hydro_SASA: {mean_val:.4f}")

if len(stable_candidates) > 0:
    print(f"\n--- STABLE CANDIDATES FOUND (Top 5) ---")
    print(stable_candidates[['Therapeutic', 'Isotype', 'Hydro_SASA', 'Net_Charge']].head(5).to_string(index=False))
else:
    print("\n--- NO CANDIDATES PASSED THE STABILITY FILTER ---")

print("="*60)
print(f"CSV SAVED: {output_csv}")
print("\nINTERACTIVE HTML FILES CREATED:")
print("  - Whole_mAbs_Hydro_Histogram.html")
print("  - Whole_mAbs_Developability_Map.html")
print("  - Whole_mAbs_Isotype_Boxplot.html")
print("="*60)