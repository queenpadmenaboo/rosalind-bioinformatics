import pandas as pd
from Bio.PDB import PDBParser
from pathlib import Path

# Load and merge data
df_p = pd.read_csv("Whole_mAb/Whole_mAb_3D_Physics_Features_Filtered.csv")
df_m = pd.read_csv("Whole_mAb/Whole_mAb_Data_Sync_Master.csv")

df_p['Therapeutic'] = df_p['Therapeutic'].str.strip().str.lower()
df_m['Therapeutic'] = df_m['Therapeutic'].str.strip().str.lower()
df = pd.merge(df_p, df_m, on='Therapeutic')

# Check for isotype conflicts
conflicts = df[df['CH1 Isotype_x'] != df['CH1 Isotype_y']]
if len(conflicts) > 0:
    print(f"WARNING: {len(conflicts)} antibodies have mismatched isotypes!")
    print(conflicts[['Therapeutic', 'CH1 Isotype_x', 'CH1 Isotype_y']])
else:
    print("No conflicts - isotypes match perfectly")

print("="*70)
print("VALIDATION 1: DO CHAIN CHARGES SUM TO TOTAL NET_CHARGE?")
print("="*70)

# Check if individual chain charges sum to total
df['Chain_Sum'] = df['Ch_A_Charge'] + df['Ch_B_Charge'] + df['Ch_C_Charge'] + df['Ch_D_Charge']
df['Charge_Diff'] = abs(df['Net_Charge'] - df['Chain_Sum'])

# Flag discrepancies > 0.5 (allowing for rounding)
issues = df[df['Charge_Diff'] > 0.5]

if len(issues) > 0:
    print(f"FOUND {len(issues)} ANTIBODIES WITH CHARGE MISMATCHES:")
    print(issues[['Therapeutic', 'Net_Charge', 'Chain_Sum', 'Charge_Diff']].head(10))
else:
    print("PASS: All chain charges sum correctly to Net_Charge")

print(f"\nMax difference: {df['Charge_Diff'].max():.2f}")
print(f"Mean difference: {df['Charge_Diff'].mean():.4f}")

print("\n" + "="*70)
print("VALIDATION 2: RECALCULATE CHARGE FROM PDB FOR SAMPLE ANTIBODIES")
print("="*70)

# Recalculate charges for a few test antibodies
CHARGE_SCALE = {
    'ARG': 1.0, 'LYS': 1.0, 'HIS': 0.1,
    'ASP': -1.0, 'GLU': -1.0,
}

PDB_DIR = Path("Whole_mAb/PDB_Output_Files_GPU_Full")
parser = PDBParser(QUIET=True)

test_antibodies = ['trastuzumab', 'pembrolizumab', 'nivolumab']

for ab_name in test_antibodies:
    pdb_file = PDB_DIR / f"{ab_name}.pdb"
    
    if not pdb_file.exists():
        print(f"WARNING: {ab_name}: PDB file not found")
        continue
    
    # Recalculate from PDB
    struct = parser.get_structure(ab_name, str(pdb_file))
    recalc_charge = 0
    chain_charges = {}
    
    for model in struct:
        for chain in model:
            c_id = chain.get_id()
            c_charge = 0
            for res in chain:
                c_charge += CHARGE_SCALE.get(res.get_resname(), 0)
            chain_charges[c_id] = c_charge
            recalc_charge += c_charge
    
    # Get stored value
    stored = df[df['Therapeutic'] == ab_name]['Net_Charge'].values
    stored_val = stored[0] if len(stored) > 0 else None
    
    if stored_val is not None:
        diff = abs(recalc_charge - stored_val)
        status = "PASS" if diff < 0.5 else "FAIL"
        print(f"{status}: {ab_name:15s}: Stored={stored_val:6.1f}, Recalculated={recalc_charge:6.1f}, Diff={diff:.2f}")
        print(f"       Chain breakdown: {chain_charges}")
    else:
        print(f"WARNING: {ab_name}: Not found in dataset")

print("\n" + "="*70)
print("VALIDATION 3: CHARGE DISTRIBUTION SANITY CHECK")
print("="*70)

print(f"\nNet_Charge statistics (all 300 antibodies):")
print(df['Net_Charge'].describe())

print("\nBy Isotype:")
print(df.groupby('CH1 Isotype_x')['Net_Charge'].agg(['mean', 'std', 'min', 'max']))

print("\n" + "="*70)