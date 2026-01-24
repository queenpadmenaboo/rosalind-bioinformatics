import os
import warnings
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings("ignore")

# --- PATHS (Targeting your verified folder) ---
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")
PDB_DIR = BASE_DIR / "PDB_Output_Files_GPU_Full"
SYNC_MASTER = BASE_DIR / "Whole_mAb_Data_Sync_Master.csv"
OUTPUT_CSV = BASE_DIR / "Whole_mAb_3D_Physics_Features_Filtered.csv"

"""
THEORY & BACKGROUND: 3D BIOPHYSICAL PROFILING (VALIDATED METRICS)

1. SOLVENT ACCESSIBLE SURFACE AREA (SASA) & THE SHRAKE-RUPLEY ALGORITHM
   - Theory: SASA measures the surface area of a biomolecule accessible to solvent. 
     A probe sphere (radius ~1.4Å, approximating water) rolls over the Van der Waals 
     surface of the protein.
   - Algorithm: Shrake-Rupley (1973) generates a mesh of points on each atom's surface 
     and tests if points remain solvent-accessible or are buried by neighbors.
   - Significance: SASA distinguishes buried vs. surface-exposed residues in the 3D 
     folded context. Surface-exposed hydrophobic residues create aggregation hotspots.

2. HYDROPHOBIC SURFACE PATCHES & AGGREGATION PROPENSITY
   - Theory: Aggregation is driven by localized hydrophobic patches, not bulk 
     hydrophobicity. Contiguous surface regions of hydrophobic residues (ALA, VAL, 
     ILE, LEU, MET, PHE, TYR, TRP) create "sticky" zones that promote intermolecular 
     association.
   - Metrics Calculated:
     * Fv_Hydro_Fraction: Hydrophobic SASA / Total SASA in variable domains (residues 1-120 
       heavy, 1-110 light) - captures antibody-specific differences
     * Fc_Hydro_Fraction: Hydrophobic SASA / Total SASA in constant domains (residues 120+ 
       heavy, 110+ light) - captures isotype-specific differences (G1 vs G2 vs G4)
     * Max_Hydrophobic_Patch: Largest contiguous hydrophobic surface area (5-residue 
       sliding window) - correlates with experimental aggregation propensity
   - Literature Basis: Sharma et al. (2014), Jain et al. (2017) - patch size matters 
     more than total surface hydrophobicity for aggregation

3. NET CHARGE & CHARGE SYMMETRY
   - Theory: At pH 7.4, net charge determines electrostatic repulsion between antibody 
     molecules (DLVO theory). Charge asymmetry between chains creates dipole moments 
     that can promote self-association.
   - Metrics Calculated:
     * Net_Charge: Total ionizable residues (ARG, LYS, HIS positive; ASP, GLU negative)
     * Heavy_Charge_Asymmetry: |Charge(Chain_A) - Charge(Chain_C)| - flags unbalanced 
       heavy chains in multispecifics
     * Light_Charge_Asymmetry: |Charge(Chain_B) - Charge(Chain_D)| - flags unbalanced 
       light chains
   - Literature Basis: Therapeutic Antibody Profiler (TAP) guidelines - charge asymmetry 
     thresholds correlate with viscosity and clearance

4. REGIONAL ANALYSIS FOR MULTISPECIFIC ANTIBODIES
   - Rationale: In multispecific H2L2 structures, each arm can have distinct properties. 
     A single high-risk arm can limit the entire molecule's developability.
   - Strategy: Calculate Fv vs. Fc metrics separately, and per-chain breakdowns (A, B, 
     C, D), to localize stability bottlenecks to specific binding domains or the Fc 
     scaffold.
"""

# --- BIOPHYSICAL SCALES ---
HYDROPHOBIC_RESIDUES = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
CHARGE_SCALE = {
    'ARG': 1.0, 'LYS': 1.0, 'HIS': 0.1,  # Positive
    'ASP': -1.0, 'GLU': -1.0,            # Negative
}

# Heavy chain: Fv = 1-120, Fc = 120+
# Light chain: Fv = 1-110, CL = 110+
FV_CUTOFF_HEAVY = 120
FV_CUTOFF_LIGHT = 110
PATCH_WINDOW = 5  # Sliding window for hydrophobic patches

def run_3D_physics_analysis():
    if not PDB_DIR.exists():
        print(f"ERROR: PDB directory not found at {PDB_DIR}")
        return

    # 1. Map available PDBs (lowercase for easy matching)
    pdb_map = {f.stem.lower(): f for f in PDB_DIR.glob("*.pdb")}
    print(f"Found {len(pdb_map)} PDB files ready for analysis.")

    # 2. Setup tools
    parser = PDBParser(QUIET=True)
    sr = ShrakeRupley()
    master_features = []

    # 3. Process every PDB found
    print(f"--- EXTRACTING PHYSICS DATA ---")
    for ab_name, pdb_path in tqdm(pdb_map.items()):
        try:
            struct = parser.get_structure(ab_name, str(pdb_path))
            sr.compute(struct, level="R") # Residue-level SASA
            
            total_sasa, hydro_sasa, total_charge = 0, 0, 0
            chain_stats = {}

            # Analyze every chain (A, B, C, D)
            for model in struct:
                for chain in model:
                    c_id = chain.get_id()
                    c_sasa, c_hydro, c_charge = 0, 0, 0
                    
                    residues = list(chain.get_residues())
                    is_heavy = c_id in ['A', 'C']
                    cutoff = FV_CUTOFF_HEAVY if is_heavy else FV_CUTOFF_LIGHT

                    fv_sasa, fv_hydro, fc_sasa, fc_hydro = 0, 0, 0, 0
                    max_patch = 0

                    for i, res in enumerate(residues):
                        res_name = res.get_resname()
                        res_num = res.get_id()[1]
                        res_sasa = getattr(res, 'sasa', 0)
                        
                        c_sasa += res_sasa
                        c_charge += CHARGE_SCALE.get(res_name, 0)
                        
                        is_hydro = res_name in HYDROPHOBIC_RESIDUES
                        
                        if is_hydro:
                            c_hydro += res_sasa

                        # Regional breakdown
                        if res_num < cutoff:  # Fv region
                            fv_sasa += res_sasa
                            if is_hydro:
                                fv_hydro += res_sasa
                        else:  # Fc/CL region
                            fc_sasa += res_sasa
                            if is_hydro:
                                fc_hydro += res_sasa

                    # Hydrophobic patch detection
                    for i in range(len(residues) - PATCH_WINDOW + 1):
                        window = residues[i:i+PATCH_WINDOW]
                        if all(r.get_resname() in HYDROPHOBIC_RESIDUES for r in window):
                            patch_sasa = sum(getattr(r, 'sasa', 0) for r in window)
                            max_patch = max(max_patch, patch_sasa)

                    chain_stats[f"Ch_{c_id}_Fv_Hydro"] = round(fv_hydro / fv_sasa, 4) if fv_sasa > 0 else 0
                    chain_stats[f"Ch_{c_id}_Fc_Hydro"] = round(fc_hydro / fc_sasa, 4) if fc_sasa > 0 else 0
                    chain_stats[f"Ch_{c_id}_Max_Patch"] = round(max_patch, 2)
                    
                    
                    # Store individual chain data
                    chain_stats[f"Ch_{c_id}_SASA"] = round(c_sasa, 2)
                    chain_stats[f"Ch_{c_id}_Charge"] = round(c_charge, 1)
                    
                    total_sasa += c_sasa
                    hydro_sasa += c_hydro
                    total_charge += c_charge

            # 4. Compile the data for this mAb
            entry = {
                "Therapeutic": ab_name.capitalize(),
                "Total_SASA": round(total_sasa, 2),
                "Net_Charge": round(total_charge, 1),
                "Hydro_SASA": round(hydro_sasa, 2),    
            }
            entry.update(chain_stats) # Adds Ch_A, Ch_B, etc.
            master_features.append(entry)

        except Exception as e:
            print(f"Error processing {ab_name}: {e}")

    # 5. Save results
    if master_features:
        df = pd.DataFrame(master_features)

        # Merge with master data for isotypes
        if SYNC_MASTER.exists():
            df_master = pd.read_csv(SYNC_MASTER)
            df['Therapeutic'] = df['Therapeutic'].str.strip().str.lower()
            df_master['Therapeutic'] = df_master['Therapeutic'].str.strip().str.lower()
            df = pd.merge(df, df_master[['Therapeutic', 'CH1 Isotype', 'VD LC']], on='Therapeutic', how='left')
            cols = ['Therapeutic', 'CH1 Isotype', 'VD LC'] + [c for c in df.columns if c not in ['Therapeutic', 'CH1 Isotype', 'VD LC']]
            df = df[cols]
        
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"\n=== SUCCESS: Physics extracted for {len(df)} mAbs ===")
        print(f"File saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    run_3D_physics_analysis()