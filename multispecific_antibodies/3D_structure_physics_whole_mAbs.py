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
SYNC_MASTER = BASE_DIR / "mAb_Data_Sync_Master.csv"
OUTPUT_CSV = BASE_DIR / "mAb_3D_Physics_Features.csv"

"""
THEORY & BACKGROUND: 3D BIOPHYSICAL PROFILING

1. SOLVENT ACCESSIBLE SURFACE AREA (SASA) & THE SHRAKE-RUPLEY ALGORITHM
   - Theory: SASA measures the surface area of a biomolecule visible to a solvent. 
     A 'Probe Sphere' (radius ~1.4Å, approximating a water molecule) rolls over 
     the Van der Waals surface of the protein.
   - Algorithm: The Shrake-Rupley algorithm (1973) generates a mesh of points 
     on the surface of each atom. Tests determine if a point remains accessible 
     to the solvent probe or buried by a neighboring atom.
   - Significance: SASA accounts for the 3D 'Folding Context.' Buried hydrophobic 
     residues remain harmless, while surface-exposed residues create potential 
     'hotspots' for aggregation.

2. THE HYDROPHOBIC EFFECT & AGGREGATION PROPENSITY
   - Background: The 'Hydrophobic Effect' drives protein clumping. Water molecules 
     form highly ordered, low-entropy 'cages' around exposed non-polar side chains 
     (ALA, VAL, ILE, LEU, MET, PHE, TYR, TRP).
   - Hydro_Ratio: To minimize this energetic cost, hydrophobic patches stick 
     together to 'hide' from the water, causing irreversible aggregation. A high 
     Hydro_Ratio (Hydro_SASA / Total_SASA) flags 'sticky' antibodies likely to 
     fail in high-concentration liquid formulations.

3. NET CHARGE & COLLOIDAL STABILITY
   - Theory: Antibodies act as 'Colloidal' particles. Solubility depends on 
     electrostatic repulsion. According to DLVO theory, identical net charges 
     induce repulsion between molecules.
   - Isoelectric Point (pI) Context: At physiological pH (~7.4), ionizable groups 
     determine the net charge.
     - Positive: ARG, LYS, HIS (partially).
     - Negative: ASP, GLU.
   - Significance: A Net_Charge near zero (the Isoelectric Point) increases 
     aggregation risk because no electrostatic 'force field' prevents 
     molecular collisions.

4. MULTISPECIFIC ARCHITECTURE (CHAIN SPECIFICITY)
   - Background: For multispecific mAbs (4-chain H2L2 structures), a single 
     'weak' chain often bottlenecks stability. 
   - Strategy: Breaking down SASA and Charge by Chain ID (A, B, C, D) 
     identifies whether stability issues localize to a specific binding arm 
     or affect the entire scaffold.
"""

# --- BIOPHYSICAL SCALES ---
HYDROPHOBIC_RESIDUES = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
CHARGE_SCALE = {
    'ARG': 1.0, 'LYS': 1.0, 'HIS': 0.1,  # Positive
    'ASP': -1.0, 'GLU': -1.0,            # Negative
}

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
                    
                    for res in chain:
                        res_name = res.get_resname()
                        res_sasa = getattr(res, 'sasa', 0)
                        
                        c_sasa += res_sasa
                        if res_name in HYDROPHOBIC_RESIDUES:
                            c_hydro += res_sasa
                        c_charge += CHARGE_SCALE.get(res_name, 0)
                    
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
                "Hydro_SASA": round(hydro_sasa, 2),
                "Hydro_Ratio": round(hydro_sasa / total_sasa, 4) if total_sasa > 0 else 0,
                "Net_Charge": round(total_charge, 1)
            }
            entry.update(chain_stats) # Adds Ch_A, Ch_B, etc.
            master_features.append(entry)

        except Exception as e:
            print(f"Error processing {ab_name}: {e}")

    # 5. Save results
    if master_features:
        df = pd.DataFrame(master_features)
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"\n=== SUCCESS: Physics extracted for {len(df)} mAbs ===")
        print(f"File saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    run_3D_physics_analysis()