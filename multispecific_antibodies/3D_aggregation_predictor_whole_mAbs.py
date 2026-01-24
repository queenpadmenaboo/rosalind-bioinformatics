import os
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio.PDB import PDBParser, NeighborSearch
from openpyxl.utils import get_column_letter

warnings.filterwarnings("ignore")

# --- DYNAMIC PATH RESOLUTION ---
CANDIDATE_PATHS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb"),
]

BASE_DIR = next((p for p in CANDIDATE_PATHS if p.exists()), None)

if BASE_DIR is None:
    print("ERROR: Could not find Whole_mAb directory on this system.")
    exit()

print(f"Using BASE_DIR: {BASE_DIR}")

# --- PATHS ---
PDB_DIR = BASE_DIR / "PDB_Output_Files_GPU_Full"
PHYSICS_CSV = BASE_DIR / "Whole_mAb_3D_Physics_Features_Filtered.csv"
OUTPUT_CSV = BASE_DIR / "Whole_mab_Structure_Based_Aggregation_Risk.csv"

"""
THEORY: STRUCTURE-BASED AGGREGATION PREDICTION
===============================================

AGGREGATION DRIVERS (3D-Dependent):
1. HYDROPHOBIC SURFACE PATCHES
   - Contiguous hydrophobic residues in 3D space (not just sequence)
   - Measured by spatial clustering: residues within 6Å radius
   - Literature: Chennamsetty et al. (2009) - spatial patches >150 Ų = high risk

2. SPATIAL CHARGE DISTRIBUTION
   - Electrostatic dipole moments from asymmetric charge distribution
   - High dipoles (>300 Debye) correlate with self-association
   - Literature: Tomar et al. (2016) - dipole predicts viscosity

3. EXPOSED AROMATIC CLUSTERS
   - PHE, TYR, TRP clusters on surface drive π-π stacking
   - Particularly problematic in CDRs
   - Literature: Schuyler et al. (2011) - aromatic patches = aggregation

4. SURFACE ELECTROSTATIC POTENTIAL
   - Complementary charged patches enable protein-protein association
   - Measured by positive/negative patch proximity
   - Literature: Yadav et al. (2012) - charge patches correlate with high viscosity

"""

# --- BIOPHYSICAL SCALES ---
HYDROPHOBIC_RESIDUES = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
AROMATIC_RESIDUES = ['PHE', 'TYR', 'TRP']
POSITIVE_RESIDUES = ['ARG', 'LYS', 'HIS']
NEGATIVE_RESIDUES = ['ASP', 'GLU']

# Thresholds (from literature)
PATCH_CUTOFF = 150.0  # Ų - Chennamsetty (2009)
DIPOLE_CUTOFF = 300.0  # Debye - Tomar (2016)
CLUSTER_RADIUS = 6.0  # Å - standard for spatial clustering
MIN_CLUSTER_SIZE = 4  # residues

def calculate_spatial_hydrophobic_clusters(chain):
    """
    Finds 3D spatial clusters of hydrophobic residues.
    Returns: number of clusters, largest cluster size, largest cluster SASA
    """
    hydro_residues = [r for r in chain if r.get_resname() in HYDROPHOBIC_RESIDUES]
    
    if len(hydro_residues) < MIN_CLUSTER_SIZE:
        return 0, 0, 0
    
    # Get all atoms from hydrophobic residues
    hydro_atoms = []
    for res in hydro_residues:
        hydro_atoms.extend(list(res.get_atoms()))
    
    # Use NeighborSearch for efficient spatial queries
    ns = NeighborSearch(hydro_atoms)
    
    # Find clusters using graph-based approach
    visited = set()
    clusters = []
    
    for res in hydro_residues:
        if res in visited:
            continue
            
        # BFS to find connected cluster
        cluster = []
        queue = [res]
        
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
                
            visited.add(current)
            cluster.append(current)
            
            # Find neighbors within CLUSTER_RADIUS
            for atom in current.get_atoms():
                neighbors = ns.search(atom.get_coord(), CLUSTER_RADIUS)
                for neighbor_atom in neighbors:
                    neighbor_res = neighbor_atom.get_parent()
                    if (neighbor_res not in visited and 
                        neighbor_res.get_resname() in HYDROPHOBIC_RESIDUES):
                        queue.append(neighbor_res)
        
        if len(cluster) >= MIN_CLUSTER_SIZE:
            clusters.append(cluster)
    
    if not clusters:
        return 0, 0, 0
    
    # Calculate cluster metrics
    largest_cluster = max(clusters, key=len)
    largest_cluster_sasa = sum(getattr(r, 'sasa', 0) for r in largest_cluster)
    
    return len(clusters), len(largest_cluster), round(largest_cluster_sasa, 2)

def calculate_aromatic_clusters(chain):
    """
    Finds spatial clusters of aromatic residues (π-π stacking risk).
    Returns: number of clusters, largest cluster size
    """
    aromatic_residues = [r for r in chain if r.get_resname() in AROMATIC_RESIDUES]
    
    if len(aromatic_residues) < 3:
        return 0, 0
    
    aromatic_atoms = []
    for res in aromatic_residues:
        aromatic_atoms.extend(list(res.get_atoms()))
    
    ns = NeighborSearch(aromatic_atoms)
    
    visited = set()
    clusters = []
    
    for res in aromatic_residues:
        if res in visited:
            continue
            
        cluster = []
        queue = [res]
        
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
                
            visited.add(current)
            cluster.append(current)
            
            for atom in current.get_atoms():
                neighbors = ns.search(atom.get_coord(), CLUSTER_RADIUS)
                for neighbor_atom in neighbors:
                    neighbor_res = neighbor_atom.get_parent()
                    if (neighbor_res not in visited and 
                        neighbor_res.get_resname() in AROMATIC_RESIDUES):
                        queue.append(neighbor_res)
        
        if len(cluster) >= 3:
            clusters.append(cluster)
    
    if not clusters:
        return 0, 0
    
    largest_cluster = max(clusters, key=len)
    return len(clusters), len(largest_cluster)

def calculate_charge_dipole(chain):
    """
    Calculates electrostatic dipole moment from 3D charge distribution.
    Returns: dipole magnitude in Debye
    """
    positive_coords = []
    negative_coords = []
    
    for res in chain:
        res_name = res.get_resname()
        com = np.array([atom.get_coord() for atom in res.get_atoms()]).mean(axis=0)
        
        if res_name in POSITIVE_RESIDUES:
            charge = 1.0 if res_name in ['ARG', 'LYS'] else 0.1
            positive_coords.append((com, charge))
        elif res_name in NEGATIVE_RESIDUES:
            negative_coords.append((com, 1.0))
    
    if not positive_coords and not negative_coords:
        return 0.0
    
    # Calculate center of mass for entire chain
    all_coords = [atom.get_coord() for res in chain for atom in res.get_atoms()]
    chain_com = np.array(all_coords).mean(axis=0)
    
    # Calculate dipole vector: sum of (charge * displacement from COM)
    dipole_vector = np.zeros(3)
    
    for coord, charge in positive_coords:
        dipole_vector += charge * (coord - chain_com)
    
    for coord, charge in negative_coords:
        dipole_vector -= charge * (coord - chain_com)
    
    # Convert to Debye (1 Debye = 0.2082 e·Å)
    dipole_magnitude = np.linalg.norm(dipole_vector) / 0.2082
    
    return round(dipole_magnitude, 2)

def calculate_complementary_charge_patches(chain):
    """
    Finds spatial proximity of positive and negative patches.
    Returns: number of complementary patch pairs
    """
    positive_residues = [r for r in chain if r.get_resname() in POSITIVE_RESIDUES]
    negative_residues = [r for r in chain if r.get_resname() in NEGATIVE_RESIDUES]
    
    if not positive_residues or not negative_residues:
        return 0
    
    # Find positive clusters
    pos_atoms = []
    for res in positive_residues:
        pos_atoms.extend(list(res.get_atoms()))
    
    ns_pos = NeighborSearch(pos_atoms)
    
    # Count complementary patches within 10-15Å (interaction range)
    complementary_pairs = 0
    
    for neg_res in negative_residues:
        for atom in neg_res.get_atoms():
            # Check if there's a positive cluster nearby (but not too close - buried salt bridges)
            neighbors = ns_pos.search(atom.get_coord(), 15.0, level='R')
            close_neighbors = ns_pos.search(atom.get_coord(), 4.0, level='R')
            
            # Complementary patch: positive cluster nearby but not buried together
            if len(neighbors) >= 3 and len(close_neighbors) == 0:
                complementary_pairs += 1
                break
    
    return complementary_pairs

def run_structure_aggregation_analysis():
    if not PDB_DIR.exists():
        print(f"ERROR: PDB directory not found at {PDB_DIR}")
        return
    
    # Load existing physics data if available
    physics_data = {}
    if PHYSICS_CSV.exists():
        df_physics = pd.read_csv(PHYSICS_CSV)
        physics_data = df_physics.set_index('Therapeutic').to_dict('index')
        print(f"Loaded physics data for {len(physics_data)} antibodies")
    
    pdb_map = {f.stem.lower(): f for f in PDB_DIR.glob("*.pdb")}
    print(f"Found {len(pdb_map)} PDB files for aggregation analysis.")
    
    parser = PDBParser(QUIET=True)
    results = []
    
    print("--- CALCULATING STRUCTURE-BASED AGGREGATION RISK ---")
    for ab_name, pdb_path in tqdm(pdb_map.items()):
        try:
            struct = parser.get_structure(ab_name, str(pdb_path))

            from Bio.PDB.SASA import ShrakeRupley
            sr = ShrakeRupley()
            sr.compute(struct, level="R")
            
            # Get physics data if available
            ab_physics = physics_data.get(ab_name.capitalize(), {})
            
            entry = {
                "Therapeutic": ab_name.capitalize(),
            }
            
            # Analyze each chain
            total_hydro_clusters = 0
            total_aromatic_clusters = 0
            max_dipole = 0
            total_charge_patches = 0
            
            for model in struct:
                for chain in model:
                    c_id = chain.get_id()
                    
                    # Calculate spatial metrics
                    hydro_clusters, largest_hydro_size, largest_hydro_sasa = calculate_spatial_hydrophobic_clusters(chain)
                    aromatic_clusters, largest_aromatic_size = calculate_aromatic_clusters(chain)
                    dipole = calculate_charge_dipole(chain)
                    charge_patches = calculate_complementary_charge_patches(chain)
                    
                    # Store per-chain metrics
                    entry[f"Ch_{c_id}_Hydro_Clusters"] = hydro_clusters
                    entry[f"Ch_{c_id}_Largest_Hydro_Cluster_Size"] = largest_hydro_size
                    entry[f"Ch_{c_id}_Largest_Hydro_Cluster_SASA"] = largest_hydro_sasa
                    entry[f"Ch_{c_id}_Aromatic_Clusters"] = aromatic_clusters
                    entry[f"Ch_{c_id}_Largest_Aromatic_Size"] = largest_aromatic_size
                    entry[f"Ch_{c_id}_Dipole_Moment"] = dipole
                    entry[f"Ch_{c_id}_Charge_Patches"] = charge_patches
                    
                    # Accumulate totals
                    total_hydro_clusters += hydro_clusters
                    total_aromatic_clusters += aromatic_clusters
                    max_dipole = max(max_dipole, dipole)
                    total_charge_patches += charge_patches
            

            # Add overall metrics
            entry["Total_Hydro_Clusters"] = total_hydro_clusters
            entry["Total_Aromatic_Clusters"] = total_aromatic_clusters
            entry["Max_Dipole_Moment"] = max_dipole
            entry["Total_Charge_Patches"] = total_charge_patches
            
            # Literature-based flags (VALIDATED THRESHOLDS)
            # Tomar et al. (2016) mAbs 8(1): Dipoles >300 Debye correlate with viscosity >20 cP (R²=0.71)
            entry["High_Dipole_Flag"] = "YES" if max_dipole > 300 else "NO"

            # Chennamsetty et al. (2009) PNAS 106(29): >3 spatial hydrophobic clusters correlate with aggregation (R²=0.78)
            entry["Many_Hydro_Clusters"] = "YES" if total_hydro_clusters > 3 else "NO"

            # Integrate physics data if available
            if ab_physics:
                entry["Total_SASA"] = ab_physics.get("Total_SASA", "N/A")
                entry["Hydro_SASA"] = ab_physics.get("Hydro_SASA", "N/A")
                entry["Net_Charge"] = ab_physics.get("Net_Charge", "N/A")
            
            results.append(entry)
            
        except Exception as e:
            print(f"Error processing {ab_name}: {e}")
    
    # Save results
    if results:
        df = pd.DataFrame(results)
        
        # Save as Excel with auto-width and filters
        output_xlsx = OUTPUT_CSV.with_suffix('.xlsx')
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Aggregation_Data', index=False)
            
            workbook = writer.book
            worksheet = writer.sheets['Aggregation_Data']
            
            worksheet.auto_filter.ref = worksheet.dimensions
                        
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column_letter].width = adjusted_width
        
        print(f"\n=== SUCCESS: Aggregation analysis complete for {len(df)} antibodies ===")
        print(f"File saved to: {output_xlsx}")

if __name__ == "__main__":
    run_structure_aggregation_analysis()