import numpy as np
from Bio import PDB
import os
import sys
import argparse
from Bio.PDB import Superimposer
from scipy.spatial import cKDTree, distance
from scipy.optimize import linear_sum_assignment
from scipy.stats import wasserstein_distance

def extract_ca_coordinates(file_path):
    """Extract Cα atom coordinates from PDB or CIF files"""
    _, file_extension = os.path.splitext(file_path)
    if file_extension.lower() == '.cif':
        parser = PDB.MMCIFParser()
        structure = parser.get_structure('structure', file_path)
        ca_coords = [atom.get_coord() for atom in structure.get_atoms() if atom.get_name() == 'CA']
        return np.array(ca_coords)
    elif file_extension.lower() == '.pdb':
        parser = PDB.PDBParser(QUIET=True, PERMISSIVE=True)
        ca_coords = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        ca_coords.append([x, y, z])
                    except ValueError:
                        continue
        return np.array(ca_coords)
    else:
        raise ValueError("Unsupported file format: {}".format(file_extension))

def calculate_metrics(ref_ca, pred_ca, thresholds=[3.0, 2.0, 1.0]):
    """Calculate evaluation metrics at multiple thresholds"""
    results = {}
    
    pos_tree = cKDTree(pred_ca)
    dist_ref_to_pred, idx_ref_to_pred = pos_tree.query(ref_ca, k=1)
    
    ref_tree = cKDTree(ref_ca)
    dist_pred_to_ref, idx_pred_to_ref = ref_tree.query(pred_ca, k=1)
    
    for threshold in thresholds:
        mask_true = dist_ref_to_pred < threshold
        TP_true = np.sum(mask_true)
        
        mask_pred = dist_pred_to_ref < threshold
        TP_pred = np.sum(mask_pred)
        
        recall = TP_true / len(ref_ca) if len(ref_ca) > 0 else 0.0
        precision = TP_pred / len(pred_ca) if len(pred_ca) > 0 else 0.0

        if precision + recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0

        if TP_true > 0:
            matched_ref = ref_ca[mask_true]
            matched_pred = pred_ca[idx_ref_to_pred[mask_true]]
            squared_diff = np.sum(np.square(matched_pred - matched_ref), axis=1)
            rmsd = np.sqrt(np.mean(squared_diff))
        else:
            rmsd = float('nan')

        results[threshold] = {
            'TP_true': TP_true, 
            'TP_pred': TP_pred, 
            'recall': recall,
            'precision': precision,
            'f1': f1,
            'rmsd': rmsd
        }
    
    return results

def chamfer_distance(true_coords, pred_coords):
    """Calculate Chamfer Distance"""
    tree_pred = cKDTree(pred_coords)
    dist1, _ = tree_pred.query(true_coords)
    term1 = np.mean(dist1)
    
    tree_true = cKDTree(true_coords)
    dist2, _ = tree_true.query(pred_coords)
    term2 = np.mean(dist2)
    
    return term1 + term2


def earth_mover_distance(true_coords, pred_coords):
    """Calculate Earth Mover's Distance"""
    n_true = len(true_coords)
    n_pred = len(pred_coords)
    
    if n_true == n_pred:
        dist_matrix = distance.cdist(true_coords, pred_coords, 'euclidean')
        
        row_ind, col_ind = linear_sum_assignment(dist_matrix)
        
        emd = dist_matrix[row_ind, col_ind].sum() / n_true
        return emd
    
   
    emd_x = wasserstein_distance(true_coords[:, 0], pred_coords[:, 0])
    emd_y = wasserstein_distance(true_coords[:, 1], pred_coords[:, 1])
    emd_z = wasserstein_distance(true_coords[:, 2], pred_coords[:, 2])
    
    return (emd_x + emd_y + emd_z) / 3

def main(true_file, pred_file):
    true_coords = extract_ca_coordinates(true_file)
    pred_coords = extract_ca_coordinates(pred_file)

    n_true = len(true_coords)
    n_pred = len(pred_coords)

    metrics_results = calculate_metrics(true_coords, pred_coords, thresholds=[3.0, 2.0, 1.0])
    
    cd = chamfer_distance(true_coords, pred_coords)
    emd = earth_mover_distance(true_coords, pred_coords)
    
    results = {
        'n_true': n_true,
        'n_pred': n_pred,
        'chamfer_distance': cd,
        'earth_mover_distance': emd
    }
    
    for threshold in [3.0, 2.0, 1.0]:
        prefix = f"{int(threshold)}A"
        threshold_results = metrics_results[threshold]
        results.update({
            f"{prefix}_correct": threshold_results['TP_true'],
            f"{prefix}_precision": threshold_results['precision'],
            f"{prefix}_recall": threshold_results['recall'],
            f"{prefix}_f1": threshold_results['f1'],
            f"{prefix}_rmsd": threshold_results['rmsd']
        })
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate evaluation metrics from PDB files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
            python CAEVAL.py -r native.pdb -p predicted.pdb -o eval_result.log
            python CAEVAL.py --reference-structure native.cif --predicted-structure predicted.cif --output-file eval_result.log
        """
    )

    parser.add_argument(
        "--reference-structure",
        "--reference",
        "-r",
        type=str,
        required=True,
        help="Path to reference/native structure file (PDB or CIF format)"
    )
    parser.add_argument(
        "--predicted-structure",
        "--predicted",
        "-p",
        type=str,
        required=True,
        help="Path to predicted structure file (PDB or CIF format)"
    )

    parser.add_argument(
        "--output-file",
        "-o",
        type=str,
        default="evaluation_results.json",
        help="Path to output result file (default: evaluation_results.json)"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    metrics = main(args.reference_structure, args.predicted_structure)

    output_fields = [
        metrics['n_true'], metrics['n_pred'],
        metrics['chamfer_distance'], metrics['earth_mover_distance']
    ]
    
    for threshold in [3.0, 2.0, 1.0]:
        prefix = f"{int(threshold)}A"
        output_fields.extend([
            metrics[f"{prefix}_correct"],
            metrics[f"{prefix}_precision"],
            metrics[f"{prefix}_recall"],
            metrics[f"{prefix}_f1"],
            metrics[f"{prefix}_rmsd"]
        ])
    
    output_line = "\t".join(str(x) for x in output_fields) + "\n"
    
    if not os.path.exists(args.output_file):
        header_fields = [
            "n_true", "n_pred",
            "chamfer_distance", "earth_mover_distance"
        ]
        
        for threshold in [3.0, 2.0, 1.0]:
            prefix = f"{int(threshold)}A"
            header_fields.extend([
                f"{prefix}_correct",
                f"{prefix}_precision",
                f"{prefix}_recall",
                f"{prefix}_f1",
                f"{prefix}_rmsd"
            ])
        
        header_line = "\t".join(header_fields) + "\n"
        with open(args.output_file, 'w') as f:
            f.write(header_line)
    
    with open(args.output_file, 'a') as f:
        f.write(output_line)
