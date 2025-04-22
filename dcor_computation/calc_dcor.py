import numpy as np
import pandas as pd
import argparse

import dcor
from multiprocessing import Pool, cpu_count

# Load dataset
def load_data(file_path):
    print(f"Loading data from {file_path}")
    df = pd.read_csv(file_path, sep="\t")  # Assuming tab-separated file
    genes = df.iloc[:, 0].values  # Extract gene names
    data = df.iloc[:, 1:].values  # Extract expression values
    return genes, data

# Compute distance correlation for one pair
def compute_dcor_pair(args):
    """
    Helper function to compute distance correlation for a single pair of genes.
    """
    i, j, data = args
    if i == j:
        return (i, j, 1.0)  # Distance correlation of a variable with itself is 1
    else:
        dcor_value = dcor.distance_correlation(data[i, :], data[j, :])
        print(f"Done calculating for pair {i} {j}...")
        return (i, j, dcor_value)

def compute_dcor_matrix(data):
    """
    Compute an NxN symmetric matrix of pairwise distance correlation (dCor) values.

    Parameters:
    - data: numpy array of shape (N, M), where N = genes, M = samples.

    Returns:
    - dcor_matrix: NxN symmetric matrix of dCor values.
    """
    N, M = data.shape
    dcor_matrix = np.zeros((N, N))  # Initialize NxN matrix

    # Generate all pairs of indices (i, j) for the upper triangle
    pairs = [(i, j, data) for i in range(N) for j in range(i, N)]

    # Use multiprocessing to compute dCor in parallel
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(compute_dcor_pair, pairs)

    # Fill the dCor matrix with the results
    for i, j, value in results:
        dcor_matrix[i, j] = dcor_matrix[j, i] = value

    return dcor_matrix

# Save matrix to CSV
def save_matrix(genes, matrix, output_file):
    df = pd.DataFrame(matrix, index=genes, columns=genes)
    df.to_csv(output_file)


def main(args):
    input_file = args.input  # Change to your file path
    output_file = args.output

    print(f"Loading data from {input_file}")
    genes, data = load_data(input_file)
    print(f"Computing dcor matrix...")
    dcor_matrix = compute_dcor_matrix(data)
    print(f"Saving matrix to file {output_file}")
    save_matrix(genes, dcor_matrix, output_file)
    print(f"Distance correlation matrix saved to {output_file}")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculating dcor adjacency.")
    parser.add_argument("--input", type=str, required=True, help="Input txt.")
    parser.add_argument("--output", type=str, required=True, help="Output csv.")
    args = parser.parse_args()
    
    main(args)