import argparse
import numpy as np
import pandas as pd

import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

from scipy.stats import pearsonr
from sklearn.decomposition import PCA

# Load data
def load_data(moduleA_file, moduleB_file, exprA_file, exprB_file, adjA_file, adjB_file):
    # Load module memberships
    moduleA = pd.read_csv(moduleA_file)
    moduleB = pd.read_csv(moduleB_file)
    
    # Load expression data
    exprA = pd.read_csv(exprA_file, index_col=0, sep="\t")
    exprB = pd.read_csv(exprB_file, index_col=0, sep="\t")
    
    # Load adjacency matrices
    adjA = np.loadtxt(adjA_file)
    adjB = np.loadtxt(adjB_file)
    
    return moduleA, moduleB, exprA, exprB, adjA, adjB

# Create a NetworkX graph
def get_top_central_genes(adj_matrix, genes_in_module, top_n=10):
    G = nx.Graph(adj_matrix)
    mapping = {i: genes_in_module[i] for i in range(len(genes_in_module))}
    G = nx.relabel_nodes(G, mapping)

    # Calculate eigenvector centrality
    eigenvector_centrality = nx.eigenvector_centrality(G)
    sorted_genes = sorted(eigenvector_centrality.items(), key=lambda x: x[1], reverse=True)[:top_n]
    return [gene for gene, _ in sorted_genes[:top_n]]


# Perform PCA on each module
def compute_pc1(data):
    pca = PCA(n_components=1)  # Extract only PC1
    pc1_scores = pca.fit_transform(data.T)  # Transpose to have samples as rows
    return pc1_scores.flatten()  # Convert to 1D array


# Calculate correlation between two modules
def calculate_correlation(exprA, exprB, genenames_A, genenames_B):
    """
    Calculate Pearson correlation between two modules.
    Filters for common genes between the two modules and calculates the correlation.
    If no common genes exist, return corrval = None, pval = 1.
    """
    exprA_top10 = exprA.loc[exprA.index[exprA.index.isin(genenames_A)]]
    exprB_top10 = exprB.loc[exprB.index[exprB.index.isin(genenames_B)]]
    pc1_A = compute_pc1(exprA_top10.transpose())
    pc1_B = compute_pc1(exprB_top10.transpose())
    assert(len(pc1_A)==len(pc1_B))
    
    # Calculate Pearson correlation
    corrval, pval = pearsonr(np.abs(pc1_A), np.abs(pc1_B))

    return corrval, pval

def extract_label(filename):
        # Define a mapping of terms
        label_map = {
            "nonsepsis": "nonsepsis",
            "preshock": "nonshock",  # Convert preshock to nonshock
            "sepsis": "sepsis",
            "shock": "shock"
        }
        """Extract the condition label from a filename and map it."""
        for key in label_map.keys():
            if key in filename:
                return label_map[key]
        return None  # Return None if no match found

# Plot and save correlation matrix as PNG
def plot_correlation_matrix(correlation_matrix, pval_matrix, output_png, xlabel="State B", ylabel="State A"):
    """
    Plot the correlation matrix as a heatmap with color intensity based on p-values.
    """
        # Function to map p-values to significance annotations
    def annotate_significance(pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return 'n.s'

    # Create annotation matrix
    annotation_matrix = pval_matrix.applymap(annotate_significance)

    # Plot the heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        correlation_matrix,
        vmin=-1, vmax=1,  # Set color range for correlation
        cmap="RdBu_r",  # Custom Red-Blue colormap 
        annot=annotation_matrix.values,  # Annotate with significance
        fmt='',  # No formatting for annotations
        cbar_kws={'label': 'Correlation'},  # Add colorbar label
        linewidths=0.5  # Add lines between cells
    )

    # Remove x and y tick labels
    # plt.xticks([])
    # plt.yticks([])

    # Improve label formatting
    plt.xticks(rotation=20, ha="right", fontsize=10)
    plt.yticks(fontsize=10)
    plt.title(f"Correlation Matrix: {ylabel} vs {xlabel}", fontsize=14)
    # plt.xlabel(xlabel, fontsize=12)
    # plt.ylabel(ylabel, fontsize=12)

    # Add color blocks for index and columns
    for i, color in enumerate(correlation_matrix.index):
        plt.gca().add_patch(plt.Rectangle((-0.5, i), 0.2, 1, color=color, transform=plt.gca().get_yaxis_transform(), clip_on=False))

    for j, color in enumerate(correlation_matrix.columns):
        plt.gca().add_patch(plt.Rectangle((j, -0.5), 1, 0.2, color=color, transform=plt.gca().get_xaxis_transform(), clip_on=False))

    # Adjust layout for better spacing
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_png, dpi=300, bbox_inches='tight')  # Save with high resolution and tight b
    print(f"Heatmap saved to {output_png}")


# Main function
def main(args):
    """
    Main function to calculate the correlation matrix between modules in SetA and SetB.
    """
    # Load arguments
    moduleA_file = args.moduleA  # Path to module membership file for Set A
    moduleB_file = args.moduleB  # Path to module membership file for Set B
    expression_file_A = args.exprA  # Path to expression data for Set A
    expression_file_B = args.exprB  # Path to expression data for Set B
    adjA_file = args.adjA  # Path to adjacency matrix for Set A
    adjB_file = args.adjB  # Path to adjacency matrix for Set B
    output_file = args.output  # Output file path

    # Load data
    print("Loading data....")
    moduleAlist, moduleBlist, exprA, exprB, adjA, adjB = load_data(
        moduleA_file, moduleB_file, expression_file_A, expression_file_B, adjA_file, adjB_file
    )

    # Get unique modules
    unique_modules_A = moduleAlist.iloc[:, 1].unique()  # Modules in Set A
    unique_modules_B = moduleBlist.iloc[:, 1].unique()  # Modules in Set B

    print("Modules in Set A:", unique_modules_A)
    print("Modules in Set B:", unique_modules_B)

    # Convert results to matrices
    correlation_matrix = pd.DataFrame(index=unique_modules_A, columns=unique_modules_B)
    pvalue_matrix = pd.DataFrame(index=unique_modules_A, columns=unique_modules_B)

    for moduleA in unique_modules_A:
        for moduleB in unique_modules_B:
            print(f"Processing {moduleA} (Set A) vs {moduleB} (Set B)")

            # Genes in Set A module
            genes_in_A = moduleAlist[moduleAlist.iloc[:, 1] == moduleA].iloc[:, 0].values

            # Calculate top 10 influential genes in Set A
            genes_in_A_index = [i for i, gene in enumerate(moduleAlist.iloc[:, 0]) if gene in genes_in_A]
            sub_adj_matrixA = adjA[np.ix_(genes_in_A_index, genes_in_A_index)]
            genenames_A = get_top_central_genes(sub_adj_matrixA, genes_in_A, top_n=10)
            print(f"Top 10 genes in {moduleA} (Set A):", genenames_A)

            # Genes in Set B module
            genes_in_B = moduleBlist[moduleBlist.iloc[:, 1] == moduleB].iloc[:, 0].values

            # Calculate top 10 influential genes in Set B
            genes_in_B_index = [i for i, gene in enumerate(moduleBlist.iloc[:, 0]) if gene in genes_in_B]
            sub_adj_matrixB = adjB[np.ix_(genes_in_B_index, genes_in_B_index)]
            genenames_B = get_top_central_genes(sub_adj_matrixB, genes_in_B, top_n=10)
            print(f"Top 10 genes in {moduleB} (Set B):", genenames_B)

            # Filter expression data for the top 10 genes
            exprA_top10 = exprA.loc[exprA.index[exprA.index.isin(genenames_A)]]
            exprB_top10 = exprB.loc[exprB.index[exprB.index.isin(genenames_B)]]

            # # Get the minimum length between the two datasets
            x = min(len(exprA_top10), len(exprB_top10))

            # Ensure both datasets have the same length (top x genes)
            exprA_top10 = exprA_top10.iloc[:x]
            exprB_top10 = exprB_top10.iloc[:x]

            # Compute PC1 scores
            pc1_A = compute_pc1(exprA_top10.transpose())
            pc1_B = compute_pc1(exprB_top10.transpose())

            # Calculate Pearson correlation
            r, pval = pearsonr(np.abs(pc1_A), np.abs(pc1_B))
            print(f"Correlation: {r}, P-value: {pval}")

            # Plot if significant
            if pval < 0.05:
                plt.figure(figsize=(6, 6))
                plt.scatter(np.abs(pc1_A), np.abs(pc1_B), color='teal', alpha=0.7, edgecolor='k')
                plt.xlabel(f'PC1 ({moduleA})', fontsize=15)
                plt.ylabel(f'PC1 ({moduleB})', fontsize=15)
                plt.title(f'{moduleA} vs. {moduleB}\nr = {r:.3f}, p = {pval:.3g}', fontsize=12)

                # Add a regression line
                m, b = np.polyfit(np.abs(pc1_A), np.abs(pc1_B), 1)
                plt.plot(np.abs(pc1_A), m * np.abs(pc1_A) + b, color='red', linestyle='--', label='Fit')
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.5)
                plt.tight_layout()
                plt.savefig(f'{output_file}-{moduleA}_vs_{moduleB}.png', dpi=300)  # Save as PNG
            else:
                print("No significant correlation (p >= 0.05), plot not generated.")

            # Append data to matrices
            correlation_matrix.loc[moduleA, moduleB] = r
            pvalue_matrix.loc[moduleA, moduleB] = pval

    # Save results
    correlation_matrix.to_csv(f"{output_file}-corr.csv")
    print(f"Correlation matrix saved to {output_file}-corr.csv")

    pvalue_matrix.to_csv(f"{output_file}-pval.csv")
    print(f"P-value matrix saved to {output_file}-pval.csv")

    # Plot and save the correlation matrix
    # xlabel = extract_label(expression_file_A)
    # ylabel = extract_label(expression_file_B)
    # plot_correlation_matrix(correlation_matrix, pvalue_matrix, f"{output_file}-serial.png", xlabel, ylabel)

    
# Run the program   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculating correlation matrix between modules")
    parser.add_argument("--moduleA", type=str, required=True, help="Path to the CSV file containing module membership for Set A")
    parser.add_argument("--moduleB", type=str, required=True, help="Path to the CSV file containing module membership for Set B")
    parser.add_argument("--exprA", type=str, required=True, help="Path to the CSV file containing the expression dataset for Set A")
    parser.add_argument("--exprB", type=str, required=True, help="Path to the CSV file containing the expression dataset for Set B")
    parser.add_argument("--adjA", type=str, required=True, help="Path to the CSV/TXT file of adjacency matrix for Set A")
    parser.add_argument("--adjB", type=str, required=True, help="Path to the CSV/TXT file of adjacency matrix for Set B")
    parser.add_argument("--output", type=str, required=True, help="Name of the file containing the correlation values")

    args = parser.parse_args()
    
    # Run the main function
    main(args)
