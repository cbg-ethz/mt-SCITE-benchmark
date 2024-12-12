import marimo as mo
import pandas as pd
import numpy as np
import os
import argparse
import shutil

def compute_normalized_likelihood(dat, dat_star, epsilon):
    """
    Compute the normalized likelihood matrix from the input data and star data.

    Parameters:
    dat (pd.DataFrame): Tree likelihood on heldout data
    dat_star (pd.DataFrame): Star tree likelihood on heldout data
    epsilon (float): Small value to avoid division by zero.

    Returns:
    pd.DataFrame: Normalized likelihood matrix.
    """
    normalised_dat = dat.copy()
    normalised_dat.iloc[1, :] = dat.iloc[1, :] / (dat.iloc[1, :] - dat_star.iloc[1, :] + epsilon)
    normalised_dat.iloc[2, :] = dat.iloc[2, :] / (dat.iloc[2, :] - dat_star.iloc[2, :] + epsilon)
    normalised_dat.iloc[3, :] = dat.iloc[3, :] / (dat.iloc[3, :] - dat_star.iloc[3, :] + epsilon)
    return normalised_dat

def format_best_error_rate(error_rate):
    """
    Format the best error rate to avoid scientific notation and trailing zeros.

    Parameters:
    error_rate (float): The error rate to format.

    Returns:
    str: Formatted error rate.
    """
    return f"{error_rate:.10f}".rstrip('0').rstrip('.')  # Remove trailing zeros and decimal point if unnecessary


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Find the best error rate and tree file.")
    parser.add_argument("--tree_scores_file", type=str, required=True, help="Path to the file, where the scores for every tree are saved.")
    parser.add_argument("--star_tree_scores_file", type=str, required=True, help="Path to the file, where the scores on the star tree are saved.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory to save the best tree file.")
    parser.add_argument("--epsilon", type=float, default=0.0001, help="Epsilon value to avoid division by zero.")
    args = parser.parse_args()

    # Load data
    dat = pd.read_csv(args.tree_scores_file, header=None)
    dat_star = pd.read_csv(args.star_tree_scores_file, header=None)

    # Compute normalized likelihood
    normalised_dat = compute_normalized_likelihood(dat, dat_star, args.epsilon)

    # Calculate error rate, mean log likelihood, and standard deviation
    error_rate = normalised_dat.iloc[0, 1:].values
    mean_log_lik = np.mean(normalised_dat.iloc[1:31, 1:], axis=0)
    sd_log_lik = np.std(normalised_dat.iloc[1:31, 1:], axis=0)

    # Find the best error rate
    max_log_lik_index = np.argmax(mean_log_lik)
    max_log_lik = mean_log_lik.iloc[max_log_lik_index]
    best_error_rate = error_rate[max_log_lik_index]
    sd_max_lik = sd_log_lik.iloc[max_log_lik_index]

    # Create data frame for best error rate
    data = {
        "best_error_rate": [best_error_rate],
        "mean_likelihood": [max_log_lik],
        "sd_likelihood": [sd_max_lik]
    }
    df = pd.DataFrame(data)

    # Determine the best tree file path
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    best_tree_file = os.path.join(output_dir, f"learned_{best_error_rate:6f}_0_map0.gv")
    new_tree_file = os.path.join(output_dir, f"best_learned_tree.gv")
    new_error_rate_file = os.path.join(output_dir, "best_error_rate.csv")
    
    # Copy the best tree file to the new file
    if os.path.exists(best_tree_file):
        shutil.copy(best_tree_file, new_tree_file)
        print(f"Copied best tree file to {new_tree_file}")

        df.to_csv(new_error_rate_file, index=False)
        print(f"Written best error rate and likelihood to {new_error_rate_file}")
        
    else:
        print(f"Best tree file not found: {best_tree_file}")

if __name__ == "__main__":
    main()
