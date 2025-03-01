import marimo

__generated_with = "0.8.3"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import os
    import argparse
    import shutil
    return argparse, mo, np, os, pd, shutil


@app.cell
def __():
    def compute_normalized_likelihood(dat, dat_star, epsilon=0.0001):
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

        # skip first row because that has the error rates and not the likelihood values
        for row in range(1, normalised_dat.shape[0]):
            normalised_dat.iloc[row, 1:] = dat.iloc[row, 1:] / (dat.iloc[row, 1:] - dat_star.iloc[row, 1:] + epsilon)
            
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
    return compute_normalized_likelihood, format_best_error_rate


@app.cell
def __():
    return


@app.cell
def __():
    # Input parameters
    n_cells = 500
    n_mutations = 10
    true_error_rate = 0.0005
    seed = 3
    initial_mut_freq = 0.01
    epsilon = 0.0001
    return (
        epsilon,
        initial_mut_freq,
        n_cells,
        n_mutations,
        seed,
        true_error_rate,
    )


@app.cell
def __(
    compute_normalized_likelihood,
    epsilon,
    initial_mut_freq,
    n_cells,
    n_mutations,
    np,
    os,
    pd,
    seed,
    shutil,
    true_error_rate,
):
    # File paths
    file_name = os.path.join(
            f"../../results/inference_output/"
            f"{n_mutations}_120_{true_error_rate}_500_{n_cells}_{initial_mut_freq}/seed_{seed}/val_scores.txt"
        )
    file_name_2 = os.path.join(
            f"../../results/inference_output/"
            f"{n_mutations}_120_{true_error_rate}_500_{n_cells}_{initial_mut_freq}/seed_{seed}/val_scores_star_trees.txt"
        )

    print(file_name)
    # Load data
    dat = pd.read_csv(file_name, header=None)
    dat_star = pd.read_csv(file_name_2, header=None)

    # Compute normalized likelihood
    normalised_dat = compute_normalized_likelihood(dat, dat_star, epsilon)

    # Calculate error rate, mean log likelihood, and standard deviation
    error_rate = normalised_dat.iloc[0, 1:].values
    # take the mean and standard deviation for each error rate, excluding the error rate itself (row 1) and the index row (col 1)
    mean_log_lik = np.mean(normalised_dat.iloc[1:normalised_dat.shape[0], 1:], axis=0)
    sd_log_lik = np.std(normalised_dat.iloc[1:normalised_dat.shape[0], 1:], axis=0)

    # Find the best error rate
    max_log_lik_index = np.argmax(mean_log_lik)
    max_log_lik = mean_log_lik.iloc[max_log_lik_index]
    best_error_rate = error_rate[max_log_lik_index]
    sd_max_lik = sd_log_lik.iloc[max_log_lik_index]


    # Output the result
    best_tree_file = f"../../results/inference_output/{n_mutations}_120_{true_error_rate}_500_{n_cells}_{initial_mut_freq}/seed_{seed}/learned_{best_error_rate:6f}_0_map0.gv"

    new_tree_file =f"../../results/inference_output/{n_mutations}_120_{true_error_rate}_500_{n_cells}_{initial_mut_freq}/seed_{seed}/best_learned_tree_{best_error_rate:6f}_0_map0.gv"

    # Copy the best tree file to the new file
    if os.path.exists(best_tree_file):
            shutil.copy(best_tree_file, new_tree_file)
            print(f"Copied best tree file to {new_tree_file}")
    else:
            print(f"Best tree file not found: {best_tree_file}")
    return (
        best_error_rate,
        best_tree_file,
        dat,
        dat_star,
        error_rate,
        file_name,
        file_name_2,
        max_log_lik,
        max_log_lik_index,
        mean_log_lik,
        new_tree_file,
        normalised_dat,
        sd_log_lik,
        sd_max_lik,
    )


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
