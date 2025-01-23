import os
import sys
import numpy as np
import subprocess

def generate_error_rates(true_error_rate):
    """Generate 20 error rates spaced logarithmically between min and max."""
    error_rate_min = true_error_rate / 10
    error_rate_max = true_error_rate * 10
    error_rates = np.round(
        np.logspace(np.log10(error_rate_min), np.log10(error_rate_max), num=20),
        decimals=6
    )
    return error_rates

def prepare_directories(input_dir, output_dir, seed_nr):
    """Prepare input and output directories for a specific seed."""
    input_file = os.path.join(input_dir, f"read_counts_{seed_nr}.npy")
    output_dir_seed = os.path.join(output_dir, f"seed_{seed_nr}")
    
    if not os.path.exists(output_dir_seed):
        os.makedirs(output_dir_seed)
    
    return input_file, output_dir_seed

def compute_mutation_probabilities(input_file, output_dir_seed, error_rates):
    """Run the mutation probability computation script for each error rate."""
    for error_rate in error_rates:
        output_file = os.path.join(output_dir_seed, f"{error_rate:.6f}.csv")
        command = [
            "python", "../src/simulation/adapt_data_for_mtscite.py",
            input_file, output_file, str(error_rate)
        ]
        try:
            subprocess.run(command, check=True)
            print(f"Successfully computed mutation probability matrix for error rate {error_rate}")
            
        except subprocess.CalledProcessError as e:
            print(f"Error running script for error rate {error_rate}: {e}")


    # create file indicating completion to snakemake
    file_path = os.path.join(output_dir_seed, "mutation_probs.done")
    with open(file_path, "w") as file:
        pass
    print(f"File created at: {file_path}")


def main():
    if len(sys.argv) != 5:
        print("Usage: python compute_mutation_probs_for_error_learning.py <true_error_rate> <input_dir> <output_dir> <seed>")
        sys.exit(1)
    
    # Parse command-line arguments
    true_error_rate = float(sys.argv[1])
    input_dir = sys.argv[2]
    output_dir = sys.argv[3]
    seed_nr = sys.argv[4]
    
    # Generate error rates
    error_rates = generate_error_rates(true_error_rate)
    
    # Process seed
    print(f"Processing seed {seed_nr}...")
    input_file, output_dir_seed = prepare_directories(input_dir, output_dir, seed_nr)
    compute_mutation_probabilities(input_file, output_dir_seed, error_rates)
    

if __name__ == "__main__":
    main()

