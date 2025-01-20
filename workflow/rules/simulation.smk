
rule simulate_data:
    output:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/tree_{seed}.gml",
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/read_counts_{seed}.npy",
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/mutation_prob_{seed}.npy"

    params:
        num_mutations=lambda wildcards: wildcards.num_mutations,
        concentration=lambda wildcards: wildcards.concentration,
        error_rate=lambda wildcards: wildcards.error_rate,
        num_reads=lambda wildcards: wildcards.num_reads,
        num_cells=lambda wildcards: wildcards.num_cells,
        seed=lambda wildcards: wildcards.seed,
        initial_mutation_freq=lambda wildcards: wildcards.initial_mutation_freq

    #conda:
    #    "merlin.yml"

    shell: """python ../src/simulation/simulate_data.py {params.num_mutations} {params.concentration} {params.error_rate} {params.num_reads} {params.num_cells} ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq} {params.seed} {params.initial_mutation_freq}"""


rule adapt_data_for_scite:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/read_counts_{seed}.npy"
    output:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/scite_input_{seed}.txt"

    shell:"""python ../src/simulation/adapt_data_for_scite.py {input}  ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}/ {wildcards.seed}"""


rule adapt_data_for_merlin:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/read_counts_{seed}.npy"
    output:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/merlin_variant_matrix_{seed}.csv"

    #conda:
    #    "merlin.yml"

    shell:"""python ../src/simulation/adapt_data_for_merlin.py {input} ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}/ {wildcards.seed}"""


rule simulate_data_with_erroneous_mutations:
    output:
        "../results/simulated_data_with_err_muts/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}_{n_err_mutations}_{n_err_samples}/tree_{seed}.gml",
        "../results/simulated_data_with_err_muts/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}_{n_err_mutations}_{n_err_samples}/read_counts_{seed}.npy",
        "../results/simulated_data_with_err_muts/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}_{n_err_mutations}_{n_err_samples}/mutation_prob_{seed}.npy"

    params:
        num_mutations=lambda wildcards: wildcards.num_mutations,
        concentration=lambda wildcards: wildcards.concentration,
        error_rate=lambda wildcards: wildcards.error_rate,
        num_reads=lambda wildcards: wildcards.num_reads,
        num_cells=lambda wildcards: wildcards.num_cells,
        seed=lambda wildcards: wildcards.seed,
        initial_mutation_freq=lambda wildcards: wildcards.initial_mutation_freq,
        n_err_mutations=lambda wildcards: wildcards.n_err_mutations,
        n_err_samples=lambda wildcards: wildcards.n_err_samples

    #conda:
    #    "merlin.yml"

    shell:""" echo ../results/simulated_data_with_err_muts/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}_{wildcards.n_err_mutations}_{wildcards.n_err_samples}
    python ../src/simulation/simulate_data.py {params.num_mutations} {params.concentration} {params.error_rate} {params.num_reads} {params.num_cells} ../results/simulated_data_with_err_muts/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}_{wildcards.n_err_mutations}_{wildcards.n_err_samples} {params.seed} {params.initial_mutation_freq} {params.n_err_mutations} {params.n_err_samples}
    """
