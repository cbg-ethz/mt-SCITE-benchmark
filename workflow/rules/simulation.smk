
rule simulate_data:
    output:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/tree_{seed}.gml",
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/read_counts_{seed}.npy"

    params:
        num_mutations=lambda wildcards: wildcards.num_mutations,
        concentration=lambda wildcards: wildcards.concentration,
        error_rate=lambda wildcards: wildcards.error_rate,
        num_reads=lambda wildcards: wildcards.num_reads,
        num_cells=lambda wildcards: wildcards.num_cells,
        seed=lambda wildcards: wildcards.seed

    shell: """python ../src/simulation/simulate_data.py {params.num_mutations} {params.concentration} {params.error_rate} {params.num_reads} {params.num_cells} ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells} {params.seed}"""


rule adapt_data_for_scite:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/read_counts_{seed}.npy"
    output:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/scite_input_{seed}.txt"

    shell:"""python ../src/simulation/adapt_data_for_scite.py ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/read_counts.npy ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/ {wildcards.seed}"""

