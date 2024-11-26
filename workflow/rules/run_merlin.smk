rule run_merlin:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/merlin_variant_matrix_{seed}.csv"
    output: "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/merlin_{seed}_clone_tree_edge_list.txt"

    resources:
        mem_mb=4000,
        runtime=120
    params:
        variant_matrix="../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/merlin_variant_matrix_{seed}.csv",
        total_matrix="../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/merlin_total_matrix_{seed}.csv"

    #conda:
    #    "merlin.yml"
            
    shell:"""python ./../software/MERLIN/src/merlin.py -v {params.variant_matrix} -t {params.total_matrix}  -o ../results/inference_output/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}/merlin_{wildcards.seed}"""
