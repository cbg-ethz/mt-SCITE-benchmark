rule run_scite:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/scite_input_{seed}.txt"
	
    resources:
        runtime=30
	
    output: "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/scite_{inf_error_rate}_{seed}_ml0.newick"
    
    shell:"""./../software/SCITE/scite -n {wildcards.num_mutations} -m {wildcards.num_cells} -r 1 -l 1000000 -fd 0.001 -ad 0.001 -max_treelist_size 1  -i {input} -o ../results/inference_output/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}_{wildcards.initial_mutation_freq}/scite_{wildcards.inf_error_rate}_{wildcards.seed}"""
