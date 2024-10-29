rule run_scite:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/scite_input_{seed}.txt"
    output: "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/scite_{seed}_ml0.newick"
    shell:"""./../software/SCITE/scite -n {wildcards.num_mutations} -m {wildcards.num_cells} -r 1 -l 1000000 -fd {wildcards.error_rate} -ad 0.000001 -max_treelist_size 1  -i ../results/simulated_data/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/scite_input.txt -o ../results/inference_output/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/scite_{wildcards.seed}"""
