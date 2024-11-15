rule generate_mutation_probability_matrices_for_fixed_error_rate:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/read_counts_{seed}.npy"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/mt_scite_mutation_prob_{inf_error_rate}_{seed}.txt"
    params:
        error_rate = "{inf_error_rate}"
        
    shell:"""python ../software/mt-SCITE/scripts/compute_mutation_probability.py {input} {output} {params.error_rate}"""


rule estimate_tree_given_mutation_probability_matrix:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/mt_scite_mutation_prob_{inf_error_rate}_{seed}.txt"
        
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/mt_scite_mutation_prob_{inf_error_rate}_{seed}_map0.newick"

    shell:"""./../software/mt-SCITE/mtscite -n {wildcards.num_mutations} -m {wildcards.num_cells} -r 1 -l 1000000 -fd 0.0001 -ad 0.0001 -s -max_treelist_size 1  -i {input}"""


##-------------------------------------------##
## Note that we are using two different ways of inferring mutation trees using mt-scite here.
## Above - check that mt-scite works for true error rate and is robust to the error rate fixed an order of magnitude higher and lower
## Below - check that mt-scite erstimates an error rate that is sufficiently close to the truth or at least does not bias tree reconstruction

rule generate_mutation_probability_matrices_for_range_of_error_rates:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/read_counts_{seed}.npy"

    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/seed_{seed}/mutation_probs.done"
    params:
        input_dir = "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/",
        output_dir = "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/"
        
        
    shell:
        """
        python ../src/learn_error_rates/compute_mutation_probs_for_error_learning.py {wildcards.error_rate} {params.input_dir} {params.output_dir} {wildcards.seed}
        touch ../results/inference_output/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/seed_{wildcards.seed}/mutation_probs.done
        """

rule perform_kcval_to_learn_error:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/seed_{seed}/mutation_probs.done"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}/seed_{seed}/val_scores.txt"

    shell:"""python ../software/mt-SCITE/scripts/cv_without_filtering.py --mtscite_bin_path ../software/mt-SCITE/ --directory {input} -o {input}"""
        
