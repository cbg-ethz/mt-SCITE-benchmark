rule generate_mutation_probability_matrices_for_fixed_error_rate:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/read_counts_{seed}.npy"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/mt_scite_mutation_prob_{inf_error_rate}_{seed}.txt"
    params:
        error_rate = "{inf_error_rate}"
        
    shell:"""python ../software/mt-SCITE/scripts/compute_mutation_probability.py {input} {output} {params.error_rate}"""


rule estimate_tree_given_mutation_probability_matrix:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/mt_scite_mutation_prob_{inf_error_rate}_{seed}.txt"
        
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/mt_scite_mutation_prob_{inf_error_rate}_{seed}_map0.newick"

    shell:"""./../software/mt-SCITE/mtscite -n {wildcards.num_mutations} -m {wildcards.num_cells} -r 1 -l 1000000 -fd 0.0001 -ad 0.0001 -s -max_treelist_size 1  -i {input}"""


##-------------------------------------------##
## Note that we are using two different ways of inferring mutation trees using mt-scite here.
## Above - check that mt-scite works for true error rate and is robust to the error rate fixed an order of magnitude higher and lower
## Below - check that mt-scite erstimates an error rate that is sufficiently close to the truth or at least does not bias tree reconstruction

rule generate_mutation_probability_matrices_for_range_of_error_rates:
    input:
        "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/read_counts_{seed}.npy"

    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/mutation_probs.done"
    params:
        input_dir = "../results/simulated_data/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/",
        output_dir = "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/"
        
        
    shell:
        """
        python ../src/learn_error_rates/compute_mutation_probs_for_error_learning.py {wildcards.error_rate} {params.input_dir} {params.output_dir} {wildcards.seed}
        touch ../results/inference_output/{wildcards.num_mutations}_{wildcards.concentration}_{wildcards.error_rate}_{wildcards.num_reads}_{wildcards.num_cells}/seed_{wildcards.seed}/mutation_probs.done
        """

rule perform_kcval_to_learn_error:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/mutation_probs.done"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores.txt"

    params:
        input_dir = "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}"

    shell:"""python ../software/mt-SCITE/scripts/cv_without_filtering.py --mtscite_bin_path ../software/mt-SCITE/mtscite --directory {params.input_dir} -o {params.input_dir}"""


rule pick_best_error_rate_given_only_tree_lik:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores.txt",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores_star_trees.txt",
        "../src/learn_error_rates/pick_best_error_rate_lik_only.py"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_learned_tree_1.gv",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_error_rate_1.csv"
    params:
        output_dir="../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/"

    shell:
        """
        python ../src/learn_error_rates/pick_best_error_rate_lik_only.py --tree_scores_file {input[0]} --star_tree_scores_file {input[1]} --output_dir {params.output_dir}
        """

        
rule pick_best_error_rate_and_tree:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores.txt",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores_star_trees.txt",
        "../src/learn_error_rates/pick_best_error_rate.py"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_learned_tree_2.gv",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_error_rate_2.csv"
    params:
        output_dir="../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/"

    shell:
        """
        python ../src/learn_error_rates/pick_best_error_rate.py --tree_scores_file {input[0]} --star_tree_scores_file {input[1]} --output_dir {params.output_dir}

        """

rule pick_best_error_rate_and_tree_diff:
    input:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores.txt",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores_star_trees.txt",
        "../src/learn_error_rates/pick_best_error_rate_diff.py"
    output:
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_learned_tree_3.gv",
        "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_error_rate_3.csv"
    params:
        output_dir="../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/"

    shell:
        """
        python ../src/learn_error_rates/pick_best_error_rate_diff.py --tree_scores_file {input[0]} --star_tree_scores_file {input[1]} --output_dir {params.output_dir}
        """


        
# This is to compare against just the division of logLik and log star
# rule pick_best_error_rate_and_tree_2:
#     input:
#         "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores.txt",
#         "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/val_scores_star_trees.txt",
#         "../src/learn_error_rates/pick_best_error_rate_2.py"
#     output:
#         "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_learned_tree_2.gv",
#         "../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/best_error_rate_2.csv"
#     params:
#         output_dir="../results/inference_output/{num_mutations}_{concentration}_{error_rate}_{num_reads}_{num_cells}_{initial_mutation_freq}/seed_{seed}/"

#     shell:
#         """
#         python ../src/learn_error_rates/pick_best_error_rate_2.py --tree_scores_file {input[0]} --star_tree_scores_file {input[1]} --output_dir {params.output_dir}

#         """




