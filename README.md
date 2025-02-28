# mt-scite-benchmark

Simulate benchmark data and perform inference using mt-SCITE, SCITE and Merlin using the snakemake workflow under workflow.

The simulation code can be found in ./src/simulation.

The functions that choose an error rate given the results of the k-fold cross validation procedure are available under ./src/learn_error_rates/

Under ./software, I stored the working versions of mt-SCITE, SCITE and Merlin.

To produce the figures of the method comparison, use the script ./src/method_comparison/marimo-plot-method-comparion.py