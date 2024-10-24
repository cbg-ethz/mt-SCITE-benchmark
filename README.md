# mt-scite-benchmark

Simulate benchmark data.

Decide on params:
1) num_reads = 500, based on dataset yf_2001, see script src/exploratory/determine_num_reads.R

2) error_rate; in Joannas paper: 8e-4 to 6e-2, thus set to: [5e-4, 5e-3, 5e-2]