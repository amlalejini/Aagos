# 2025-10-29 - Architecture enumeration experiment

Goal: enumerate all meaningfully different genetic architectures comprising four, four-bit genes.
For each architecture, evolve 10 independent populations for 10,000 generations in a static environment, disallowing architecture-altering mutations.

## Planning

This experiment is based on 2020-05-22--enum-architectures experiment.

Configuration

- NUM_GENES=4
- GENE_SIZE=4
- GENOME_SIZE=16

- Bit flip rates: 0.003, 0.03, 0.3
- Locked genetic architectures (x65)

- Hand-designed architectures
  - Enumerate all possible architectures (for small gene count, gene size, and genome size)
    - Characterize distribution of graph properties
    - Evolve each for 10k generations in random static environment (x10)
    - Look at landscape (graph properties vs. fitness)

NOTE: before running this experiment, be sure to unzip the appropriate genome architecture archive: `genome-architectures/num-genes_4__gene-size_4__genome-length_16`