#!/bin/bash
cd ..
make clean
make Aagos
./Aagos -POP_SIZE 1 -NUM_BITS 10 -NUM_GENES 5 -GENE_SIZE 3 -MAX_GENS 10 -PRINT_INTERVAL 1 -STATISTICS_INTERVAL 1 -SNAPSHOT_INTERVAL 5 -SEED 6 > ./Tests/correct.txt