#!/bin/bash
ml HyperQueue

# CPU 1
hq alloc add slurm \
    --time-limit 1h \
    --worker-start-cmd "ml palma/2024a" \
    --name cpu-pool \
    -- \
        --partition=zen5 \
        --cpus-per-task=16 \
        --mem=16G

# CPU 2
hq alloc add slurm \
    --time-limit 3h \
    --worker-start-cmd "ml palma/2024a" \
    --name cpu-pool \
    -- \
        --partition=zen4 \
        --cpus-per-task=16 \
        --mem=20G

# CPU 3
hq alloc add slurm \
    --time-limit 5h \
    --worker-start-cmd "ml palma/2024a" \
    --name cpu-pool \
    -- \
        --partition=normal \
        --cpus-per-task=16 \
        --mem=32G

