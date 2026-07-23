#!/bin/bash
ml HyperQueue

# CPU 
hq alloc add slurm \
    --time-limit 3h \
    --worker-start-cmd "ml palma/2024a uv" \
    -- \
        --partition=zen5 \
        --cpus-per-task=16 \
        --mem=32G \

# CPU 3
hq alloc add slurm \
    --time-limit 4h \
    --worker-start-cmd "ml palma/2024a uv" \
    -- \
        --partition=zen4 \
        --cpus-per-task=16 \
        --mem=32G






