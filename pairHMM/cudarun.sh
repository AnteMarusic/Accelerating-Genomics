#!/bin/bash

PROGRAM1="pairHMM"
INPUT_FILE="test_set/10s.in"
OUTPUT_FILE="test_set/10s.out"

# Compile the cuda program
nvcc $FLAGS $PROGRAM1.cu -o $PROGRAM1
if [ $? -eq 0 ]; then
    echo "Compiled $PROGRAM1.cu successfully"
else
    echo "Failed to compile $PROGRAM1.cu"
    exit 1
fi

./$PROGRAM1 $INPUT_FILE $OUTPUT_FILE
if [ $? -ne 0 ]; then
    echo "Failed to execute $PROGRAM1"
    exit 1
fi

echo "All tasks completed successfully."
