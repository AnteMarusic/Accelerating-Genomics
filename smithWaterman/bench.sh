#!/bin/bash

PROGRAM1="hipvers"
INPUT_FILE="input.txt"
GENERATOR="generator"
FLAGS="-O3" 

# Compile the hipvers program
hipcc $FLAGS $PROGRAM1.cpp -o $PROGRAM1
if [ $? -eq 0 ]; then
    echo "Compiled $PROGRAM1.cpp successfully"
else
    echo "Failed to compile $PROGRAM1.cpp"
    exit 1
fi

# Loop to generate input files and run the hipvers program
for i in 64 128 256 512 1024; do
    echo "$GENERATOR.py $i $i"
    python3 $GENERATOR.py $i $i

    if [ $? -ne 0 ]; then
        echo "Failed to execute $GENERATOR.py with LEN=$LEN"
        exit 1
    fi

    for block_size in 32 64 128 256 512 1024; do
        # Define output file names
        OUTPUT1="${PROGRAM1}_${i}_${block_size}.txt"

        echo "Running: $PROGRAM1 with LEN=$i and BLOCK_SIZE=$block_size"
        echo "./$PROGRAM1 $INPUT_FILE $OUTPUT1 $block_size"
        ./$PROGRAM1 $INPUT_FILE $OUTPUT1 $block_size

        if [ $? -ne 0 ]; then
            echo "Failed to execute $PROGRAM1 with LEN=$i and BLOCK_SIZE=$block_size"
            exit 1
        fi
    done
done

echo "All tasks completed successfully."

   

