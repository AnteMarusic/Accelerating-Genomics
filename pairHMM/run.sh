#!/bin/bash
#PROGRAM1="antidiagsPairHMM"
PROGRAM1="pairHMMmatrix"
#PROGRAM1="testantid"
#INPUT_FILE="test_set/test.in"
#OUTPUT_FILE="test_set/test.out"
INPUT_FILE="test_set/10s.in"
OUTPUT_FILE="test_set/10s.out"
FLAGS= "-Wall -Wextra -Werror -std=c99 -pedantic -O3 -lm"
MATH="-lm"

gcc $FLAGS $PROGRAM1.c -o $PROGRAM1.exe $MATH
if [ $? -eq 0 ]; then
    echo "Compiled $PROGRAM1.c successfully"
else
    echo "Failed to compile $PROGRAM1.c"
    exit 1
fi

echo "running: $PROGRAM1"
./$PROGRAM1".exe" $INPUT_FILE $OUTPUT_FILE

