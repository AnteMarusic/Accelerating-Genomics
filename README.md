## Goal
This project aims at accelerating with CUDA Smith Waterman and PairHMM algorithms. Smith Waterman has been implemented without the traceback part, while a complete version of PairHMM has been implemented. PairHMM has a bug in the C that I wasn't able to discover and that carried on to the CUDA version, however the acceleration of PairHMM is coherent with the C code as they outputs the same results on the test cases.

## Usage
I wrote some bash (.sh) files that can be found in the two folders. Those can be runned and will execute the implementation on the selected hardcoded input file.