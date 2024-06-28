#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <hip/hip_runtime.h>

//DEBUGGING MACROS______________________________________________________________
// Comment out the following line to disable the DEBUG macro
#undef DEBUG
//#define DEBUG

#ifdef DEBUG
#define PRINT printf
#else
#define PRINT // macros
#endif

#define CHECK(call)                                                            \
{                                                                              \
    const hipError_t error = call;                                            \
    if (error != hipSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                hipGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

//used to find the execution time of the program
inline double seconds(){
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

//CONSTANTS_____________________________________________________________________
#define MAX_LINE_LENGTH 10000
#define SCORE_MATCH 1
#define SCORE_MISMATCH (-1)
#define SCORE_OPEN_GAP (-3)
#define SCORE_EXTEND_GAP (-1)
#define P 1
#define Q 2
#define D 3


/*
The convention used in the matrixes is the following:
the X is the horizontal axis and the Y is the vertical axis,
the origin is in the top left corner of the matrix.
so the indexes are ix for the columns and iy for the rows.
Example of how the sequences are in the matrixes
the number of columns is always less than the number of rows
s1 = CG
s2 = CCGA
this two sequences are in the matrixes in this way
    _CG
_
C   000
C   000
G   000
A   000

ny is the number of rows and nx is the number of columns
int this example
nx is the size of s1 + 1
ny is the size of s2 + 1
*/

//FUNCTIONS_____________________________________________________________________
//used to print the three antidiagonals. Used for debugging
__device__ void a_print(int *a, int len, int nx) {
    int k = 1;
    int q = 1;
    PRINT("\n____________________________\n");
    for(int i = 0; i < len; i ++) {
        a[i] == INT_MIN ? printf("-∞,") : printf("%d,", a[i]);
        if (i + 1 == (nx) * k) {
            if (i + 1 == (nx * 3) * q) {
                PRINT("\n");
                q++;
            }
            PRINT("|");
            k ++;
        }
    }
    PRINT("\n____________________________\n");
    return;
}

//used to print the array. Used for debugging
__device__ void array_print(int *a, int len) {
    PRINT("\n____________________________\n");
    for(int i = 0; i<len; i++) {
        PRINT("%d ", a[i]);
    }
    PRINT("\n____________________________\n");
}

//functions to sum and find the max of two integers handling correctly the case of negative infinity (INT_MIN)
__device__ int max_of_two_integers(int a, int b) {
    return (a > b) ? a : b;
}

__device__ int sum_with_infinity (int a, int b) {
    return (a == INT_MIN || b == INT_MIN) ? INT_MIN : a + b;
}

__device__ int max_of_four_integers (int a, int b, int c, int d) {
    return max_of_two_integers (max_of_two_integers (a, b), max_of_two_integers(c, d));
}

//used to reverse a sequence
void reverse_array(char arr[], int size) {
  // Loop through half the array
  for (int i = 0; i < size / 2; i++) {
    // Swap elements at opposite positions
    char temp = arr[i];
    arr[i] = arr[size - 1 - i];
    arr[size - 1 - i] = temp;
  }
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
//Only three antidiagonals at a time are saved. 
//This function is used to get the value of a coordinate that is currently saved.
__device__ int m_get(int *val, int iy, int ix, int mat, int *three_antidiagonals, int nx, int ny) {
    int antid_dim = nx;
    int offset_mat = 0;
    int index = 0;
    int cur_antid_num = ix + iy;
    int offset_antid = cur_antid_num % 3;  //three is the num of antid saved at the same time
    //check the input
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        //printf("[m_get thread %d] out of bounds set iy: %d, ix: %d\n",threadIdx.x, iy, ix);
        return -1;
    }
    //compute the index of the required value
    switch (mat) {
        case P :
            offset_mat = 0;
            break;
        case Q :
            offset_mat = antid_dim*3;
            break;
        case D :
            offset_mat = antid_dim*2*3;
            break;
        default :
            //printf("[m_get thread %d] mat doesn't exist\n", threadIdx.x);
    }
    
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }
    //get the value
    *val = three_antidiagonals[index];
    //PRINT("[m_get thread %d] getting val: %d, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", threadIdx.x, *val, iy, ix, mat, index);
    return 1;
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
//Only three antidiagonals at a time are saved. 
//This function is used to set the value of a coordinate that is currently saved.
__device__ int m_set(int *val, int iy, int ix, int mat, int *three_antidiagonals, int nx, int ny) {
    int antid_dim = nx;
    int offset_mat = 0;
    int index = 0;
    int cur_antid_num = ix + iy;
    int offset_antid = cur_antid_num % 3;  //three is the num of antid saved at the same time

    //check the input
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        //printf("[m_set thread %d] out of bounds set iy: %d, ix: %d\n", threadIdx.x, iy, ix);
        return -1;
    }
    //compute the index of the required value
    switch (mat) {
        case P :
            offset_mat = 0;
            break;
        case Q :
            offset_mat = antid_dim*3;
            break;
        case D :
            offset_mat = antid_dim*2*3;
            break;
        default :
            //printf("[m_set thread %d] mat doesn't exist\n", threadIdx.x);
    }
    
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }

    //PRINT("[m_set thread %d] setting val: %d, iy:%d, ix:%d, at index: %d on matrix %d\n", threadIdx.x, *val, iy, ix, index, mat);
    //set the value
    three_antidiagonals[index] = *val;
    return 1;
}

//KERNEL_______________________________________________________________________
__global__ void alignGPU(char **sequences, size_t *sequences_len, int sequences_num, int *d_result, int *d_error, int three_antidiagonals_size, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
    int bid = blockIdx.x + blockIdx.y + blockIdx.z; //block index
    int tid = threadIdx.x + threadIdx.y + threadIdx.z; //thread index
    int grid_size = gridDim.x*gridDim.y*gridDim.z;
    int num_of_threads = blockDim.x*blockDim.y*blockDim.z;
    char *d_sx = sequences[bid*2];
    char *d_sy = sequences[bid*2+1];
    //sequences are accessed based on the bid
    int sx_len = sequences_len[bid*2];
    int sy_len =sequences_len[bid*2+1];
    int nx = sx_len+1;
    int ny = sy_len+1;
    extern __shared__ int s[]; //declare the whole shared memory
    int *three_antidiagonals = s; //part of shared memory to store the three antidiagonals
    int *shared_max = &s[three_antidiagonals_size]; //part of shared memory to store the max of each thread
    int antid_num = nx + ny - 1;
    int yr_phase = nx - 1;//yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //with each iteration. red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    int o_phase = ny - nx + 1;//orange phase is the one where the cur_antid_dim is the max possible
    int err;
    int val = 1;
    int minus_infty = INT_MIN;
    int zero=0;
    int iy = 0, ix = 0;
    int cur_antid_dim;
    int gen_iy = 0; //general indexes to keep track of the starting point of the antidiagonal on the y axis
    int gen_ix = 0; //general indexes to keep track of the starting point of the antidiagonal on the x axis
    int temp1,temp2,temp3;

    //safety check to avoid accessing unwanted memory areas if the gridsize has not been correctly set.
    if (bid >= sequences_num) {*d_error = -1; return;} 

    //debugging
    if(bid == 0 && tid == 0) {
        //PRINT("bid: %d, tid: %d\n", bid, tid);
        array_print(d_result, grid_size);
    }
    if(tid==0){
        //PRINT("[bid: %d, tid: %d] sequences[%d] = %s\nsequences[%d] = %s\n", bid, tid, bid*2, sequences[bid*2], bid*2+1, sequences[bid*2+1]);
    }

    //initialize the maximmum value of each thread
    shared_max[threadIdx.x] = 0;
    for(int cur_antid_num = 0; cur_antid_num < antid_num; cur_antid_num  ++) {
        //compute dimension of the current antidiagonal
        if(cur_antid_num <= yr_phase) { //at index yr_phase starts the o_phase
            //PRINT("[kernel %d] growing phase\n", tid);
            cur_antid_dim = cur_antid_num + 1;
        }
        else if(cur_antid_num >= yr_phase + o_phase) {
            //PRINT("[kernel %d] decreasing phase\n", tid);
            cur_antid_dim = cur_antid_dim - 1;
        }
        else {
            //PRINT("[kernel] const phase\n");
        }
        //compute the start of ix for the current antidiagonal
        if(cur_antid_num >= yr_phase + o_phase) {
            gen_ix ++;
        }
        iy = gen_iy - tid;
        ix = gen_ix + tid;
        while(true) {
            //check if the indexes are out of bounds
            if (ix < 0 || iy < 0 || ix >= nx || iy >= ny) {
                //PRINT("out of bounds iy = %d, ix = %d\n", iy, ix);
                break;
            }
            //set the first row
            if (iy == 0) {
                err = m_set(&minus_infty, iy, ix, P, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, Q, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
            }
            //set the first col
            else if(ix==0){
                err = m_set(&zero, iy, ix, P, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&minus_infty, iy, ix, Q, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
            }
            //dynamic programming
            else {
                err = m_get(&temp1, iy-1, ix, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy-1, ix, P, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                int b = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&b, iy, ix, P, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}

                err = m_get(&temp1, iy, ix-1, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy, ix-1, Q, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                int a = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&a, iy, ix, Q, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}

                err = m_get(&temp1, iy, ix, P, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy, ix, Q, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp3, iy-1, ix-1, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                //accessing the y string is reversed to decrease the number of accesses to global memory
                temp3 = d_sy[sy_len - iy] == d_sx[ix - 1] ? temp3 + SCORE_MATCH : temp3 + SCORE_MISMATCH;
                val = max_of_four_integers(temp1, temp2, temp3, 0);
                err = m_set(&val, iy, ix, D, three_antidiagonals, nx, ny);
                if(err<1) {*d_error = -1; return;}
                //update maximum
                shared_max[tid] = val > shared_max[tid] ? val : shared_max[tid];
            }
            //update the matrix indexes
            iy -= num_of_threads;
            ix += num_of_threads;
            __syncthreads();
        }
        //move the starting point in case we are in the growing part of the matrix
        if (gen_iy < ny - 1) {
            gen_iy ++;
        }
    }

    //COMPUTE MAX between threads
    for (int s = num_of_threads / 2; s > 0; s >>= 1) {
    if (tid < s && tid + s < num_of_threads) {
        shared_max[tid] = fmaxf(shared_max[tid], shared_max[tid + s]);
    }
    __syncthreads();
    }
    // The first thread in the block writes the result to global memory
    if (tid == 0) {
        atomicMax(&d_result[bid], shared_max[0]);
    }
}


//the main function reads the input file and then calls the kernel
//no error handling is present in case the strings are longer than the maximum allowed length
//no error handling in case the structure of the file is wrong.
//the expected structure is the first line is the number of sequences and then the sequences are in pairs
//and are as many as the number of sequences specified in the first line.
int main(int argc, char *argv[]) {
    char line1[MAX_LINE_LENGTH];
    char line2[MAX_LINE_LENGTH];
    int longest_antidiagonal_len = 0; //len of the longest antidiagonal
    int num_of_sequences; //number of lines in the file
    int sx_len, sy_len;
    char **d_sequences; //pointer to the sequences in the device
    char **h_sequences; //pointer to the sequences in the host
    int *d_result; //pointer to the result in the device
    int *h_result;
    int result_len;
    int *d_error; //used to check if the kernel has been executed correctly MAYBE HAS TO BE AN ARRAY
    int *h_error;
    size_t *h_sequences_len;
    size_t *d_sequences_len;
    int three_antidiagonals_size;
    int nBytes;
    double iStart, iElaps;
    char *output_file_path;
    int block_size;
    FILE *input_file;
    FILE *output_file;



     // set up device
    int dev = 1;
    hipDeviceProp_t deviceProp;
    CHECK(hipGetDeviceProperties(&deviceProp, dev));
    printf("[main] Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(hipSetDevice(dev));
    //check the main arguments
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_file_path> <output_file_path> <block_size>\n", argv[0]);
        return 1;
    }
    output_file_path = argv[2];
    block_size = atoi(argv[3]);
    //open the file
    input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    if (fgets(line1, sizeof(line1), input_file) == NULL) {
        printf("file is empty");
        return 1; // Exit loop if EOF or error
    }
    //get the number of sequences in the file (stored as the first line of the file)
    num_of_sequences = atoi(line1);
    printf("num_of_sequences: %d\n", num_of_sequences);
    //allocate the memory for the sequences, the result and the error
    h_sequences = (char**)malloc(num_of_sequences * sizeof(char*));
    h_sequences_len = (size_t *)malloc(num_of_sequences * sizeof(size_t));
    result_len = num_of_sequences/2; //number of results (since we have two sequences for each alignment is the half of the sequences)
    h_result = (int *)malloc(result_len*sizeof(int));
    h_error = (int *)malloc(sizeof(int));
    
    //read the sequences from the file and store them in the device memory
    for (int i = 0; i < num_of_sequences; i+=2) {
        PRINT("%d\n", i);
        // Read the first line
        if (fgets(line1, sizeof(line1), input_file) == NULL) {
            break; // Exit loop if EOF or error
        }
        // Read the second line
        if (fgets(line2, sizeof(line2), input_file) == NULL) {
            break;
        }
        //PRINT("[host] line1: %s\n line2: %s\n", line1, line2);
        longest_antidiagonal_len = strlen(line1) > longest_antidiagonal_len ? strlen(line1) : longest_antidiagonal_len;
        longest_antidiagonal_len = strlen(line2) > longest_antidiagonal_len ? strlen(line2) : longest_antidiagonal_len;
        if(strlen(line1) > strlen(line2)){
            sx_len = strlen(line2);
            sy_len = strlen(line1);
            h_sequences_len[i] = sx_len;
            h_sequences_len[i+1] = sy_len;
            reverse_array(line1, sy_len);
            CHECK(hipMalloc(&h_sequences[i], sx_len * sizeof(char)));
            CHECK(hipMalloc(&h_sequences[i+1], sy_len * sizeof(char)));
            CHECK(hipMemcpy(h_sequences[i], line2, sx_len * sizeof(char), hipMemcpyHostToDevice));
            CHECK(hipMemcpy(h_sequences[i+1], line1, sy_len * sizeof(char), hipMemcpyHostToDevice));
        }
        else{
            sx_len = strlen(line1);
            sy_len = strlen(line2);
            h_sequences_len[i] = sx_len;
            h_sequences_len[i+1] = sy_len;
            reverse_array(line2, sy_len);
            CHECK(hipMalloc(&h_sequences[i], sx_len * sizeof(char)));
            CHECK(hipMalloc(&h_sequences[i+1], sy_len * sizeof(char)));
            CHECK(hipMemcpy(h_sequences[i], line1, sx_len * sizeof(char), hipMemcpyHostToDevice));
            CHECK(hipMemcpy(h_sequences[i+1], line2, sy_len * sizeof(char), hipMemcpyHostToDevice));
        }
    }
    CHECK(hipMalloc((void **)&d_sequences, num_of_sequences * sizeof(char *)));
    hipMemcpy(d_sequences, h_sequences, sizeof(char*) * num_of_sequences, hipMemcpyHostToDevice);
    CHECK(hipMalloc((void **)&d_sequences_len, num_of_sequences * sizeof(size_t*)));
    hipMemcpy(d_sequences_len, h_sequences_len, sizeof(size_t*) * num_of_sequences, hipMemcpyHostToDevice);
    CHECK(hipMalloc((void **)&d_error, sizeof(int))); //maybe has to be an array
    CHECK(hipMalloc(&d_result, result_len*sizeof(int)));

    
    three_antidiagonals_size = longest_antidiagonal_len * 3 * 3;
    nBytes = three_antidiagonals_size * sizeof(int);
    //la grandezza della shared memory va data per blocco o in generale?
    //todo inserire una parte che controlla quanto spazio abbiamo in memoria e quando si satura chiama il kernel e
    //se rimangono allineamenti lo richiama una seconda volta.
    dim3 block(block_size);
    dim3 grid(num_of_sequences/2);
    int num_of_threads = block.x*block.y*block.z;
    printf("[main] block_size: %d\n", block.x*block.y*block.z);
    printf("[main] grid_size: %d\n", grid.x*grid.y*grid.z);
    iStart = seconds();
    alignGPU<<<grid, block, nBytes+ num_of_threads*sizeof(int)>>>(d_sequences, d_sequences_len, num_of_sequences, d_result, d_error, three_antidiagonals_size, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP);
    CHECK(hipDeviceSynchronize());
    CHECK(hipGetLastError());
    CHECK(hipMemcpy(h_result, d_result, result_len * sizeof(int), hipMemcpyDeviceToHost));
    CHECK(hipMemcpy(h_error, d_error, sizeof(int), hipMemcpyDeviceToHost));

    iElaps = seconds() - iStart;
    printf("elapsed %f\n", iElaps);
    
    
    output_file = fopen(output_file_path, "a");
    if (output_file == NULL) {
        perror("Error opening file");
        return 1;
    }

    //print the results
    for(int i = 0; i < result_len; i++) {
        fprintf(output_file, "Score: %d\n", h_result[i]);
    }
    // Close the file
    if (fclose(output_file) != 0) {
        perror("Error closing file");
        return 1;
    }


    //free the global memory
    for(int i = 0; i < num_of_sequences; i++){
        hipFree(h_sequences[i]); //host not device memory
        //TODO: check error
    }
    hipFree(d_sequences);
    hipFree(d_result);
    hipFree(d_error);

    //free host memory
    free(h_sequences);
    free(h_result);
    free(h_error);
    fclose(input_file);
    CHECK(hipDeviceReset());
    return 0;
}




/*cambiare in size_t le dimensioni invece che int (code readability)
    fare il controllo se c'è abbastanza memoria per mettere tutto nella GPU. Se no mettere solo una parte.
*/
