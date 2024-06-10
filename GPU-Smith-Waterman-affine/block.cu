#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <string.h> // Include the string.h header for strcspn function
#include <limits.h>

#undef DEBUG
//#define DEBUG

#ifdef DEBUG
#define PRINT printf
#else
#define PRINT // macros
#endif

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

inline double seconds()
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#define NEGATIVE_INFINITY INT_MIN  // Minimum value for an int
#define MAX_LINE_LENGTH 10000
#define SCORE_MATCH 1
#define SCORE_MISMATCH (-1)
#define SCORE_OPEN_GAP (-3)
#define SCORE_EXTEND_GAP (-1)
#define P 1
#define Q 2
#define D 3


/*
Example of how the strings are in the matrixes
s1 = CG
s2 = CCGA
nx is the size s2 + 1
    _CCGA nx
_   
C    0000
G    0000
ny
*/

__device__ void a_print(int *a, int len, int nx) {
    int k = 1;
    int q = 1;
    PRINT("\n____________________________\n");
    for(int i = 0; i < len; i ++) {
        a[i] == NEGATIVE_INFINITY ? printf("-∞,") : printf("%d,", a[i]);
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

__device__ void array_print(int *a, int len) {
    PRINT("\n____________________________\n");
    for(int i = 0; i<len; i++) {
        PRINT("%d ", a[i]);
    }
    PRINT("\n____________________________\n");
}

__device__ size_t device_strlen(const char* str) {
    size_t len = 0;
    while (str[len] != '\0') {
        len++;
    }
    return len;
}


__device__ int max_of_two_integers(int a, int b) {
    return (a > b) ? a : b;
}

__device__ int sum_with_infinity (int a, int b) {
    return (a == NEGATIVE_INFINITY || b == NEGATIVE_INFINITY) ? NEGATIVE_INFINITY : a + b;
}

__device__ int max_of_four_integers (int a, int b, int c, int d) {
    return max_of_two_integers (max_of_two_integers (a, b), max_of_two_integers(c, d));
}

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
__device__ int m_get(int *val, int iy, int ix, int mat, int *antidiags, int nx, int ny) {
    int antid_dim = nx;
    int offset_antid = 0;
    int offset_mat = 0;
    int index = 0;
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        printf("[m_get thread %d] out of bounds set iy: %d, ix: %d\n",threadIdx.x, iy, ix);
        return -1;
    }
    /*
    if (ix + iy != cur_antid_num) {
        printf("[m_get] iy: %d, ix: %d not in mem\n", iy, ix);
        return -1;
    }
    */
    
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
            printf("[m_get thread %d] mat doesn't exist\n", threadIdx.x);
    }
    int cur_antid_num = ix + iy;
    offset_antid = cur_antid_num % 3;//three is the num of antid saved at the same time
    
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }
    *val = antidiags[index];
    PRINT("[m_get thread %d] getting val: %d, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", threadIdx.x, *val, iy, ix, mat, index);
    return 1;
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
__device__ int m_set(int *val, int iy, int ix, int mat, int *antidiags, int nx, int ny) {
    int antid_dim = nx;
    int offset_antid = 0;
    int offset_mat = 0;
    int index = 0;

    //check the input
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        printf("[m_set thread %d] out of bounds set iy: %d, ix: %d\n", threadIdx.x, iy, ix);
        return -1;
    }
    /*if (ix + iy != cur_antid_num) {
        printf("[m_set] iy:%d, ix:%d not in mem\n", iy, ix);
        return -1;
    }*/
    
    //compute the offset_mat
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
            printf("[m_set thread %d] mat doesn't exist\n", threadIdx.x);
    }
    int cur_antid_num = ix + iy;
    //three is the num of antid saved at the same time
    offset_antid = cur_antid_num % 3;
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }

    PRINT("[m_set thread %d] setting val: %d, iy:%d, ix:%d, at index: %d on matrix %d\n", threadIdx.x, *val, iy, ix, index, mat);
    antidiags[index] = *val;
    return 1;
}


__global__ void alignGPU(char **sequences, int sequences_num, int *d_result, int *d_error, int antidiags_size, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
    int bid = blockIdx.x + blockIdx.y + blockIdx.z;
    int tid = threadIdx.x + threadIdx.y + threadIdx.z;
    int grid_size = gridDim.x+gridDim.y+gridDim.z;
    if(bid == 0 && tid == 0) {
        PRINT("bid: %d, tid: %d\n", bid, tid);
        array_print(d_result, grid_size);
    }
    //safety check to avoid accessing unwanted memory areas if the gridsize has not been correctly set.
    if (bid >= sequences_num) {*d_error = -1; return;} 
    char *d_sx = sequences[bid*2];
    char *d_sy = sequences[bid*2+1];
    if(tid==0){
        PRINT("[bid: %d, tid: %d] sequences[%d] = %s\nsequences[%d] = %s\n", bid, tid, bid*2, sequences[bid*2], bid*2+1, sequences[bid*2+1]);
    }
    int sx_len = device_strlen(sequences[bid*2]);
    int sy_len = device_strlen(sequences[bid*2+1]);
    int nx = sx_len+1;
    int ny = sy_len+1;
    extern __shared__ int s[];
    int *antidiags = s;
    int *shared_max = &s[antidiags_size];
    //__shared__ int max; //check if I can do this. How to collect the max value between all threads? I can do it at the end I think
    int num_of_threads = blockDim.x*blockDim.y*blockDim.z;
    //printf("num_of_threads: %d\n", num_of_threads);
    int antid_num = nx + ny - 1;
    int yr_phase = nx - 1;//yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //with each iteration. red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    int o_phase = ny - nx + 1;//orange phase is the one where the cur_antid_dim is the max possible
    int err;
    int val = 1;
    int minus_infty = NEGATIVE_INFINITY;
    int zero=0;
    int iy = 0, ix = 0;
    int cur_antid_dim;
    int gen_iy = 0;
    int gen_ix = 0;
    int temp1,temp2,temp3;
    //if (threadIdx.x == 0) max = 0;
    shared_max[threadIdx.x] = 0;
    for(int cur_antid_num = 0; cur_antid_num < antid_num; cur_antid_num  ++) {
        //compute dimension of the current antidiagonal
        if(cur_antid_num <= yr_phase) { //at index yr_phase starts the o_phase
            //printf("[kernel %d] growing phase\n", threadIdx.x);
            cur_antid_dim = cur_antid_num + 1;
        }
        else if(cur_antid_num >= yr_phase + o_phase) {
            //printf("[kernel %d] decreasing phase\n", threadIdx.x);
            cur_antid_dim = cur_antid_dim - 1;
        }
        else {
            //printf("[kernel] const phase\n");
        }
        //printf("[kernel %d] cur_antid_num: %d, cur_antid_dim: %d\n", threadIdx.x, cur_antid_num, cur_antid_dim);
        //compute the start of ix for the current antidiagonal
        if(cur_antid_num >= yr_phase + o_phase) {
            gen_ix ++;
        }
        iy = gen_iy - threadIdx.x;
        ix = gen_ix + threadIdx.x;
        while(true) {
            /*
            if (ix + gen_ix >= cur_antid_dim || iy - gen_iy < 0) {
                printf("out of bounds iy = %d, ix = %d\n", iy, ix);
                break;
            }
            */
            if (ix < 0 || iy < 0 || ix >= nx || iy >= ny) {
                PRINT("out of bounds iy = %d, ix = %d\n", iy, ix);
                break;
            }
            
            //set the first row
            if (iy == 0) {
                err = m_set(&minus_infty, iy, ix, P, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, Q, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, D, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
            }
            //set the first col
            else if(ix==0){
                err = m_set(&zero, iy, ix, P, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&minus_infty, iy, ix, Q, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_set(&zero, iy, ix, D, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
            }
            //dynamic programming
            else {
                err = m_get(&temp1, iy-1, ix, D, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy-1, ix, P, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                int b = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&b, iy, ix, P, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}

                err = m_get(&temp1, iy, ix-1, D, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy, ix-1, Q, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                int a = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&a, iy, ix, Q, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}

                err = m_get(&temp1, iy, ix, P, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp2, iy, ix, Q, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                err = m_get(&temp3, iy-1, ix-1, D, antidiags, nx, ny);
                if(err<1) {*d_error = -1; return;}
                //temp3 = d_sy[iy - 1] == d_sx[ix - 1] ? temp3 + SCORE_MATCH : temp3 + SCORE_MISMATCH;
                //memorise the string in opposite direction to increase locality
                temp3 = d_sy[sy_len - iy] == d_sx[ix - 1] ? temp3 + SCORE_MATCH : temp3 + SCORE_MISMATCH;
                val = max_of_four_integers(temp1, temp2, temp3, 0);
                err = m_set(&val, iy, ix, D, antidiags, nx, ny);
                shared_max[threadIdx.x] = val > shared_max[threadIdx.x] ? val : shared_max[threadIdx.x];
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
        /*
        if(threadIdx.x == 0) {
            printf("cur_antd_num: %d\n", cur_antid_num);
            a_print(antidiags, antidiags_size, nx);
        }
        */
        
    }

    //COMPUTE MAX
    for (int s = blockDim.x * blockDim.y / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s && threadIdx.x + s < blockDim.x * blockDim.y) {
        shared_max[threadIdx.x] = fmaxf(shared_max[threadIdx.x], shared_max[threadIdx.x + s]);
    }
    __syncthreads();
    }
    // The first thread in the block writes the result to global memory
    if (tid == 0) {
        atomicMax(&d_result[bid], shared_max[0]);
    }
    //if(threadIdx.x == 0) printf("max: %d\n", max);
    //*d_result = max;
}

__global__ void test(char **sequences, int num_of_sequences, int *d_result) {
    int bid = blockIdx.x + blockIdx.y + blockIdx.z;
    int tid = threadIdx.x + threadIdx.y + threadIdx.z;
    if(bid<num_of_sequences/2 && tid==0) {
        printf("block: %d, thread: %d\n[device] sequences[%d]: %s\n [device] sequences[%d]: %s\n\n", bid, tid, bid*2, sequences[bid*2], bid*2+1, sequences[bid*2+1]);
        d_result[bid] = bid;
    }
}


//todo maybe are necessary controls on the length of the input
//todo what happens if _ are present in the input string?
//int main()
int main(int argc, char *argv[]) {
     // set up device
    int dev = 1;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("[main] Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <file_path>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line1[MAX_LINE_LENGTH];
    char line2[MAX_LINE_LENGTH];
    int max_antidiag = 0;//len of the longest antidiagonal

    if (fgets(line1, sizeof(line1), file) == NULL) {
        printf("file is empty");
        return 1; // Exit loop if EOF or error
    }
    int line_num = atoi(line1); //get the number of lines of the file
    printf("line_num: %d\n", line_num);
    int sx_len, sy_len;
    char **d_sequences;
    char **h_sequences = (char**)malloc(line_num * sizeof(char*));
    int *d_result;
    int result_len = line_num/2;
    CHECK(cudaMalloc(&d_result, result_len*sizeof(int)));
    int *h_result = (int *)malloc(result_len*sizeof(int));
    for (int i = 0; i < line_num; i+=2) {
        PRINT("%d\n", i);
        // Read the first line
        if (fgets(line1, sizeof(line1), file) == NULL) {
            break; // Exit loop if EOF or error
        }
        // Read the second line
        if (fgets(line2, sizeof(line2), file) == NULL) {
            // Print the first line only if EOF or error on the second line
            printf("%s", line1);
            break;
        }
        PRINT("[host] line1: %s\n line2: %s\n", line1, line2);
        max_antidiag = strlen(line1) > max_antidiag ? strlen(line1) : max_antidiag;
        max_antidiag = strlen(line2) > max_antidiag ? strlen(line2) : max_antidiag;
        if(strlen(line1) > strlen(line2)){
            sx_len = strlen(line2);
            sy_len = strlen(line1);
            reverse_array(line1, sy_len);
            CHECK(cudaMalloc(&h_sequences[i], sx_len * sizeof(char)));
            CHECK(cudaMalloc(&h_sequences[i+1], sy_len * sizeof(char)));
            CHECK(cudaMemcpy(h_sequences[i], line2, sx_len * sizeof(char), cudaMemcpyHostToDevice));
            CHECK(cudaMemcpy(h_sequences[i+1], line1, sy_len * sizeof(char), cudaMemcpyHostToDevice));
        }
        else{
            sx_len = strlen(line1);
            sy_len = strlen(line2);
            reverse_array(line2, sy_len);
            CHECK(cudaMalloc(&h_sequences[i], sx_len * sizeof(char)));
            CHECK(cudaMalloc(&h_sequences[i+1], sy_len * sizeof(char)));
            CHECK(cudaMemcpy(h_sequences[i], line1, sx_len * sizeof(char), cudaMemcpyHostToDevice));
            CHECK(cudaMemcpy(h_sequences[i+1], line2, sy_len * sizeof(char), cudaMemcpyHostToDevice));
        }
    }
    CHECK(cudaMalloc((void **)&d_sequences, line_num * sizeof(char *)));
    cudaMemcpy(d_sequences, h_sequences, sizeof(char*) * line_num, cudaMemcpyHostToDevice);

    int *d_error;
    CHECK(cudaMalloc((void **)&d_error, sizeof(int)));
    int max_antidiags_size = max_antidiag * 3 * 3;
    int nBytes = max_antidiags_size * sizeof(int);
    dim3 block(32);
    dim3 grid(line_num/2);
    int num_of_threads = block.x*block.y*block.z;
    printf("[main] block_size: %d\n", block.x*block.y*block.z);
    printf("[main] grid_size: %d\n", grid.x*grid.y*grid.z);
    double iStart = seconds();
    //la grandezza della shared memory va data per blocco o in generale?
    //todo inserire una parte che controlla quanto spazio abbiamo in memoria e quando si satura chiama il kernel e
    //se rimangono allineamenti lo richiama una seconda volta.
    alignGPU<<<grid, block, nBytes+ num_of_threads*sizeof(int)>>>(d_sequences, line_num, d_result, d_error, max_antidiags_size, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP);
    CHECK(cudaDeviceSynchronize());
    CHECK(cudaGetLastError());
    CHECK(cudaMemcpy(h_result, d_result, result_len * sizeof(int), cudaMemcpyDeviceToHost));
    int *h_error = (int *)malloc(sizeof(int));
    CHECK(cudaMemcpy(h_error, d_error, sizeof(int), cudaMemcpyDeviceToHost));

    for(int i = 0; i < result_len; i++) {
        printf("Score: %d\n", h_result[i]);
    }
    double iElaps = seconds() - iStart;
    printf("elapsed %f\n", iElaps);
    //write a kernel that simply prints the strings
    //write the strings in an array of strings
    //allocate the string and then a pointer to the string. Is it possible?
    for(int i = 0; i < line_num; i++){
        cudaFree(h_sequences[i]); //host not device memory
        //TODO: check error
    }
    cudaFree(d_sequences);
    free(h_sequences);
    fclose(file);
    CHECK(cudaDeviceReset());
    return 0;
}


/*
run commands
nvcc cur.cu -o cur
./cur input.txt
*/

//THERE IS SOME KIND OF MISTAKE IN THE HANDLING OF THE THREADS.
//TOO MANY THREADS ARE ALLOCATED
