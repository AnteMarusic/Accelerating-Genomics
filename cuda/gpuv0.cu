#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <string.h> // Include the string.h header for strcspn function
#include <limits.h>

#define TRUE 1
#define FALSE 0

//#undef DEBUG
#define DEBUG

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
#define MAX_SIZE 10000
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
    printf("\n____________________________\n");
    for(int i = 0; i < len; i ++) {
        a[i] == NEGATIVE_INFINITY ? printf("-∞,") : printf("%d,", a[i]);
        if (i + 1 == (nx) * k) {
            if (i + 1 == (nx * 3) * q) {
                printf("\n");
                q++;
            }
            printf("|");
            k ++;
        }
    }
    printf("\n____________________________\n");
    return;
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
    printf("[m_get thread %d] getting val: %d, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", threadIdx.x, *val, iy, ix, mat, index);
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

    printf("[m_set thread %d] setting val: %d, iy:%d, ix:%d, at index: %d on matrix %d\n", threadIdx.x, *val, iy, ix, index, mat);
    antidiags[index] = *val;
    return 1;
}


__global__ void alignGPU(char *d_sx, char *d_sy, int *d_result, int *d_error, int antidiags_size, int sx_len, int sy_len, int nx, int ny, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
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
                printf("out of bounds iy = %d, ix = %d\n", iy, ix);
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
        if(threadIdx.x == 0) {
            printf("cur_antd_num: %d\n", cur_antid_num);
            a_print(antidiags, antidiags_size, nx);
        }
    }

    //COMPUTE MAX
    for (int s = blockDim.x * blockDim.y / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s && threadIdx.x + s < blockDim.x * blockDim.y) {
        shared_max[threadIdx.x] = fmaxf(shared_max[threadIdx.x], shared_max[threadIdx.x + s]);
    }
    __syncthreads();
    }
    // The first thread in the block writes the result to global memory
    if (threadIdx.x == 0) {
        atomicMax(d_result, shared_max[0]);
    }
    //if(threadIdx.x == 0) printf("max: %d\n", max);
    //*d_result = max;
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

    if (argc != 3) {
        printf("Usage: %s <string1> <string2>\n", argv[0]);
        return 1;
    }
    char *h_sx = argv[1];
    char *h_sy = argv[2];
    int sx_len = strlen(h_sx);
    int sy_len = strlen(h_sy);

    if (sx_len > sy_len) {
        char *tmp_sx = (char *)malloc(sx_len+1*sizeof(int));
        char *tmp_sy = (char *)malloc(sy_len+1*sizeof(int));
        strncpy(tmp_sx, h_sx, sx_len);
        strncpy(tmp_sy, h_sy, sy_len);
        h_sy = (char *)malloc(sx_len+1*sizeof(int));
        h_sx = (char *)malloc(sy_len+1*sizeof(int));
        strncpy(h_sx, tmp_sy, sy_len);
        strncpy(h_sy, tmp_sx, sx_len);
        int tmp = sx_len;
        sx_len = sy_len;
        sy_len = tmp;
        free(tmp_sx);
        free(tmp_sy);
    }
    
    printf("[main] sx: %s, xy: %s, of size %d, %d respectively\n", h_sx, h_sy, sx_len, sy_len);
    int nx = sx_len+1;
    int ny = sy_len+1;
    printf("[main] nx: %d, ny: %d\n", nx, ny);
    reverse_array(h_sy, sy_len);
    printf("[main] sx: %s, xy: %s, of size %d, %d respectively\n", h_sx, h_sy, sx_len, sy_len);

    double iStart;
    double iElaps;
    char *d_sx;
    char *d_sy;
    int *d_result;
    int *d_error;
    CHECK(cudaMalloc((void **)&d_sx, sx_len));
    CHECK(cudaMalloc((void **)&d_sy, sy_len));
    CHECK(cudaMalloc((void **)&d_result, sizeof(int)));
    CHECK(cudaMalloc((void **)&d_error, sizeof(int)));
    CHECK(cudaMemcpy(d_sx, h_sx, sx_len, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_sy, h_sy, sy_len, cudaMemcpyHostToDevice));
    int antidiags_size = nx * 3 * 3;
    int nBytes = antidiags_size * sizeof(int);

    dim3 block(32);
    dim3 grid(1);
    //printf("[main] GPU invocation\n");
    iStart = seconds();
    int num_of_threads = block.x*block.y*block.z;
    alignGPU <<<grid, block, nBytes+ num_of_threads*sizeof(int)>>>(d_sx, d_sy, d_result, d_error, antidiags_size, sx_len, sy_len, nx, ny, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP); 
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    printf("[main] alignGPU <<<(%d,%d), (%d,%d), %d>>> elapsed %f sec\n", grid.x, grid.y, block.x, block.y, nBytes, iElaps);
    printf("elapsed %f\n", iElaps);
    // check kernel error
    CHECK(cudaGetLastError());

    //Copy data back from GPU memory to CPU memory.
    int *h_result = (int *)malloc(sizeof(int));
    CHECK(cudaMemcpy(h_result, d_result, sizeof(int), cudaMemcpyDeviceToHost));
    printf("Score: %d\n", *h_result);
    int *h_error = (int *)malloc(sizeof(int));
    CHECK(cudaMemcpy(h_error, d_error, sizeof(int), cudaMemcpyDeviceToHost));
    //printf("[main] Error: %d\n", *h_error);

    // free device global memory
    CHECK(cudaFree(d_result));
    CHECK(cudaFree(d_sy));
    CHECK(cudaFree(d_sx));

    // reset device
    CHECK(cudaDeviceReset());
    return 0;
}
