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

void a_print(int *a, int len, int nx) {
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


int max_of_two_integers(int a, int b) {
    return (a > b) ? a : b;
}

int sum_with_infinity (int a, int b) {
    return (a == NEGATIVE_INFINITY || b == NEGATIVE_INFINITY) ? NEGATIVE_INFINITY : a + b;
}

int max_of_four_integers (int a, int b, int c, int d) {
    return max_of_two_integers (max_of_two_integers (a, b), max_of_two_integers(c, d));
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
int m_get(int *val, int iy, int ix, int mat, int *antidiags, int cur_antid_num, int nx, int ny) {
    int antid_dim = nx;
    int offset_antid = 0;
    int offset_mat = 0;
    int index = 0;
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        PRINT("[m_get] out of bounds set iy: %d, ix: %d\n", iy, ix);
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
            printf("[m_get] mat doesn't exist\n");
    }
    cur_antid_num = ix + iy;
    offset_antid = cur_antid_num % 3;//three is the num of antid saved at the same time
    
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }
    *val = antidiags[index];
    printf("[m_get] getting val: %d, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", *val, iy, ix, mat, index);
    return 1;
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
int m_set(int *val, int iy, int ix, int mat, int *antidiags, int cur_antid_num, int nx, int ny) {
    int antid_dim = nx;
    int offset_antid = 0;
    int offset_mat = 0;
    int index = 0;

    //check the input
    if (ix >= nx || iy >= ny|| ix < 0 || iy < 0) {
        printf("[m_set] out of bounds set iy: %d, ix: %d\n", iy, ix);
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
            printf("[m_set] mat doesn't exist\n");
    }
    cur_antid_num = ix + iy;
    //three is the num of antid saved at the same time
    offset_antid = cur_antid_num % 3;
    if(cur_antid_num > ny-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }

    printf("[m_set] setting val: %d, iy:%d, ix:%d, at index: %d on matrix %d\n", *val, iy, ix, index, mat);
    antidiags[index] = *val;
    return 1;
}

//todo maybe are necessary controls on the length of the input
//todo what happens if _ are present in the input string?

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <string1> <string2>\n", argv[0]);
        return 1;
    }
    char *sx = argv[1];
    char *sy = argv[2];
    int sx_len = strlen(sx);
    int sy_len = strlen(sy);

    printf("sx: %s, sy: %s\n", sx, sy);

    if (sx_len > sy_len) {
        char *tmp_sx = (char *)malloc(sx_len+1*sizeof(int));
        char *tmp_sy = (char *)malloc(sy_len+1*sizeof(int));
        strncpy(tmp_sx, sx, sx_len);
        strncpy(tmp_sy, sy, sy_len);
        sy = (char *)malloc(sx_len+1*sizeof(int));
        sx = (char *)malloc(sy_len+1*sizeof(int));
        strncpy(sx, tmp_sy, sy_len);
        strncpy(sy, tmp_sx, sx_len);
        int tmp = sx_len;
        sx_len = sy_len;
        sy_len = tmp;
        free(tmp_sx);
        free(tmp_sy);
    }

    printf("sx: %s, sy: %s\n", sx, sy);

    
    //printf("[main] sx: %s, xy: %s, of size %d, %d respectively\n", sx, sy, sx_len, sy_len);
    int nx = sx_len+1;
    int ny = sy_len+1;
    //printf("[main] nx: %d, ny: %d\n", nx, ny);
    //allocate an array to fit three antidiag of max size min_len for the three matrixes
    int antidiags_size = nx * 3 * 3;
    int nBytes = antidiags_size * sizeof(int);
    int *antidiags = (int *) malloc (nBytes);

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
    double iStart;
    double iElaps;
   
    int gen_iy = 0;
    int gen_ix = 0;
    int temp1,temp2,temp3,max=0;
    iStart = seconds();
    for(int cur_antid_num = 0; cur_antid_num  < antid_num; cur_antid_num  ++) {
        if(cur_antid_num  <= yr_phase) { //at index yr_phase starts the o_phase
            //printf("[main] growing phase\n");
            cur_antid_dim = cur_antid_num + 1;
        }
        else if(cur_antid_num >= yr_phase + o_phase) {
            //printf("[main] decreasing phase\n");
            cur_antid_dim = cur_antid_dim - 1;
        }
        else {
            //printf("[main] const phase\n");
        }
        //printf("[main] cur_antid_num: %d, cur_antid_dim: %d\n", cur_antid_num, cur_antid_dim);
        if(cur_antid_num >= yr_phase + o_phase) {
            gen_ix ++;
        }
        iy = gen_iy;
        ix = gen_ix;
        for (int j = 0; j < cur_antid_dim; j ++) {
            //set the first row
            if (iy == 0) {
                err = m_set(&minus_infty, iy, ix, P, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_set(&zero, iy, ix, Q, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_set(&zero, iy, ix, D, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
            }
            //set the first col
            else if(ix==0){
                err = m_set(&zero, iy, ix, P, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_set(&minus_infty, iy, ix, Q, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_set(&zero, iy, ix, D, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
            }
            //dynamic programming
            else {
                err = m_get(&temp1, iy-1, ix, D, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_get(&temp2, iy-1, ix, P, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                int b = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&b, iy, ix, P, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;

                err = m_get(&temp1, iy, ix-1, D, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_get(&temp2, iy, ix-1, Q, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                int a = max_of_two_integers(sum_with_infinity(temp1, SCORE_OPEN_GAP + SCORE_EXTEND_GAP), sum_with_infinity(temp2, SCORE_EXTEND_GAP));
                err = m_set(&a, iy, ix, Q, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;

                err = m_get(&temp1, iy, ix, P, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_get(&temp2, iy, ix, Q, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                err = m_get(&temp3, iy-1, ix-1, D, antidiags, cur_antid_num, nx, ny);
                if(err<1) return -1;
                printf("sy[iy-1]:%c, sx[ix-1]:%c\n",sy[iy - 1], sx[ix - 1]);
                temp3 = sy[iy - 1] == sx[ix - 1] ? temp3 + SCORE_MATCH : temp3 + SCORE_MISMATCH;
                val = max_of_four_integers(temp1, temp2, temp3, 0);
                err = m_set(&val, iy, ix, D, antidiags, cur_antid_num, nx, ny);
                max = val > max ? val : max;
            }
            //update the matrix indexes
            iy --;
            ix ++;
        }
        //move the starting point in case we are in the growing part of the matrix
        if (gen_iy < ny - 1) {
            gen_iy ++;
        }
        printf("cur_antd_num: %d\n", cur_antid_num);
        a_print(antidiags, antidiags_size, nx);
    }
    iElaps = seconds() - iStart;
    printf("elapsed %f\n", iElaps);
    printf("Score: %d\n", max);
    return 0;
}
