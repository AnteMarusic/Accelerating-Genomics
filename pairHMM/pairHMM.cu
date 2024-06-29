#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>
#include <sys/time.h>

#define MAX_READ_LEN 1000
#define OFFSET 33.0

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

//used to find the execution time of the program
inline double seconds(){
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}


// function to get the time of day in seconds
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

void antidiags_print(double *antidiags, int min_len) {
    int k = 1;
    int q = 1;
    printf("\n____________________________\n");
    for(int i = 0; i < (min_len+1)*3*3; i ++) {
        printf("%e,", antidiags[i]);
        if (i + 1 == (min_len+1) * k) {
            if (i + 1 == ((min_len+1) * 3) * q) {
                printf("\n\n");
                q++;
            }
            printf("|");
            k ++;
        }
    }
    printf("\n____________________________\n");
    return;
}


//returns positive values in case of success and the value is in val
//returns negative number in case of failure
__device__ int m_get(double *val, int iy, int ix, double *mat, int read_len, int haplotype_len) {
    int min_len = read_len < haplotype_len ? read_len : haplotype_len;
    int max_antid_dim = min_len + 1;
    int index = 0;
    int cur_antid_num = ix + iy;
    int offset_antid = cur_antid_num % 3;

    //check the input
    if (ix >= haplotype_len+1 || iy >= read_len+1 || ix < 0 || iy < 0) {
        printf("[m_get] out of bounds set iy: %d, ix: %d\n", iy, ix);
        return -1;
    }
    if(cur_antid_num >= read_len+1){
        index = (read_len+1) - 1 - iy + offset_antid*max_antid_dim;
    }
    else {
        index = ix + offset_antid*max_antid_dim;
    }
    *val = mat[index];
    //printf("[m_get] getting val: %e, iy:%d, ix:%d, at index: %d\n", *val, iy, ix, index);
    return 1;
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
__device__ int m_set(double val, int iy, int ix, double *mat, int read_len, int haplotype_len) {
    int min_len = read_len < haplotype_len ? read_len : haplotype_len;
    int max_antid_dim = min_len + 1;
    int index = 0;
    int cur_antid_num = ix + iy;
    int offset_antid = cur_antid_num % 3;

    //check the input
    if (ix >= haplotype_len+1 || iy >= read_len+1 || ix < 0 || iy < 0) {
        printf("[m_set] out of bounds set iy: %d, ix: %d\n", iy, ix);
        return -1;
    }

    if(cur_antid_num >= read_len+1){
        index = (read_len+1) - 1 - iy + offset_antid*max_antid_dim;
    }
    else {
        index = ix + offset_antid*max_antid_dim;
    }

    //  if(cur_antid_num > ny-1){
    //     index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    // }
    // else {
    //     index = ix + offset_antid*antid_dim + offset_mat;
    // }
    //printf("[m_set] setting val: %e, iy:%d, ix:%d, at index: %d\n", val, iy, ix, index);

    mat[index] = val;
    return 1;
}

//input analysis
void partition_read(char *line, char *read, double *Qr, double *Qi, double *Qd, double *Qg, int len) {
    char m[len+1], ins[len+1], d[len+1], g[len+1];
    sscanf(line, "%s %s %s %s %s", read, m, ins, d, g);

    for (int i = 0; i < len; i++) {
        Qr[i] = pow(10.0, -(m[i] - OFFSET) * 0.1);
        Qi[i] = pow(10.0, -(ins[i] - OFFSET) * 0.1);
        Qd[i] = pow(10.0, -(d[i] - OFFSET) * 0.1);
        Qg[i] = pow(10.0, -(g[i] - OFFSET) * 0.1);
    }
}

__device__ double p(char read, char haplotype, double read_quality) {
    return (read == haplotype || read == 'N' || haplotype == 'N') ? 1 - read_quality : read_quality;
}

__device__ double mm(double insert_quality, double delete_quality) {
    return 1 - (insert_quality + delete_quality);
}

void cleanup(double *antidiags, char *read, double *Qr, double *Qi, double *Qd, double *Qg, char **haplotypes, int num_haplotypes) {
    if (antidiags != NULL) {
        free(antidiags);
    }
    if (Qr != NULL) {
        free(Qr);
    }
    if (Qi != NULL) {
        free(Qi);
    }
    if (Qd != NULL) {
        free(Qd);
    }
    if (Qg != NULL) {
        free(Qg);
    }
    if (haplotypes != NULL) {
        for (int i = 0; i < num_haplotypes; i++) {
            if (haplotypes[i] != NULL) free(haplotypes[i]);
        }
        free(haplotypes);
    }
    if (read != NULL) {
        free(read);
    }
}

__host__ void freeCudaAll(char **d_haplotypes, char **d_reads, double **d_Qr, double **d_Qi, double **d_Qd, double **d_Qg, double *d_result, char **h_haplotypes, char **h_reads, double **h_Qr, double **h_Qi, double **h_Qd, double **h_Qg, int num_of_aligmments, int num_haplotypes, int num_read){
    for (int i = 0; i < num_haplotypes; i ++) {
        cudaFree(h_haplotypes[i]);
    }
    cudaFree(d_haplotypes);
    for (int i = 0; i < num_read; i ++) {
        cudaFree(h_reads[i]);
        cudaFree(h_Qr[i]);
        cudaFree(h_Qi[i]);
        cudaFree(h_Qd[i]);
        cudaFree(h_Qg[i]);
    }
    cudaFree(d_reads);
    cudaFree(d_Qr);
    cudaFree(d_Qi);
    cudaFree(d_Qd);
    cudaFree(d_Qg);
    cudaFree(d_result);
    d_haplotypes = NULL;
    d_reads = NULL;
    d_result = NULL;
    d_Qr = NULL;
    d_Qi = NULL;
    d_Qd = NULL;
    d_Qg = NULL;
}

__global__ void printTest() {
    printf("Hello from the GPU!\n");
}

__global__ void printREADS(char **reads, int num_read) {
    for (int i = 0; i < num_read; i ++) {
        printf("%s\n", reads[i]);
    }   
    printf("finished printing\n");
}

__global__ void setRESULT(double *result, int num_of_aligmments) {
    for (int i = 0; i < num_of_aligmments; i ++) {
        result[i] = i;
    }   
    printf("set all values in result\n");
}

__global__ void alignPairs(char **haplotypes, char **reads, int num_of_aligmments, int num_read, int num_haplotypes) {
    int bid = blockIdx.x + blockIdx.y + blockIdx.z; //block index
    int tid = threadIdx.x + threadIdx.y + threadIdx.z; //thread index
    int grid_size = gridDim.x*gridDim.y*gridDim.z;
    int num_of_threads = blockDim.x*blockDim.y*blockDim.z;
    int haplotype_index = -1;
    int read_index = -1;
    if (bid < num_of_aligmments) {
        read_index = bid % num_read;
        haplotype_index = bid/num_read;
        printf("align read[%d]: %s, haplotype[%d]: %s\n", read_index, reads[read_index], haplotype_index, haplotypes[haplotype_index]);
    } 

}
//
__global__ void align1thread(char **haplotypes, char **reads, double *result, int num_of_aligmments, int num_read, int num_haplotypes,
int *haplotypes_len, int *reads_len, double **array_Qr, double **array_Qi, double **array_Qd, double **array_Qg) {
    int bid = blockIdx.x + blockIdx.y + blockIdx.z; //block index
    int tid = threadIdx.x + threadIdx.y + threadIdx.z; //thread index
    int grid_size = gridDim.x*gridDim.y*gridDim.z;
    int num_of_threads = blockDim.x*blockDim.y*blockDim.z;
    if (bid >= num_of_aligmments) {
        return;
    }
    int haplotype_index = bid/num_read;
    int read_index = bid % num_read;
    extern __shared__ double antidiags[];
    int read_len = reads_len[read_index];
    int haplotype_len = haplotypes_len[haplotype_index];
    int antid_num = (read_len+1) + (haplotype_len+1) - 1;
    int min_len = read_len < haplotype_len ? read_len : haplotype_len;
    int max_len = read_len > haplotype_len ? read_len : haplotype_len;
    int yr_phase = min_len; //yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //with each iteration. if the matrix has for example less rows than cols the y phase will last number of rows iterations.
    //red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    int antidiags_size = (min_len+1) * 3 * 3;
    double *M = antidiags;
    double *X = antidiags + (min_len+1)*3;
    double *Y = antidiags + 2*(min_len+1)*3;
    int o_phase = (max_len+1) - (min_len+1);//orange phase is the one where the cur_antid_dim is the max possible
    int i = 0, j = 0, k = 0;
    int cur_antid = 0;
    int gen_i = 0;
    int gen_j = 0;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double value = DBL_MAX/16 / (double)haplotype_len;
    int antid_dim = 1;
    double likelihood = 0;
    char *R = reads[read_index];
    char *H = haplotypes[haplotype_index];
    double *Qr = array_Qr[read_index];
    double *Qi = array_Qi[read_index];
    double *Qd = array_Qd[read_index];
    double *Qg = array_Qg[read_index];



    //printf("align read[%d]: %s, haplotype[%d]: %s\n", read_index, reads[read_index], haplotype_index, haplotypes[haplotype_index]);
    // int antid_num = nx + ny - 1;
    //     int yr_phase = nx - 1;//yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //     //with each iteration. red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    //     int o_phase = ny - nx + 1;//orange phase is the one where the cur_antid_dim is the max possible

        // Declare 2D arrays to store the matrices

    for (cur_antid = 0; cur_antid < antid_num; cur_antid++) {
        // printf("cur_antid: %d, antid_dim: %d\n",cur_antid, antid_dim);
        for (k = 0; k < antid_dim; k++) {
            // printf("i:%d, j:%d\n",i,j);
            if (i == 0) {
                m_set(0, i, j, M, read_len, haplotype_len);
                m_set(0, i, j, X, read_len, haplotype_len);
                m_set(value, i, j, Y, read_len, haplotype_len);
            }
            else if (j == 0) {
                m_set(0, i, j, M, read_len, haplotype_len);
                m_set(0, i, j, X, read_len, haplotype_len);
                m_set(0, i, j, Y, read_len, haplotype_len);
            }
            else {
                //non so se va bene dargli l'indirizzo di temp che è una variabile locale
                m_get(&temp1, i-1, j-1, M, read_len, haplotype_len);
                m_get(&temp2, i-1, j-1, X, read_len, haplotype_len);
                m_get(&temp3, i-1, j-1, Y, read_len, haplotype_len);
                m_set(p(R[i-1], H[j-1], Qr[i-1]) * (mm(Qi[i-1], Qd[i-1]) * temp1 + (1 - Qg[i-1]) * (temp2 + temp3)), i, j, M, read_len, haplotype_len);

                m_get(&temp1, i-1, j, M, read_len, haplotype_len);
                m_get(&temp2, i-1, j, X, read_len, haplotype_len);
                m_set(temp1 * Qi[i-1] + temp2 * Qg[i-1], i, j, X, read_len, haplotype_len);

                m_get(&temp1, i, j-1, M, read_len, haplotype_len);
                m_get(&temp2, i, j-1, Y, read_len, haplotype_len);
                m_set(temp1 * Qd[i-1] + temp2 * Qg[i-1], i, j, Y, read_len, haplotype_len);
            }
            // antidiags_print(M, min_len);

            //compute likelihood
            if (i == read_len){
                // printf("LIKELYHOOD i:%d, j:%d\n",i,j);
                m_get(&temp1, i, j, M, read_len, haplotype_len);
                m_get(&temp2, i, j, X, read_len, haplotype_len);
                // printf("temp1: %e, temp2: %e\n",temp1, temp2);
                likelihood += temp1 + temp2;
            }
            if (i > 0) { //this should always be true
                i--;
            }
            if (j < haplotype_len) { //+2?
                j++;
            }
        }

        
        if(cur_antid < yr_phase) { //at index yr_phase starts the o_phase
            // printf("[main] growing phase\n");
            antid_dim ++;
        }
        else if(cur_antid >= yr_phase + o_phase) {
            // printf("[main] decreasing phase\n");
            antid_dim --;
        }
        else {
            // printf("[main] const phase\n");
        }
        if (gen_i == read_len) {
            gen_j++;
        }
        if (gen_i < read_len) {
            gen_i++;
        }
        i = gen_i;
        j = gen_j;
    }
    result[bid] = log10(likelihood) - log10(DBL_MAX/16);
}

int main(int argc, const char *argv[]) {
    //variables:
    FILE *file_stream_h;
    FILE *file_stream_r;
    FILE *file_stream_out;

    // set up device
    int dev = 1;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("[main] Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    // Check if the correct number of arguments is provided
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file_r> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Open the files with the provided names
    file_stream_r = fopen(argv[1], "r");
    if (file_stream_r == NULL) {
        perror("Error opening input file_r");
        return EXIT_FAILURE;
    }

    file_stream_h = fopen(argv[1], "r");
    if (file_stream_h == NULL) {
        perror("Error opening input file_h");
        fclose(file_stream_r);
        return EXIT_FAILURE;
    }

    file_stream_out = fopen(argv[2], "w");
    if (file_stream_out == NULL) {
        perror("Error opening output file");
        fclose(file_stream_r);
        fclose(file_stream_h);
        return EXIT_FAILURE;
    }

    if (file_stream_r == NULL || file_stream_h == NULL || file_stream_out == NULL) {
        fprintf(stderr, "Error opening files.\n");
        return EXIT_FAILURE;
    }

    int num_read = 0;
    int num_haplotypes = 0;
    int len_read = 0;
    int len_haplotype = 0;
    int i, j;
    int iteration = 1;

    char line[MAX_READ_LEN*5+1];

    double *lh = (double*)malloc(sizeof(double));
    double *Qr = NULL;
    double *Qi = NULL;
    double *Qd = NULL;
    double *Qg = NULL;
    char *read = NULL;

    int min_len = 0;
    int antidiags_size = 0;
    int nBytes = 0;
    double iStart, iElaps;
    int num_of_aligmments = 0;
    int three_antidiagonals_size = 0;

    char **h_haplotypes = NULL;
    char **h_reads = NULL;
    double **h_Qr = NULL;
    double **h_Qi = NULL;
    double **h_Qd = NULL;
    double **h_Qg = NULL;
    double *h_result = NULL;
    int *h_haplotypes_len = NULL;
    int *h_reads_len = NULL;
    char **d_haplotypes = NULL;
    char **d_reads = NULL;
    double **d_Qr = NULL;
    double **d_Qi = NULL;
    double **d_Qd = NULL;
    double **d_Qg = NULL;
    double *d_result = NULL;
    int *d_haplotypes_len = NULL;
    int *d_reads_len = NULL;
    int longest_antidiagonal_len = 0;
 
    //while file contains line
    while (1) {
        printf("#batch: %d\n", iteration);

        //set num_read and num_haplotypes
        if (fgets(line, sizeof(line), file_stream_r) == NULL) {
            break;
        }
        sscanf(line, "%d %d", &num_read, &num_haplotypes);
        fgets(line, sizeof(line), file_stream_h); // Skip the line with the haplotype lengths

        //Allocate memory for haplotypes array
        h_haplotypes = (char**)malloc(num_haplotypes * sizeof(char*));
        if (h_haplotypes == NULL) {
            fprintf(stderr, "Memory allocation failed for haplotypes array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_reads = (char**)malloc(num_read * sizeof(char*));
        if (h_reads == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        //allocate the memory for the sequences, the result
        num_of_aligmments = num_read*num_haplotypes; //number of results (since we have two sequences for each alignment is the half of the sequences)
        h_result = (double *)malloc(num_of_aligmments*sizeof(double));
        if (h_result == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_Qr = (double**)malloc(num_read * sizeof(double*));
        if (h_Qr == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_Qi = (double**)malloc(num_read * sizeof(double*));
        if (h_Qi == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_Qd = (double**)malloc(num_read * sizeof(double*));
        if (h_Qd == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_Qg = (double**)malloc(num_read * sizeof(double*));
        if (h_Qg == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        int *h_haplotypes_len = NULL;
        h_haplotypes_len = (int *)malloc(num_haplotypes * sizeof(int*));
        if (h_haplotypes_len == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }
        h_reads_len = (int *)malloc(num_read * sizeof(int*));     
        if (h_reads_len == NULL) {
            fprintf(stderr, "Memory allocation failed for sequences array\n");
            // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
            return EXIT_FAILURE;
        }

        //put the stream to the first haplotype
        for (i = 0; i < num_read; i++) {
            fgets(line, sizeof(line), file_stream_h);
        }
        // //get all the haplotypes
        // for (i = 0; i < num_haplotypes; i++) {
        //     if (fgets(line, sizeof(line), file_stream_h) == NULL) {
        //         fprintf(stderr, "Error reading haplotypes.\n");
        //         cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
        //         return EXIT_FAILURE;
        //     }
        //     line[strcspn(line, "\n")] = '\0'; // Remove the newline character
        //     haplotypes[i] = (char*)malloc((sizeof(line) + 1) * sizeof(char));
        //     if (haplotypes[i] == NULL) {
        //         fprintf(stderr, "Error allocating memory for haplotype %d.\n", i);
        //         cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
        //         return EXIT_FAILURE;
        //     }
        //     strncpy(haplotypes[i], line, sizeof(line));
        // }

        //get all the haplotypes
        for (i = 0; i < num_haplotypes; i++) {
            if (fgets(line, sizeof(line), file_stream_h) == NULL) {
                fprintf(stderr, "Error reading haplotypes.\n");
                // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            CHECK(cudaMalloc(&h_haplotypes[i], (strlen(line) + 1) * sizeof(char)));
            CHECK(cudaMemcpy(h_haplotypes[i], line, (strlen(line) + 1) * sizeof(char), cudaMemcpyHostToDevice));
            h_haplotypes_len[i] = strlen(line);
            longest_antidiagonal_len = longest_antidiagonal_len < strlen(line) ? strlen(line) : longest_antidiagonal_len;
        }


        //read all the reads
        for (i = 0; i < num_read; i++) {
            if (fgets(line, sizeof(line), file_stream_r) == NULL) {
                fprintf(stderr, "Error reading reads.\n");
                // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            len_read = (strlen(line) - 4) / 5;
            longest_antidiagonal_len = longest_antidiagonal_len < len_read ? len_read : longest_antidiagonal_len;
            h_reads_len[i] = len_read;

            //allocation vectors
            Qr = (double *)malloc(len_read * sizeof(double));
            Qi = (double *)malloc(len_read * sizeof(double));
            Qd = (double *)malloc(len_read * sizeof(double));
            Qg = (double *)malloc(len_read * sizeof(double));
            read = (char *)malloc((len_read + 1) * sizeof(char));
            if (Qr == NULL || Qi == NULL || Qd == NULL || Qg == NULL) {
                fprintf(stderr, "Error allocating memory for quality scores.\n");
                // cleanup(antidiags, read, Qr, Qi, Qd, Qg, h_haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }

            //set values to quality arrays
            partition_read(line, read, Qr, Qi, Qd, Qg, len_read);
            //allocation matrices
            // M = (double **)malloc((len_read + 1) * sizeof(double *));
            // X = (double **)malloc((len_read + 1) * sizeof(double *));
            // Y = (double **)malloc((len_read + 1) * sizeof(double *));

            CHECK(cudaMalloc(&h_reads[i], len_read * sizeof(char)));
            CHECK(cudaMemcpy(h_reads[i], read, len_read * sizeof(char), cudaMemcpyHostToDevice));
            CHECK(cudaMalloc(&h_Qr[i], len_read * sizeof(double)));
            CHECK(cudaMemcpy(h_Qr[i], Qr, len_read * sizeof(double), cudaMemcpyHostToDevice));
            CHECK(cudaMalloc(&h_Qi[i], len_read * sizeof(double)));
            CHECK(cudaMemcpy(h_Qi[i], Qi, len_read * sizeof(double), cudaMemcpyHostToDevice));
            CHECK(cudaMalloc(&h_Qd[i], len_read * sizeof(double)));
            CHECK(cudaMemcpy(h_Qd[i], Qd, len_read * sizeof(double), cudaMemcpyHostToDevice));
            CHECK(cudaMalloc(&h_Qg[i], len_read * sizeof(double)));
            CHECK(cudaMemcpy(h_Qg[i], Qg, len_read * sizeof(double), cudaMemcpyHostToDevice));

            free(Qr);
            free(Qi);
            free(Qd);
            free(Qg);
            free(read);
        }

        //allocate the arrays in the device
        CHECK(cudaMalloc((void **)&d_haplotypes, num_read * sizeof(char *)));
        cudaMemcpy(d_haplotypes, h_haplotypes, sizeof(char*) * num_haplotypes, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_reads, num_read * sizeof(char *)));
        cudaMemcpy(d_reads, h_reads, sizeof(char*) * num_read, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_Qr, num_read * sizeof(double *)));
        cudaMemcpy(d_Qr, h_Qr, sizeof(double *) * num_read, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_Qi, num_read * sizeof(double *)));
        cudaMemcpy(d_Qi, h_Qi, sizeof(double *) * num_read, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_Qd, num_read * sizeof(double *)));
        cudaMemcpy(d_Qd, h_Qd, sizeof(double *) * num_read, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_Qg, num_read * sizeof(double *)));
        cudaMemcpy(d_Qg, h_Qg, sizeof(double *) * num_read, cudaMemcpyHostToDevice);
        //allocate the result
        CHECK(cudaMalloc((void **)&d_result, num_of_aligmments * sizeof(double)));
        CHECK(cudaMalloc((void **)&d_haplotypes_len, num_haplotypes * sizeof(int)));
        cudaMemcpy(d_haplotypes_len, h_haplotypes_len, sizeof(int) * num_haplotypes, cudaMemcpyHostToDevice);
        CHECK(cudaMalloc((void **)&d_reads_len, num_read * sizeof(int)));
        cudaMemcpy(d_reads_len, h_reads_len, sizeof(int) * num_read, cudaMemcpyHostToDevice);


        

        iStart = seconds();
        //la grandezza della shared memory va data per blocco o in generale?
        //todo inserire una parte che controlla quanto spazio abbiamo in memoria e quando si satura chiama il kernel e
        //se rimangono allineamenti lo richiama una seconda volta.
        // dim3 block(32);
        // dim3 grid(num_of_aligmments);
        // int num_of_threads = block.x*block.y*block.z;
        // printf("[main] block_size: %d\n", block.x*block.y*block.z);
        // printf("[main] grid_size: %d\n", grid.x*grid.y*grid.z);
        dim3 grid(num_of_aligmments);
        three_antidiagonals_size = longest_antidiagonal_len * 3 * 3;
        nBytes = three_antidiagonals_size * sizeof(double);
        iStart = seconds();

        align1thread<<<grid,1, nBytes>>>(d_haplotypes, d_reads, d_result, num_of_aligmments, num_read, num_haplotypes, d_haplotypes_len, d_reads_len,
            d_Qr, d_Qi, d_Qd, d_Qg);
        CHECK(cudaDeviceSynchronize());
        CHECK(cudaGetLastError());
        CHECK(cudaMemcpy(h_result, d_result, num_of_aligmments * sizeof(double), cudaMemcpyDeviceToHost));
        for (i = 0; i < num_of_aligmments; i ++) {
            //printf("%e\n", h_result[i]);
            fprintf(file_stream_out, "%f\n", h_result[i]);
        }
        printf("finished printing\n");
        //freeeeee
        freeCudaAll(d_haplotypes, d_reads, d_Qr, d_Qi, d_Qd, d_Qg, d_result, h_haplotypes, h_reads, h_Qr, h_Qi, h_Qd, h_Qg, num_of_aligmments, num_haplotypes, num_read);
        //put the stream pointing to the reads at the beginning of the next iteration
        for (i = 0; i < num_haplotypes; i++) {
            if (fgets(line, sizeof(line), file_stream_r) == NULL) {
                fprintf(stderr, "Error repositioning read file stream.\n");
                return EXIT_FAILURE;
            }
        }

        iteration++;
    }

    //close files
    fclose(file_stream_h);
    fclose(file_stream_r);
    fclose(file_stream_out);
    CHECK(cudaDeviceReset());

    return EXIT_SUCCESS;
}
