#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

#define MAX_READ_LEN 1000
#define OFFSET 33.0


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
int m_get(double *val, int iy, int ix, double *mat, int read_len, int haplotype_len) {
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
int m_set(double val, int iy, int ix, double *mat, int read_len, int haplotype_len) {
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

double p(char read, char haplotype, double read_quality) {
    return (read == haplotype || read == 'N' || haplotype == 'N') ? 1 - read_quality : read_quality;
}

double mm(double insert_quality, double delete_quality) {
    return 1 - (insert_quality + delete_quality);
}

//PairHMM algo
void pairHMM(double *likelihood, double *M, double *X, double *Y, char *R, char *H, int read_len, int haplotype_len, double *Qr, double *Qi, double *Qd, double *Qg) {
    int antid_num = (read_len+1) + (haplotype_len+1) - 1;
    int min_len = read_len < haplotype_len ? read_len : haplotype_len;
    int max_len = read_len > haplotype_len ? read_len : haplotype_len;
    int yr_phase = min_len; //yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //with each iteration. if the matrix has for example less rows than cols the y phase will last number of rows iterations.
    //red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
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

    // int antid_num = nx + ny - 1;
    //     int yr_phase = nx - 1;//yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //     //with each iteration. red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    //     int o_phase = ny - nx + 1;//orange phase is the one where the cur_antid_dim is the max possible

        // Declare 2D arrays to store the matrices
    double **M_matrix = (double **)malloc((read_len + 1) * sizeof(double *));
    double **X_matrix = (double **)malloc((read_len + 1) * sizeof(double *));
    double **Y_matrix = (double **)malloc((read_len + 1) * sizeof(double *));
    for (int i = 0; i <= read_len; i++) {
        M_matrix[i] = (double *)malloc((haplotype_len + 1) * sizeof(double));
        X_matrix[i] = (double *)malloc((haplotype_len + 1) * sizeof(double));
        Y_matrix[i] = (double *)malloc((haplotype_len + 1) * sizeof(double));
    }

    for (cur_antid = 0; cur_antid < antid_num; cur_antid++) {
        // printf("cur_antid: %d, antid_dim: %d\n",cur_antid, antid_dim);
        for (k = 0; k < antid_dim; k++) {
            // printf("i:%d, j:%d\n",i,j);
            if (i == 0) {
                m_set(0, i, j, M, read_len, haplotype_len);
                m_set(0, i, j, X, read_len, haplotype_len);
                m_set(value, i, j, Y, read_len, haplotype_len);


                // Save values to the 2D arrays
                M_matrix[i][j] = 0;
                X_matrix[i][j] = 0;
                Y_matrix[i][j] = value;
            }
            else if (j == 0) {
                m_set(0, i, j, M, read_len, haplotype_len);
                m_set(0, i, j, X, read_len, haplotype_len);
                m_set(0, i, j, Y, read_len, haplotype_len);
                //int m_set(double *val, int iy, int ix, double *mat, int read_len, int haplotype_len) {

                // Save values to the 2D arrays
                M_matrix[i][j] = 0;
                X_matrix[i][j] = 0;
                Y_matrix[i][j] = 0;
            }
            else {
                //non so se va bene dargli l'indirizzo di temp che è una variabile locale
                m_get(&temp1, i-1, j-1, M, read_len, haplotype_len);
                m_get(&temp2, i-1, j-1, X, read_len, haplotype_len);
                m_get(&temp3, i-1, j-1, Y, read_len, haplotype_len);
                double m_val = p(R[i - 1], H[j - 1], Qr[i - 1]) * (mm(Qi[i - 1], Qd[i - 1]) * temp1 + (1 - Qg[i - 1]) * (temp2 + temp3));
                m_set(p(R[i-1], H[j-1], Qr[i-1]) * (mm(Qi[i-1], Qd[i-1]) * temp1 + (1 - Qg[i-1]) * (temp2 + temp3)), i, j, M, read_len, haplotype_len);

                m_get(&temp1, i-1, j, M, read_len, haplotype_len);
                m_get(&temp2, i-1, j, X, read_len, haplotype_len);
                double x_val = temp1 * Qi[i - 1] + temp2 * Qg[i - 1];
                m_set(temp1 * Qi[i-1] + temp2 * Qg[i-1], i, j, X, read_len, haplotype_len);

                m_get(&temp1, i, j-1, M, read_len, haplotype_len);
                m_get(&temp2, i, j-1, Y, read_len, haplotype_len);
                double y_val = temp1 * Qd[i - 1] + temp2 * Qg[i - 1];
                m_set(temp1 * Qd[i-1] + temp2 * Qg[i-1], i, j, Y, read_len, haplotype_len);


                // Save values to the 2D arrays
                M_matrix[i][j] = m_val;
                X_matrix[i][j] = x_val;
                Y_matrix[i][j] = y_val;
            }
            // antidiags_print(M, min_len);

            //compute likelihood
            if (i == read_len){
                // printf("LIKELYHOOD i:%d, j:%d\n",i,j);
                m_get(&temp1, i, j, M, read_len, haplotype_len);
                m_get(&temp2, i, j, X, read_len, haplotype_len);
                // printf("temp1: %e, temp2: %e\n",temp1, temp2);
                *likelihood += temp1 + temp2;
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
    *likelihood = log10(*likelihood) - log10(DBL_MAX/16);

    // // Print the computed matrices
    // printf("Matrix M:\n");
    // for (i = 0; i <= read_len; i++) {
    //     for (j = 0; j <= haplotype_len; j++) {
    //         printf("%e ", M_matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("Matrix X:\n");
    // for (i = 0; i <= read_len; i++) {
    //     for (j = 0; j <= haplotype_len; j++) {
    //         printf("%e ", X_matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("Matrix Y:\n");
    // for (i = 0; i <= read_len; i++) {
    //     for (j = 0; j <= haplotype_len; j++) {
    //         printf("%e ", Y_matrix[i][j]);
    //     }
    //     printf("\n");
    // }

}

// //final likelihood
// double likelihood(double **M, double **X, int read_len, int haplotype_len) {
//     double l = 0;

//     for (int i = 1; i <= haplotype_len; i++) {
//         l += (M[read_len][i] + X[read_len][i]);
//     }
//     return (double)(log10(l) - log10(DBL_MAX/16));
// }

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


int main(int argc, const char *argv[]) {
    //variables:
    FILE *file_stream_h;
    FILE *file_stream_r;
    FILE *file_stream_out;

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
    double *antidiags = NULL;
    double *M = NULL;
    double *X = NULL;
    double *Y = NULL;

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
        char **haplotypes = (char**)malloc(num_haplotypes * sizeof(char*));
        if (haplotypes == NULL) {
            fprintf(stderr, "Memory allocation failed for haplotypes array\n");
            return EXIT_FAILURE;
        }

        //put the stream to the first haplotype
        for (i = 0; i < num_read; i++) {
            fgets(line, sizeof(line), file_stream_h);
        }
        //get all the haplotypes
        for (i = 0; i < num_haplotypes; i++) {
            if (fgets(line, sizeof(line), file_stream_h) == NULL) {
                fprintf(stderr, "Error reading haplotypes.\n");
                cleanup(antidiags, read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            haplotypes[i] = (char*)malloc((sizeof(line) + 1) * sizeof(char));
            if (haplotypes[i] == NULL) {
                fprintf(stderr, "Error allocating memory for haplotype %d.\n", i);
                cleanup(antidiags, read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }
            strncpy(haplotypes[i], line, sizeof(line));
        }


        //read all the reads and execute pairHMM for each couple (read, haplotype)
        for (i = 0; i < num_read; i++) {
            if (fgets(line, sizeof(line), file_stream_r) == NULL) {
                fprintf(stderr, "Error reading reads.\n");
                cleanup(antidiags, read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            len_read = (strlen(line) - 4) / 5;

            //allocation vectors
            Qr = (double *)malloc(len_read * sizeof(double));
            Qi = (double *)malloc(len_read * sizeof(double));
            Qd = (double *)malloc(len_read * sizeof(double));
            Qg = (double *)malloc(len_read * sizeof(double));
            read = (char *)malloc((len_read + 1) * sizeof(char));
            if (Qr == NULL || Qi == NULL || Qd == NULL || Qg == NULL) {
                fprintf(stderr, "Error allocating memory for quality scores.\n");
                cleanup(antidiags, read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);
                return EXIT_FAILURE;
            }

            //set values to quality arrays
            partition_read(line, read, Qr, Qi, Qd, Qg, len_read);
            //allocation matrices
            // M = (double **)malloc((len_read + 1) * sizeof(double *));
            // X = (double **)malloc((len_read + 1) * sizeof(double *));
            // Y = (double **)malloc((len_read + 1) * sizeof(double *));

            //call the algo
            for (j = 0; j < num_haplotypes; j++) {
                len_haplotype = strlen(haplotypes[j]);
                min_len = len_read < len_haplotype ? len_read : len_haplotype;
                //allocate an array to fit three antidiag of max size min_len+1 for the three matrixes
                antidiags_size = (min_len+1) * 3 * 3;
                nBytes = antidiags_size * sizeof(double);
                antidiags = (double *) malloc (nBytes);
                if (antidiags == NULL) {
                    fprintf(stderr, "Error allocating memory for matrices.\n");
                    cleanup(antidiags, read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);
                    return EXIT_FAILURE;
                }
                memset(antidiags, 0, nBytes);
                M = antidiags;
                X = antidiags + (min_len+1)*3;
                Y = antidiags + 2*(min_len+1)*3;

                //call the algo
                pairHMM(lh, M, X, Y, read, haplotypes[j], len_read, len_haplotype, Qr, Qi, Qd, Qg);
                printf("%f\n", *lh);
                //print likelihood
                fprintf(file_stream_out, "%f\n", *lh);

                //free the memory
                free(antidiags);
                antidiags = NULL;
                M = NULL;
                X = NULL;
                Y = NULL;
            }

            free(Qr);
            free(Qi);
            free(Qd);
            free(Qg);
            free(read);
            Qr = NULL;
            Qi = NULL;
            Qd = NULL;
            Qg = NULL;
            read = NULL;
        }

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

    return EXIT_SUCCESS;
}
