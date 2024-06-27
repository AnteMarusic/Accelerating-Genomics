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
void pairHMM(double **M, double **X, double **Y, char *R, char *H, int read_len, int haplotype_len, double *Qr, double *Qi, double *Qd, double *Qg) {
    //initialization:
    double value = DBL_MAX/16 / (double)haplotype_len;
    for (int i = 0; i <= haplotype_len; i++) {
        Y[0][i] = value;
    }
    
    //algo:
    for (int i = 1; i <= read_len; i++) {
        for (int j = 1; j <= haplotype_len; j++) {
            M[i][j] = p(R[i-1], H[j-1], Qr[i-1]) * (mm(Qi[i-1], Qd[i-1]) * M[i-1][j-1] + (1 - Qg[i-1]) * (X[i-1][j-1] + Y[i-1][j-1]));
            X[i][j] = M[i-1][j] * Qi[i-1] + X[i-1][j] * Qg[i-1];
            Y[i][j] = M[i][j-1] * Qd[i-1] + Y[i][j-1] * Qg[i-1];
        }
    }   
}

//final likelihood
double likelihood(double **M, double **X, int read_len, int haplotype_len) {
    double l = 0;

    for (int i = 1; i <= haplotype_len; i++) {
        l += (M[read_len][i] + X[read_len][i]);
    }
    return (double)(log10(l) - log10(DBL_MAX/16));
}

void cleanup(double **M, double **X, double **Y, int len_read, double *Qr, double *Qi, double *Qd, double *Qg, char **haplotypes, int num_haplotypes) {
    if (M != NULL) {
        for (int i = 0; i <= len_read; i++) {
            if (M[i] != NULL) free(M[i]);
        }
        free(M);
    }
    if (X != NULL) {
        for (int i = 0; i <= len_read; i++) {
            if (X[i] != NULL) free(X[i]);
        }
        free(X);
    }
    if (Y != NULL) {
        for (int i = 0; i <= len_read; i++) {
            if (Y[i] != NULL) free(Y[i]);
        }
        free(Y);
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
    int i, j, k;
    int iteration = 1;

    char read[MAX_READ_LEN+1]; 
    char line[MAX_READ_LEN*5+1];
    
    double lh = 0;
    double **M = NULL;
    double **X = NULL;
    double **Y = NULL;
    double *Qr = NULL;
    double *Qi = NULL;
    double *Qd = NULL;
    double *Qg = NULL;

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
                cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            haplotypes[i] = (char*)malloc((sizeof(line) + 1) * sizeof(char));
            if (haplotypes[i] == NULL) {
                fprintf(stderr, "Error allocating memory for haplotype %d.\n", i);
                cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                return EXIT_FAILURE;
            }
            strncpy(haplotypes[i], line, sizeof(line));
        }
        

        //read all the reads and execute pairHMM for each couple (read, haplotype)
        for (i = 0; i < num_read; i++) {
            if (fgets(line, sizeof(line), file_stream_r) == NULL) {
                fprintf(stderr, "Error reading reads.\n");
                cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                return EXIT_FAILURE;
            }
            line[strcspn(line, "\n")] = '\0'; // Remove the newline character
            len_read = (strlen(line) - 4) / 5;

            //allocation vectors
            Qr = (double *)malloc(len_read * sizeof(double));
            Qi = (double *)malloc(len_read * sizeof(double));
            Qd = (double *)malloc(len_read * sizeof(double));
            Qg = (double *)malloc(len_read * sizeof(double));
            if (Qr == NULL || Qi == NULL || Qd == NULL || Qg == NULL) {
                fprintf(stderr, "Error allocating memory for quality scores.\n");
                cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                return EXIT_FAILURE;
            }

            //set values to quality arrays
            partition_read(line, read, Qr, Qi, Qd, Qg, len_read);
            //allocation matrices
            M = (double **)malloc((len_read + 1) * sizeof(double *));
            X = (double **)malloc((len_read + 1) * sizeof(double *));
            Y = (double **)malloc((len_read + 1) * sizeof(double *));
            if (M == NULL || X == NULL || Y == NULL) {
                fprintf(stderr, "Error allocating memory for matrices.\n");
                cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                return EXIT_FAILURE;
            }
            
            //call the algo
            for (j = 0; j < num_haplotypes; j++) {
                len_haplotype = strlen(haplotypes[j]);
                for (k = 0; k <= len_read; k++) {
                    M[k] = (double *)malloc((len_haplotype + 1) * sizeof(double));
                    X[k] = (double *)malloc((len_haplotype + 1) * sizeof(double));
                    Y[k] = (double *)malloc((len_haplotype + 1) * sizeof(double));
                    if (M[k] == NULL || X[k] == NULL || Y[k] == NULL) {
                        fprintf(stderr, "Error allocating memory for matrix rows.\n");
                        cleanup(M, X, Y, len_read, Qr, Qi, Qd, Qg, haplotypes, num_haplotypes);     
                        return EXIT_FAILURE;
                    }
                    memset(M[k], 0, (len_haplotype + 1) * sizeof(double));
                    memset(X[k], 0, (len_haplotype + 1) * sizeof(double));
                    memset(Y[k], 0, (len_haplotype + 1) * sizeof(double));
                }
                pairHMM(M, X, Y, read, haplotypes[j], len_read, len_haplotype, Qr, Qi, Qd, Qg);
                //compute likelihood
                lh = likelihood(M, X, len_read, len_haplotype);
                fprintf(file_stream_out, "%f\n", lh);

                // printf("Matrix M:\n");
                // for (i = 0; i <= len_read; i++) {
                //     for (j = 0; j <= len_haplotype; j++) {
                //         printf("%e ", M[i][j]);
                //     }
                //     printf("\n");
                // }
                // printf("Matrix X:\n");
                // for (i = 0; i <= len_read; i++) {
                //     for (j = 0; j <= len_haplotype; j++) {
                //         printf("%e ", X[i][j]);
                //     }
                //     printf("\n");
                // }
                // printf("Matrix Y:\n");
                // for (i = 0; i <= len_read; i++) {
                //     for (j = 0; j <= len_haplotype; j++) {
                //         printf("%e ", Y[i][j]);
                //     }
                //     printf("\n");
                // }

                //free the memory
                for (k = 0; k <= len_read; k++) {
                    free(M[k]);
                    free(X[k]);
                    free(Y[k]);
                    M[k] = NULL;
                    X[k] = NULL;
                    Y[k] = NULL;
                }
            }
            
            free(M);
            free(X);
            free(Y);
            free(Qr);
            free(Qi);
            free(Qd);
            free(Qg);
            M = NULL;
            X = NULL;
            Y = NULL;
            Qr = NULL;
            Qi = NULL;
            Qd = NULL;
            Qg = NULL;
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
