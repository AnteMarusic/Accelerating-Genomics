#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Include the string.h header for strcspn function
#include <math.h>
#include <limits.h>
#include <float.h>

//#undef DEBUG
#define DEBUG
#ifdef DEBUG
#define PRINT printf
#else
#define PRINT // macros
#endif

#define OFFSET 33  // Typical offset for quality scores in FASTQ format
#define M 1
#define X 2
#define Y 3
#define MAX_LINE_LENGTH 10000

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

int *convert_quality_string(const char *quality_string, size_t length) {
    int *quality_array = (int *)malloc(length * sizeof(int));
    for (size_t i = 0; i < length; i++) {
        quality_array[i] = (int)quality_string[i] - OFFSET;
    }
    return quality_array;
}

double epsylon(int read_quality) {
    return pow(10, -read_quality/10);
}
double mx(int insert_quality) {
    return pow(10, -insert_quality/10);
}
double my (int delete_quality) {
    return pow(10, -delete_quality/10);
}
double p(char read, char haplotype, int read_quality){
    return read == haplotype? 1-epsylon(read_quality)*read_quality : epsylon(read_quality)*read_quality/3;
}
double mm(int insert_quality, int delete_quality) {
    return 1 - (mx(insert_quality) + my(delete_quality));
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
int m_get(double *val, int iy, int ix, int mat, double *antidiags, int cur_antid_num, int nx, int ny) {
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
        case M :
            offset_mat = 0;
            break;
        case X :
            offset_mat = antid_dim*3;
            break;
        case Y :
            offset_mat = antid_dim*2*3;
            break;
        default :
            printf("[m_get] mat doesn't exist\n");
    }
    cur_antid_num = ix + iy;
    offset_antid = cur_antid_num % 3;//three is the num of antid saved at the same time

    if(cur_antid_num > nx<ny ? ny-1:nx-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }
    *val = antidiags[index];
    //PRINT("[m_get] getting val: %f, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", *val, iy, ix, mat, index);
    PRINT("[m_get] getting val: %.2e, at iy:%d, ix:%d, on matrix: %d, at index: %d\n", *val, iy, ix, mat, index);
    return 1;
}

//returns positive values in case of success and the value is in val
//returns negative number in case of failure
int m_set(double *val, int iy, int ix, int mat, double *antidiags, int cur_antid_num, int nx, int ny) {
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
        case M :
            offset_mat = 0;
            break;
        case X :
            offset_mat = antid_dim*3;
            break;
        case Y :
            offset_mat = antid_dim*2*3;
            break;
        default :
            printf("[m_set] mat doesn't exist\n");
    }
    cur_antid_num = ix + iy;
    //three is the num of antid saved at the same time
    offset_antid = cur_antid_num % 3;
    if(cur_antid_num > nx<ny ? ny-1:nx-1){
        index = ny - 1 - iy + offset_antid*antid_dim + offset_mat;
    }
    else {
        index = ix + offset_antid*antid_dim + offset_mat;
    }

    //PRINT("[m_set] setting val: %f, iy:%d, ix:%d, at index: %d on matrix %d\n", *val, iy, ix, index, mat);
    PRINT("[m_set] setting val: %.2e iy:%d, ix:%d, at index: %d on matrix %d\n", *val, iy, ix, index, mat);
    antidiags[index] = *val;
    return 1;
}



int main(int argc, char *argv[]) {
    //read read_quality insert_quality delete_quality
    if (argc != 2) {
        printf("Usage: %s <input_file_path>\n", argv[0]);
        return 1;
    }
    // Open the file for reading
    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening file: %s\n", argv[1]);
        return 1;
    }
    char *input_string;
    char *haplotype;

    // Allocate memory for the strings
    input_string = (char *)malloc(MAX_LINE_LENGTH * sizeof(char));
    haplotype = (char *)malloc(MAX_LINE_LENGTH * sizeof(char));
    if (input_string == NULL || haplotype == NULL) {
        fprintf(stderr, "Error allocating memory\n");
        fclose(fp);
        return 1;
    }

    // Read the first line
    if (fgets(input_string, MAX_LINE_LENGTH, fp) == NULL) {
        fprintf(stderr, "Error reading first line\n");
        free(input_string);
        free(haplotype);
        fclose(fp);
        return 1;
    }

    // Remove the trailing newline character (if present)
    size_t len = strlen(input_string);
    if (len > 0 && input_string[len - 1] == '\n') {
        input_string[len - 1] = '\0';
    }

    // Read the second line
    if (fgets(haplotype, MAX_LINE_LENGTH, fp) == NULL) {
        fprintf(stderr, "Error reading second line\n");
        free(input_string);
        free(haplotype);
        fclose(fp);
        return 1;
    }

    // Remove the trailing newline character (if present) for haplotype
    len = strlen(haplotype);
    if (len > 0 && haplotype[len - 1] == '\n') {
        haplotype[len - 1] = '\0';
    }

    // Close the file
    fclose(fp);

    // Print the read strings (optional)
    printf("Input String: %s\n", input_string);
    printf("Haplotype: %s\n", haplotype);

    //parameters obtained by the sequencing machine
    double xx = 0.1;
    double yy = 0.1;
    double gm = 0.9;

    // Extract parts of the string
    char *read = strtok(input_string, " ");
    char *read_quality_str = strtok(NULL, " ");
    char *insert_quality_str = strtok(NULL, " ");
    char *delete_quality_str = strtok(NULL, " ");

    // Lengths of the quality strings
    size_t read_quality_length = strlen(read_quality_str);
    size_t insert_quality_length = strlen(insert_quality_str);
    size_t delete_quality_length = strlen(delete_quality_str);
    int haplotype_len = strlen(haplotype);
    // Convert quality strings to int arrays
    int *read_quality = convert_quality_string(read_quality_str, read_quality_length);
    int *insert_quality = convert_quality_string(insert_quality_str, insert_quality_length);
    int *delete_quality = convert_quality_string(delete_quality_str, delete_quality_length);
    /*
    // Output results for demonstration
    printf("Read: %s\n", read);
    printf("Read Quality: ");
    for (size_t i = 0; i < read_quality_length; i++) {
        printf("%d ", read_quality[i]);
    }
    printf("\nInsert Quality: ");
    for (size_t i = 0; i < insert_quality_length; i++) {
        printf("%d ", insert_quality[i]);
    }
    printf("\nDelete Quality: ");
    for (size_t i = 0; i < delete_quality_length; i++) {
        printf("%d ", delete_quality[i]);
    }
    printf("\nepsylon: ");
    for (size_t i = 0; i < read_quality_length; i++) {
        printf("%f ", epsylon(read_quality[i]));
    }
    printf("\nmy: ");
    for (size_t i = 0; i < delete_quality_length; i++) {
        printf("%f ", my(delete_quality[i]));
    }
    printf("\n");
    */

    char *sx = haplotype;
    char *sy = read;
    int sx_len = strlen(sx);
    int sy_len = strlen(sy);
    printf("[main] sx: %s, xy: %s, of size %d, %d respectively\n", sx, sy, sx_len, sy_len);
    int nx = sx_len+1;
    int ny = sy_len+1;
    printf("[main] nx: %d, ny: %d\n", nx, ny);
    //allocate an array to fit three antidiag of max size min_len for the three matrixes
    int antidiags_size = nx<ny ? nx * 3 * 3 : ny * 3 * 3;
    int nBytes = antidiags_size * sizeof(double);
    double *antidiags = (double *) malloc (nBytes);
    int antid_num = nx + ny - 1;
    int yr_phase = nx<ny ? nx - 1 : ny - 1;//yellow is the phase where the cur_antid_dim less then the max possible and is growing
    //with each iteration. red is the same but cur_antid_dim is decreasing. red and yellow have the same size.
    int o_phase = nx<ny ? ny - nx + 1: nx - ny + 1;//orange phase is the one where the cur_antid_dim is the max possible
    int err;
    double val;
    int iy = 0, ix = 0;
    int cur_antid_dim;
    int gen_iy = 0;
    int gen_ix = 0;
    double temp1,temp2,temp3, result = 0;
    for(int cur_antid_num = 0; cur_antid_num  < antid_num; cur_antid_num  ++) {

        if (cur_antid_num <= yr_phase) { //at index yr_phase starts the o_phase
            printf("[main] growing phase\n");
            cur_antid_dim = cur_antid_num + 1;
        } else if (cur_antid_num >= yr_phase + o_phase) {
            printf("[main] decreasing phase\n");
            cur_antid_dim = cur_antid_dim - 1;
        } else {
            printf("[main] const phase\n");
        }
        printf("[main] cur_antid_num: %d, cur_antid_dim: %d\n", cur_antid_num, cur_antid_dim);
        if (cur_antid_num >= ny) {
            gen_ix++;
        }
        iy = gen_iy;
        ix = gen_ix;
        for (int j = 0; j < cur_antid_dim; j++) {
            printf("iy: %d, ix: %d\n", iy, ix);
            //set the first row
            if (iy == 0) {
                val = 0;
                err = m_set(&val, iy, ix, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = 0;
                err = m_set(&val, iy, ix, X, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = pow(2, 1020) / haplotype_len;
                err = m_set(&val, iy, ix, Y, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
            } else if (ix == 0) {
                val = 0;
                err = m_set(&val, iy, ix, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = 0;
                err = m_set(&val, iy, ix, X, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = 0;
                err = m_set(&val, iy, ix, Y, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
            }
                //dynamic programming
            else {
                err = m_get(&temp1, iy - 1, ix - 1, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                err = m_get(&temp2, iy - 1, ix - 1, X, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                err = m_get(&temp3, iy - 1, ix - 1, Y, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                PRINT("sy[iy-1]:%c, sx[ix-1]:%c\n", sy[iy - 1], sx[ix - 1]);
                val = p(read[iy], haplotype[ix], read_quality[iy]) * (
                        mm(insert_quality[iy], delete_quality[iy]) * temp1 + gm * (temp2 + temp3)
                );
                err = m_set(&val, iy, ix, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;

                err = m_get(&temp1, iy, ix - 1, X, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                err = m_get(&temp2, iy, ix - 1, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = xx * temp1 + mx(insert_quality[iy]) * temp2;
                err = m_set(&val, iy, ix, X, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;

                err = m_get(&temp1, iy - 1, ix, Y, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                err = m_get(&temp2, iy - 1, ix, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                val = yy * temp1 + mx(insert_quality[iy]) * temp2;
                err = m_set(&val, iy, ix, Y, antidiags, cur_antid_num, nx, ny);
            }
            if (iy == ny - 1) {
                err = m_get(&temp1, iy, ix, M, antidiags, cur_antid_num, nx, ny);
                if (err < 1) return -1;
                err = m_get(&temp2, iy, ix, X, antidiags, cur_antid_num, nx, ny); //se metto Y smette di darmi NAN
                if (err < 1) return -1;
                result += temp1 + temp2;
            }
            //update the matrix indexes
            iy--;
            ix++;
        }
        //move the starting point in case we are in the growing part of the matrix
        if (gen_iy < ny - 1) {
            gen_iy++;
        }
        //printf("cur_antd_num: %d\n", cur_antid_num);
        //a_print(antidiags, antidiags_size, nx);
    }
    printf("result: %f\n", result);

    result = log10(result)-log10(pow(2,1020));
    //result = log10(result)-log10(DBL_MAX/16); Beatrice usa questo

    printf("result: %f\n", result);
    // Free allocated memory
    free(read_quality);
    free(insert_quality);
    free(delete_quality);
    free(input_string);
    free(haplotype);


    return 0;
}
