#ifndef HEADER_FILE_H
#define HEADER_FILE_H
#define TRUE 1
#define FALSE 0

#undef DEBUG
//define DEBUG

#ifdef DEBUG
#define PRINT printf
#else
#define PRINT // macros
#endif

#endif

