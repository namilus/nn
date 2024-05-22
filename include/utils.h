#ifndef UTILS_H
#define UTILS_H

#include <emacs-module.h>

/* Utility functions for converting emacs arguments to their c
   equivalents */
double extract_double(emacs_env *env, emacs_value arg);
void extract_matrix_shape(size_t* shape, emacs_env* env, emacs_value arg);
void extract_matrix(double* matrix, emacs_env *env, emacs_value arg, size_t shape[]);


#endif
