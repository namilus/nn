#include "utils.h"

double extract_double(emacs_env *env, emacs_value arg){
  /* Given an emacs value `arg' that represents either an integer or a
     float, return the corresponding double. If it is neither, the
     value 0.0 is returned */
  emacs_value type = env->type_of(env, arg);
  double result = 0.0;
  /* If the value is an integer, cast it to a double */
  if (env->eq(env, type, env->intern(env, "integer"))) {
      result = (double) env->extract_integer(env, arg);
  }
  else if (env->eq(env, type, env->intern(env, "float"))){
      result = env->extract_float(env, arg);
  }
  return result;
}

void extract_matrix_shape(size_t* shape, emacs_env* env, emacs_value arg){
  /* Makes the assumption that the vector represented by `arg' is
     non-ragged e.g. all the rows have the same number of elements */

  size_t nrows = env->vec_size(env, arg);

  emacs_value row1 = env->vec_get(env, arg, 0);
  size_t ncols = env->vec_size(env, row1);

  shape[0] = nrows;
  shape[1] = ncols;
}

void extract_matrix(double* matrix, emacs_env *env, emacs_value arg, size_t shape[]){
  /* Assumes `arg' is a non-ragged vector of vectors e.g a vector
     containing vectors all of the same length */
  size_t nrows = shape[0];
  size_t ncols = shape[1];

  emacs_value row, ij;

  for(int i=0; i < nrows; i++){
    row = env->vec_get(env, arg, i);
    for(int j=0; j < ncols; j++){
      ij = env->vec_get(env, row, j);
      *(matrix + (ncols*i + j)) = extract_double(env, ij);
    }
  }
}


