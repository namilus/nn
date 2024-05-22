#include <stdlib.h>
#include <emacs-module.h>
#include "utils.h"
/* Declare mandatory GPL symbol.  */
int plugin_is_GPL_compatible;

static emacs_value shape(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Assumes the input matrix is non-ragged */
  emacs_value vector = env->intern(env, "vector");

  size_t shape[2];
  extract_matrix_shape(shape, env, args[0]);

  emacs_value result[2] = {env->make_integer(env, shape[0]), env->make_integer(env, shape[1])};
  return env->funcall(env, vector, 2, result);

}


static emacs_value random_matrix(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Create a matrix of uniform random values between 1 and
     -1. Arguments are nrows and ncols */
  emacs_value vector = env->intern(env, "vector");
  size_t nrows = (size_t) env->extract_integer(env, args[0]);
  size_t ncols = (size_t) env->extract_integer(env, args[1]);

  emacs_value rows[nrows];
  emacs_value row_i[ncols];
  for(size_t i=0; i<nrows; i++){
    for(size_t j=0; j<ncols; j++){
      bool sign = rand() % 2 == 0;
      double value = (double) rand() / RAND_MAX;
      double ij = sign ? value : -1 * value;
      row_i[j] = env->make_float(env, ij);
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }
  return env->funcall(env, vector, nrows, rows);
}


static emacs_value transpose(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Assumes the input matrix is non-ragged */
  emacs_value vector = env->intern(env, "vector");

  size_t shape[2];
  extract_matrix_shape(shape, env, args[0]);
  
  double matrix[shape[0]][shape[1]];
  extract_matrix((double*) matrix, env, args[0], shape);

  emacs_value rows[shape[1]];

  emacs_value row_i[shape[0]];

  for(size_t i=0; i<shape[1]; i++){
    for(size_t j=0; j<shape[0]; j++){
      row_i[j] = env->make_float(env, matrix[j][i]);
    }
    rows[i] = env->funcall(env, vector, shape[0], row_i);
  }

  return env->funcall(env, vector, shape[1], rows);
                                      
}

static emacs_value add(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Variadic matrix addition. Works by first "loading" the first
     matrix, and then subsequently adding the corresponding values to
     each other from the further argument matrices. */
  emacs_value vector = env->intern(env, "vector");

  /* Load the first matrix and get its shape */
  size_t m1_shape[2];
  extract_matrix_shape(m1_shape, env, args[0]);

  double m1[m1_shape[0]][m1_shape[1]];
  extract_matrix((double*) m1, env, args[0], m1_shape);

  /* Array that will hold the final result */
  double result[m1_shape[0]][m1_shape[1]];

  /* Variable that will hold the m-th matrix that we're adding to
     `result' */
  double matrix_m[m1_shape[0]][m1_shape[1]];

  /* Loop through the remaining nargs - 1 argument matrices add sum
     them up cumulatively */
  for(int m=1; m < nargs; m++){

    extract_matrix((double*) matrix_m, env, args[m], m1_shape);

    for(int i=0; i < m1_shape[0]; i++){
      for(int j=0; j < m1_shape[1]; j++){

        double ij;
        if(m == 1){
          ij = m1[i][j] + matrix_m[i][j];
        }
        else {
          ij = result[i][j] + matrix_m[i][j];
        }
        result[i][j] = ij;
      }
    }
  }

  /* Convert result array into an emacs vector of vectors */
  emacs_value rows[m1_shape[0]];

  for(int i=0; i < m1_shape[0]; i++){
    emacs_value row_i[m1_shape[1]];
    for(int j=0; j < m1_shape[1]; j++){
      row_i[j] = env->make_float(env, result[i][j]);
    }
    rows[i] = env->funcall(env, vector, m1_shape[1], row_i);
  }

  return env->funcall(env, vector, m1_shape[0], rows);
}


static emacs_value subtract(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Variadic matrix subtraction. Works by first "loading" the first
     matrix, and then subsequently subtracting the corresponding
     values to each other from the further argument matrices. */
  emacs_value vector = env->intern(env, "vector");

  size_t m1_shape[2];
  extract_matrix_shape(m1_shape, env, args[0]);


  double m1[m1_shape[0]][m1_shape[1]];
  extract_matrix((double*) m1, env, args[0], m1_shape);
    
  double result[m1_shape[0]][m1_shape[1]];

  double matrix_m[m1_shape[0]][m1_shape[1]];
  for(int m=1; m < nargs; m++){

    extract_matrix((double*) matrix_m, env, args[m], m1_shape);

    for(int i=0; i < m1_shape[0]; i++){
      for(int j=0; j < m1_shape[1]; j++){

        double ij;
        if(m == 1){
          ij = m1[i][j] - matrix_m[i][j];
        }
        else {
          ij = result[i][j] - matrix_m[i][j];
        }

        result[i][j] = ij;
      }
    }
  }

  /* Convert result array into an emacs vector of vectors */
  emacs_value rows[m1_shape[0]];

  for(int i=0; i < m1_shape[0]; i++){
    emacs_value row_i[m1_shape[1]];
    for(int j=0; j < m1_shape[1]; j++){
      row_i[j] = env->make_float(env, result[i][j]);
    }
    rows[i] = env->funcall(env, vector, m1_shape[1], row_i);
  }

  return env->funcall(env, vector, m1_shape[0], rows);
}


static emacs_value hadamard(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Variadic matrix element-wise multiplication. Works by first
     "loading" the first matrix, and then subsequently multiplying the
     corresponding values to each other from the further argument
     matrices. */
  emacs_value vector = env->intern(env, "vector");

  size_t m1_shape[2];
  extract_matrix_shape(m1_shape, env, args[0]);


  double m1[m1_shape[0]][m1_shape[1]];
  extract_matrix((double*) m1, env, args[0], m1_shape);
    
  double result[m1_shape[0]][m1_shape[1]];

  double matrix_m[m1_shape[0]][m1_shape[1]];
  for(int m=1; m < nargs; m++){

    extract_matrix((double*) matrix_m, env, args[m], m1_shape);

    for(int i=0; i < m1_shape[0]; i++){
      for(int j=0; j < m1_shape[1]; j++){

        double ij;
        if(m == 1){
          ij = m1[i][j] * matrix_m[i][j];
        }
        else {
          ij = result[i][j] * matrix_m[i][j];
        }

        result[i][j] = ij;
      }
    }
  }

  /* Convert result array into an emacs vector of vectors */
  emacs_value rows[m1_shape[0]];

  for(int i=0; i < m1_shape[0]; i++){
    emacs_value row_i[m1_shape[1]];
    for(int j=0; j < m1_shape[1]; j++){
      row_i[j] = env->make_float(env, result[i][j]);
    }
    rows[i] = env->funcall(env, vector, m1_shape[1], row_i);
  }

  return env->funcall(env, vector, m1_shape[0], rows);
}


static emacs_value scalar_multiply(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Multiply matrix elements by a given scalar */
  emacs_value vector = env->intern(env, "vector");

  double scalar = extract_double(env, args[0]);

  size_t shape[2];
  extract_matrix_shape(shape, env, args[1]);
  
  double matrix[shape[0]][shape[1]];
  extract_matrix((double*) matrix, env, args[1], shape);

  emacs_value rows[shape[0]];

  emacs_value row_i[shape[1]];

  for(size_t i=0; i<shape[0]; i++){
    for(size_t j=0; j<shape[1]; j++){
      row_i[j] = env->make_float(env, scalar * matrix[i][j]);
    }
    rows[i] = env->funcall(env, vector, shape[1], row_i);
  }

  return env->funcall(env, vector, shape[0], rows);
                                      
}

static emacs_value matmul(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Variadic matrix multiplication. Works by first "loading" the
     first two argument matrices, multiplying them, and then
     subsequently multiplying the result with the next argument
     matrix. This is done until we've exhausted the argument list */
  emacs_value vector = env->intern(env, "vector");


  /* Extract the first matrix and its shape */
  size_t m1_shape[2];
  extract_matrix_shape(m1_shape, env, args[0]);

  double m1[m1_shape[0]][m1_shape[1]];
  extract_matrix((double*) m1, env, args[0], m1_shape);

  /* Declare pointers to the result of the multiplication of the
     previous 2 argument matrices. Initially, this just points to the
     first matrix */
  size_t* m_prev_shape = m1_shape;
  double* m_prev = (double*) m1;


  /* Will hold the result of the multiplication of the previous result
     matrix  */
  double* result_m;
  result_m = NULL;

  for(int m=1; m < nargs; m++){

    /* Get the shape and load the m-th matrix */
    size_t matrix_m_shape[2];
    extract_matrix_shape(matrix_m_shape, env, args[m]);
    double matrix_m[matrix_m_shape[0]][matrix_m_shape[1]];
    extract_matrix((double*) matrix_m, env, args[m], matrix_m_shape);

    /* Construct the shape of the resulting matrix. Generally,
       multiplying a m x n matrix with a n x p matrix results in a m x
       p matrix */
    size_t result_m_shape[] = {m_prev_shape[0], matrix_m_shape[1]};
        

    /* We want the final result to persist outside the loop so we can
       construct an emacs vector. Therefore we dynamically allocate
       and free at the end */
    result_m = malloc(sizeof(double) * result_m_shape[0] * result_m_shape[1]);

    size_t shared = matrix_m_shape[0];
    
    for(int i=0; i < m_prev_shape[0]; i++){
      for(int j=0; j < matrix_m_shape[1]; j++){
        double ij = 0;
        for(int k=0; k < shared; k++){
          ij += *(m_prev + (m_prev_shape[1]*i) + k) * matrix_m[k][j];
        }
        result_m[result_m_shape[1]*i + j] = ij;
      }
    }
    /* Free each intermediate result matrix. If m == 1, `m_prev' is
       just the argument m1 matrix, and so we don't call free on
       that. */
    if (m > 1) {
      free(m_prev);
    }
    m_prev_shape = result_m_shape;
    m_prev = result_m;

  }

  /* Convert result array into an emacs vector of vectors */
  emacs_value rows[m1_shape[0]];

  for(int i=0; i < m_prev_shape[0]; i++){
    emacs_value row_i[m_prev_shape[1]];
    for(int j=0; j < m_prev_shape[1]; j++){
      row_i[j] = env->make_float(env, m_prev[m_prev_shape[1]*i + j]);
    }
    rows[i] = env->funcall(env, vector, m_prev_shape[1], row_i);
  }

  /* We're done with it now, so we can free it */
  free(result_m);

  return env->funcall(env, vector, m_prev_shape[0], rows);
}

static void bind_function(emacs_env *env, const char *name, emacs_value Sfun){
  /* Set the function cell of the symbol named NAME to SFUN using
     the 'fset' function.  */

  /* Convert the strings to symbols by interning them */
  emacs_value Qfset = env->intern (env, "fset");
  emacs_value Qsym = env->intern (env, name);

  /* Prepare the arguments array */
  emacs_value args[] = { Qsym, Sfun };

  /* Make the call (2 == nb of arguments) */
  env->funcall (env, Qfset, 2, args);
}

static void provide(emacs_env *env, const char *feature){
  /* call 'provide' with FEATURE converted to a symbol */

  emacs_value Qfeat = env->intern (env, feature);
  emacs_value Qprovide = env->intern (env, "provide");
  emacs_value args[] = { Qfeat };
  env->funcall (env, Qprovide, 1, args);
}



int emacs_module_init(struct emacs_runtime *ert){
  /* Compat checks */
  if (ert->size < sizeof (*ert))
    return 1;

  emacs_env *env = ert->get_environment (ert);
  if (env->size < sizeof (*env))
    return 2;

  emacs_value e_matrix_random = env->make_function(env,
                                                   2,
                                                   2,
                                                   random_matrix,
                                                   "Return a random matrix of provided shape.",
                                                   NULL);

  bind_function (env, "matrix-random", e_matrix_random);



  emacs_value e_matrix_subtract = env->make_function(env,
                                                     2,
                                                     -2,
                                                     subtract,
                                                     "Subtract matrices in argument order and return the result.",
                                                     NULL);

  bind_function (env, "matrix-subtract", e_matrix_subtract);


  emacs_value e_matrix_hadamard = env->make_function(env,
                                                     2,
                                                     -2,
                                                     hadamard,
                                                     "Multiply matrices element-wise in argument order and return the result.",
                                                     NULL);

  bind_function (env, "matrix-hadamard", e_matrix_hadamard);
  

  emacs_value e_matrix_add = env->make_function(env,
                                                 2,
                                                 -2,
                                                 add,
                                                 "Return the sum of 2 or more matrices.",
                                                 NULL);
  bind_function(env, "matrix-add", e_matrix_add);


  emacs_value e_matrix_transpose = env->make_function(env,
                                                      1,
                                                      1,
                                                      transpose,
                                                      "Return the transpose of the argument matrix",
                                                      NULL);
  bind_function(env, "matrix-transpose", e_matrix_transpose);


  emacs_value e_matrix_scalar_multiply = env->make_function(env,
                                                      2,
                                                      2,
                                                      scalar_multiply,
                                                      "Multiply scalar and argument matrix",
                                                      NULL);
  bind_function(env, "matrix-scalar-mul", e_matrix_scalar_multiply);  


  emacs_value e_matrix_shape = env->make_function(env,
                                                  1,
                                                  1,
                                                  shape,
                                                  "Return the shape of the argument matrix. Assumes that all the rows have the same number of elements.",
                                                  NULL);
  bind_function(env, "matrix-shape", e_matrix_shape);



  emacs_value e_matrix_matmul = env->make_function(env,
                                                   1,
                                                   -2,
                                                   matmul,
                                                   "Multiply 2 or more matrices in argument order and return the result.",
                                                   NULL);
  bind_function(env, "matrix-matmul", e_matrix_matmul);




  provide(env, "matrix");

  /* loaded successfully */
  return 0;
}
