#include <stdlib.h>
#include <emacs-module.h>
#include "utils.h"
#include "math.h"

/* Declare mandatory GPL symbol.  */
int plugin_is_GPL_compatible;


static emacs_value softmax(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Apply softmax along a given axis   */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t axis = 1;              /* Default axis to sum over is 1 */
  if(nargs == 2) {
    axis = (size_t) env->extract_integer(env, args[1]);
  }

  /*  \sum_i exp(a_ij) */
  double axis_exp_sums[m_shape[axis]];

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  emacs_value rows[nrows];
  emacs_value row_i[ncols]; 

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){
      double exp_ij = exp(m[i][j]);

      /* calculate the denominator, depending on what axis we need to
         sum across */

      if(axis == 0 && j == 0) {
        /* Only need to calculate the sum across the row once, so do
           it just once when j = 0 */
        double exp_sum_i = 0;
        for(size_t k = 0; k < ncols; k++){
          exp_sum_i += exp(m[i][k]);
        }
        axis_exp_sums[i] = exp_sum_i;
      }
      else if(axis == 1 && i == 0) {
        /* Only need to calculate the sum along the column once, so do
           it just once when i = 0 */
        double exp_sum_j = 0;
        for(size_t k = 0; k < nrows; k++){
          exp_sum_j += exp(m[k][j]);
        }
        axis_exp_sums[j] = exp_sum_j;
      }

      /* Calculate the softmaxed value for ij, depending on the
         axis */
      if(axis == 0){
        row_i[j] = env->make_float(env, exp_ij / axis_exp_sums[i]);
      }
      else if(axis == 1){
        row_i[j] = env->make_float(env, exp_ij / axis_exp_sums[j]);
      }
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }

  return env->funcall(env, vector, nrows, rows);
}


static emacs_value reduce_sum(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Calculate the sum along a given axis   */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t axis = (size_t) env->extract_integer(env, args[1]);

  /*  \sum_i exp(a_ij) */
  emacs_value axis_sums[m_shape[axis]];

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){

      if(axis == 0 && j == 0) {
        /* Only need to calculate the sum across the row once, so do
           it just once when j = 0 */
        double sum_i = 0;
        for(size_t k = 0; k < ncols; k++){
          sum_i += m[i][k];
        }
        axis_sums[i] = env->make_float(env, sum_i);
      }
      else if(axis == 1 && i == 0) {
        /* Only need to calculate the sum along the column once, so do
           it just once when i = 0 */
        double sum_j = 0;
        for(size_t k = 0; k < nrows; k++){
          sum_j += m[k][j];
        }
        axis_sums[j] = env->make_float(env, sum_j);
      }
    }
  }
  if (axis == 0){
    emacs_value rows[nrows];
    for(size_t i=0; i < nrows; i++){
      emacs_value row_i_arg[] = {axis_sums[i]};
      rows[i] = env->funcall(env, vector, 1, row_i_arg);
    }
    return env->funcall(env, vector, nrows, rows);
    
  }
  else if(axis == 1){
    emacs_value final_row = env->funcall(env, vector, ncols, axis_sums);
    emacs_value final_row_arg[] = {final_row};
    return env->funcall(env, vector, 1, final_row_arg);
  }
  return env->make_integer(env, 42); // :?
}

static emacs_value reduce_mean(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Calculate the sum along a given axis   */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t axis = (size_t) env->extract_integer(env, args[1]);

  /*  \sum_i exp(a_ij) */
  emacs_value axis_sums[m_shape[axis]];

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){

      if(axis == 0 && j == 0) {
        /* Only need to calculate the sum across the row once, so do
           it just once when j = 0 */
        double sum_i = 0;
        for(size_t k = 0; k < ncols; k++){
          sum_i += m[i][k];
        }
        axis_sums[i] = env->make_float(env, sum_i / ncols);
      }
      else if(axis == 1 && i == 0) {
        /* Only need to calculate the sum along the column once, so do
           it just once when i = 0 */
        double sum_j = 0;
        for(size_t k = 0; k < nrows; k++){
          sum_j += m[k][j];
        }
        axis_sums[j] = env->make_float(env, sum_j / nrows);
      }
    }
  }
  if (axis == 0){
    emacs_value rows[nrows];
    for(size_t i=0; i < nrows; i++){
      emacs_value row_i_arg[] = {axis_sums[i]};
      rows[i] = env->funcall(env, vector, 1, row_i_arg);
    }
    return env->funcall(env, vector, nrows, rows);
    
  }
  else if(axis == 1){
    emacs_value final_row = env->funcall(env, vector, ncols, axis_sums);
    emacs_value final_row_arg[] = {final_row};
    return env->funcall(env, vector, 1, final_row_arg);
  }
  return env->make_integer(env, 42); // :?
}





static emacs_value relu(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Apply ReLU activation */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  emacs_value rows[nrows];
  emacs_value row_i[ncols]; 

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){
      double relu_ij = m[i][j] > 0 ? m[i][j] : 0.0;
      row_i[j] = env->make_float(env, relu_ij);
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }

  return env->funcall(env, vector, nrows, rows);
}


static emacs_value heaviside(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Apply Heaviside step function */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  emacs_value rows[nrows];
  emacs_value row_i[ncols]; 

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){
      double heaviside_ij = m[i][j] > 0 ? 1.0 : 0.0;
      row_i[j] = env->make_float(env, heaviside_ij);
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }

  return env->funcall(env, vector, nrows, rows);
}

static emacs_value natural_log(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Apply natural logarithm */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  emacs_value rows[nrows];
  emacs_value row_i[ncols]; 

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){
      double log_ij = log(m[i][j]);
      row_i[j] = env->make_float(env, log_ij);
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }

  return env->funcall(env, vector, nrows, rows);
}

static emacs_value identity(emacs_env* env, ptrdiff_t nargs, emacs_value args[], void *data){
  /* Apply identity e.g. do nothing :P   */
  emacs_value vector = env->intern(env, "vector");

  /* Extract the matrix and its shape */
  size_t m_shape[2];
  extract_matrix_shape(m_shape, env, args[0]);

  double m[m_shape[0]][m_shape[1]];
  extract_matrix(m, env, args[0], m_shape);

  size_t nrows = m_shape[0];
  size_t ncols = m_shape[1];

  emacs_value rows[nrows];
  emacs_value row_i[ncols]; 

  for(size_t i=0; i < nrows; i++){
    for(size_t j=0; j < ncols; j++){
      row_i[j] = env->make_float(env, m[i][j]);
    }
    rows[i] = env->funcall(env, vector, ncols, row_i);
  }

  return env->funcall(env, vector, nrows, rows);

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


  emacs_value e_ops_softmax = env->make_function(env,
                                                 1,
                                                 2,
                                                 softmax,
                                                 "Apply softmax along a given axis. Default axis is 1.",
                                                 NULL);
  bind_function(env, "ops-softmax", e_ops_softmax);


  emacs_value e_ops_reduce_sum = env->make_function(env,
                                                 2,
                                                 2,
                                                 reduce_sum,
                                                 "Sum the elements along a given axis.",
                                                 NULL);
  bind_function(env, "ops-reduce-sum", e_ops_reduce_sum);


  emacs_value e_ops_reduce_mean = env->make_function(env,
                                                 2,
                                                 2,
                                                 reduce_mean,
                                                 "Calculate the mean of the elements along a given axis.",
                                                 NULL);
  bind_function(env, "ops-reduce-mean", e_ops_reduce_mean);
  
  emacs_value e_ops_relu = env->make_function(env,
                                                 1,
                                                 2,
                                                 relu,
                                                 "Apply ReLU activation",
                                                 NULL);
  bind_function(env, "ops-relu", e_ops_relu);

  emacs_value e_ops_heaviside = env->make_function(env,
                                              1,
                                              1,
                                              heaviside,
                                              "Apply Heaviside activation",
                                              NULL);
  bind_function(env, "ops-heaviside", e_ops_heaviside);  


  emacs_value e_ops_identity = env->make_function(env,
                                                 1,
                                                 2,
                                                 identity,
                                                 "Apply identity.",
                                                 NULL);
  bind_function(env, "ops-identity", e_ops_identity);

  emacs_value e_ops_natural_log = env->make_function(env,
                                                 1,
                                                 2,
                                                 natural_log,
                                                 "Apply natural logarithm.",
                                                 NULL);
  bind_function(env, "ops-log", e_ops_natural_log);



  provide(env, "ops");
  /* loaded successfully */
  return 0;
}
