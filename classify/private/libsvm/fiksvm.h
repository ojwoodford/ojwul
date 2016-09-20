#ifndef _LIBFIKSVM_H
#define _LIBFIKSVM_H

#include "svm.h"
#include "mex.h"
#include "matrix.h"

//change this to single or double precision
typedef float Double;
static mxClassID mxDOUBLE_CLASS_ID = mxSINGLE_CLASS;

// typedef double Double;
// static mxClassID mxDOUBLE_CLASS_ID = mxDOUBLE_CLASS;

/* keep the support vectors in the sorted order */
class fiksvm_exact_classifier{
 public:
  fiksvm_exact_classifier(const struct svm_model* model);
  fiksvm_exact_classifier(const mxArray *matlab_struct);
  const char* model_to_matlab(mxArray *plhs[]);
  Double predict(const Double* x);
  ~fiksvm_exact_classifier();

 private:
  void precompute_per_dim(const struct svm_model *model,int dim,const Double * dense_sv);
  int feat_dim; 
  int num_sv; 
  Double rho;
  Double** sort_sv;
  Double** cum_alpha;
  Double** cum_alphax;

};

/* keep an approximation of the model */
class fiksvm_approx_classifier{
public:
  fiksvm_approx_classifier(const struct svm_model* model, int num_bins);
  fiksvm_approx_classifier(const mxArray *matlab_struct);
  const char* model_to_matlab(mxArray *plhs[]);
  Double pwl_predict(const Double* x);
  Double pwc_predict(const Double* x);
  ~fiksvm_approx_classifier();
  
 private:
  void precompute_per_dim(const struct svm_model *model,int dim,const Double* dense_sv);
  void compute_bin_indx(int &indx, Double dim_value, int dim_indx);
  void compute_lrbin_indx(int &l_indx, int &r_indx, Double &alpha, Double dim_value, int dim_indx);
  int feat_dim; 
  int num_sv; 
  int num_bins; 

  Double  rho;
  Double* min_sv;
  Double* max_sv;

  //rank interpolation
  Double*  a;
  Double*  b;

  //value of the decision function
  Double** h;
};
#endif
