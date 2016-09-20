#include "fiksvm.h"
#include "svm.h"
#include "mex.h"
#include <algorithm>
#include <list>
#include <limits>

typedef float M;
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif 

static const Double MIN_SV_PRECISION=1e-10;

// matlab struct fields for exact classifier
static const int FIKSVM_EXACT_NRF = 6;
static const char *fiksvm_exact_field_names[] = {
  "feat_dim",
  "num_sv",
  "rho",
  "sort_sv",
  "cum_alpha",
  "cum_alphax"
};

//matlab struct for the approximate classifier
static const int FIKSVM_APPROX_NRF = 9;
static const char * fiksvm_approx_field_names[] = {
  "feat_dim",
  "num_sv",
  "num_bins",
  "rho",
  "min_sv",
  "max_sv",
  "a",
  "b",
  "h"
};

using namespace std;
//wrapper for value and indx
class val_indx{
public:
  Double val; int indx;
  val_indx(Double val, int indx){
    this->val = val; this->indx = indx;
  }
  bool operator<(const val_indx &other) const{
    return this->val < other.val;
  }
};

//returns the index which sorts vals in increasing order
void sort_index(int *index, const Double *vals, int dim){
  val_indx *data = (val_indx*)malloc(sizeof(val_indx)*dim);
  for(int i = 0; i < dim;i++)
    data[i] = val_indx(vals[i],i);
  sort(data,data+dim);
  for(int i = 0; i < dim; i++)
    index[i] = data[i].indx;
  free(data);
};

//allocates 2d array of size m x n
Double  ** alloc_2d_array(int m, int n){
  Double ** buffer    = (Double**)malloc(m*sizeof(Double*));
  for(int i = 0; i < m ; i++)
    buffer[i] = (Double*)malloc(n*sizeof(Double));
  //zero out the array
  for(int i = 0; i < m ; i++)
    for(int j = 0; j < n; j++)
      buffer[i][j]=0;
  return buffer;
}
void free_2d_array(Double ** buffer, int m, int n){
  for(int i = 0; i < m ;i++)
    free(buffer[i]);
  free(buffer);
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
fiksvm_exact_classifier::fiksvm_exact_classifier(const mxArray *matlab_struct){
}
fiksvm_exact_classifier::fiksvm_exact_classifier(const struct svm_model *model){
  num_sv=model->l; 
  rho = model->rho[0];
  //compute the feat_dim
  feat_dim = -1;
  for(int i = 0; i < num_sv;i++){
    svm_node *px = model->SV[i];
    while(px->index !=-1){
      feat_dim = max(px->index,feat_dim);
      ++px;
    }
  }
  Double ** dense_sv   = alloc_2d_array(feat_dim,num_sv);
  sort_sv    = alloc_2d_array(feat_dim,num_sv);
  cum_alphax = alloc_2d_array(feat_dim,num_sv+1);
  cum_alpha  = alloc_2d_array(feat_dim,num_sv+1);

  //copy the support vectors
  for(int i = 0; i < num_sv; i++){
    svm_node *px = model->SV[i];
    while(px->index !=-1){
      dense_sv[px->index-1][i]=px->value;
      ++px;
    }
  }
  //loop over dimensions and compute various features
  for(int i = 0 ; i < feat_dim; i++)
    precompute_per_dim(model,i,dense_sv[i]);
  //free up memory
  free_2d_array(dense_sv,feat_dim,num_sv);
};
void fiksvm_exact_classifier::precompute_per_dim(const struct svm_model *model, 
						 int dim, const Double * dense_sv){
  int * indx = (int*)malloc(sizeof(int)*num_sv);
  sort_index(indx, dense_sv, num_sv);
  cum_alpha[dim][0]  = 0;
  cum_alphax[dim][0] = 0;
  for(int i = 0; i < num_sv; i++){
    int ik = indx[i];
    Double alpha = model->sv_coef[0][ik];
    Double x     = dense_sv[ik];
    sort_sv[dim][i]= x;
    cum_alpha[dim][i+1]  = cum_alpha[dim][i] +  alpha;
    cum_alphax[dim][i+1] = cum_alphax[dim][i] + alpha*x;
  }
  free(indx);
}


const char* fiksvm_exact_classifier::model_to_matlab(mxArray *plhs[]){
  Double *ptr;
  mxArray *return_model, **rhs;
  int out_id = 0;
  rhs = (mxArray **)mxMalloc(sizeof(mxArray *)*FIKSVM_EXACT_NRF);
  
  //feat_dim
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = feat_dim;
  out_id++;

  //num_sv
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = num_sv;
  out_id++;

  //rho
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = rho;
  out_id++;

  //sort_sv
  int rhs_indx;
  rhs[out_id] = mxCreateNumericMatrix(feat_dim,num_sv,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  rhs_indx = 0;
  for(int j = 0; j < num_sv ; j++)
    for(int i = 0; i < feat_dim; i++)
      ptr[rhs_indx++] =  sort_sv[i][j];
  out_id++;
      
  //cum_alpha
  rhs[out_id] = mxCreateNumericMatrix(feat_dim,num_sv+1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  rhs_indx = 0;
  
  //matlab arrays are row first indices
  for(int j = 0; j < num_sv+1 ; j++)
    for(int i = 0; i < feat_dim; i++)
      ptr[rhs_indx++] =  cum_alpha[i][j];

  out_id++;

  //cum_alphax
  rhs[out_id] = mxCreateNumericMatrix(feat_dim,num_sv+1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  rhs_indx = 0;
  for(int j = 0; j < num_sv+1 ; j++)
    for(int i = 0; i < feat_dim; i++)
      ptr[rhs_indx++] =  cum_alphax[i][j];
  out_id++;

  
  /* Create a struct matrix contains NUM_OF_RETURN_FIELD fields */
  return_model = mxCreateStructMatrix(1, 1, FIKSVM_EXACT_NRF, fiksvm_exact_field_names);
  
  /* Fill struct matrix with input arguments */
  for(int i = 0; i < FIKSVM_EXACT_NRF; i++)
    mxSetField(return_model,0,fiksvm_exact_field_names[i],mxDuplicateArray(rhs[i]));
  /* return */
  plhs[0] = return_model;
  mxFree(rhs);
  return NULL;
}
Double fiksvm_exact_classifier::predict(const Double* x){
  return -1;
}

// destructor
fiksvm_exact_classifier::~fiksvm_exact_classifier(){
  free_2d_array(sort_sv,feat_dim,num_sv);
  free_2d_array(cum_alpha,feat_dim,num_sv+1);
  free_2d_array(cum_alphax,feat_dim,num_sv+1);
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
fiksvm_approx_classifier::fiksvm_approx_classifier(const mxArray *matlab_struct){
  int i,j,num_of_fields;
  Double *ptr;
  int id = 0;

  mxArray **rhs;
  num_of_fields = mxGetNumberOfFields(matlab_struct);
  rhs = (mxArray **) mxMalloc(sizeof(mxArray *)*num_of_fields);
  
  for(i=0;i<num_of_fields;i++)
    rhs[i] = mxGetFieldByNumber(matlab_struct, 0, i);

  //feature dimension
  ptr = (Double *)mxGetPr(rhs[id++]);
  feat_dim = (int)ptr[0];

  //number of sup vectors
  ptr = (Double *)mxGetPr(rhs[id++]);
  num_sv = (int)ptr[0];

  //number of bins to approximate
  ptr = (Double *)mxGetPr(rhs[id++]);
  num_bins = (int)ptr[0];

  //model->rho
  ptr = (Double *)mxGetPr(rhs[id++]);
  rho = ptr[0];

  // allocate memory for each dimension
  h = alloc_2d_array(feat_dim, num_bins+1);
  a = (Double *)malloc(feat_dim*sizeof(Double));
  b = (Double *)malloc(feat_dim*sizeof(Double));
  min_sv = (Double *)malloc(feat_dim*sizeof(Double));
  max_sv = (Double *)malloc(feat_dim*sizeof(Double));
  
  // min_sv
  ptr = (Double *)mxGetPr(rhs[id++]);
  for(i=0; i<feat_dim;i++)
    min_sv[i] = ptr[i];

  // max_sv
  ptr = (Double *)mxGetPr(rhs[id++]);
  for(i=0; i<feat_dim;i++)
    max_sv[i] = ptr[i];

  // a
  ptr = (Double *)mxGetPr(rhs[id++]);
  for(i=0; i<feat_dim;i++)
    a[i] = ptr[i];

  // b
  ptr = (Double *)mxGetPr(rhs[id++]);
  for(i=0; i<feat_dim;i++)
    b[i] = ptr[i];

  // h
  ptr = (Double *)mxGetPr(rhs[id++]);
  int rhs_indx = 0;
  for(j = 0; j < num_bins+1 ; j++)
    for(i = 0; i < feat_dim; i++)
      h[i][j] = ptr[rhs_indx++];
}
fiksvm_approx_classifier::fiksvm_approx_classifier(const struct svm_model * model, int num_bins){
  num_sv=model->l; 
  //compute the feat_dim
  feat_dim = -1;
  for(int i = 0; i < num_sv;i++){
    svm_node *px = model->SV[i];
    while(px->index !=-1){
      feat_dim = max(px->index,feat_dim);
      ++px;
    }
  }
  this->num_bins = num_bins;
  rho = model->rho[0];
  
  // allocate memory for each dimension
  h = alloc_2d_array(feat_dim, num_bins+1);
  a = (Double *)malloc(feat_dim*sizeof(Double));
  b = (Double *)malloc(feat_dim*sizeof(Double));
  min_sv = (Double *)malloc(feat_dim*sizeof(Double));
  max_sv = (Double *)malloc(feat_dim*sizeof(Double));

  Double ** dense_sv   = alloc_2d_array(feat_dim,num_sv);
  //copy the support vectors
  for(int i = 0; i < num_sv; i++){
    svm_node *px = model->SV[i];
    while(px->index !=-1){
      dense_sv[px->index-1][i]=px->value;
      ++px;
    }
  }
  for(int i = 0; i < feat_dim; i++){
    precompute_per_dim(model, i, dense_sv[i]);
  }
  free_2d_array(dense_sv,feat_dim,num_sv);
}

// precompute the per dimension features
void fiksvm_approx_classifier::precompute_per_dim(const struct svm_model* model, int dim, const Double * dense_sv){
  Double min_val= numeric_limits<Double>::max();
  Double max_val=-numeric_limits<Double>::max();
  for(int i = 0; i < num_sv; i++){
    if(dense_sv[i] >= max_val)
      max_val = dense_sv[i];
    if(dense_sv[i] <= min_val)
      min_val = dense_sv[i];
  }
  min_sv[dim] = min_val;
  max_sv[dim] = max_val;
  if( (max_val - min_val) < MIN_SV_PRECISION) {
    // make them all zeros
    a[dim] = 0; b[dim] = 0;
  }else{
    //set the coefs for linear interpolation
    Double step_size = (max_val - min_val)/num_bins;
    Double x_val, h_val;
    a[dim] = num_bins/(max_val - min_val);
    b[dim] = -min_val*num_bins/(max_val - min_val);
    for(int i = 0; i < num_bins + 1; i++){
      x_val = min_val + i*step_size;
      h_val = 0;
      for(int j = 0; j < num_sv ; j++){
	Double alpha = model->sv_coef[0][j];
	h_val += alpha*min(x_val, dense_sv[j]);
      }
      h[dim][i] = h_val;
    }
  }
}

const char* fiksvm_approx_classifier::model_to_matlab(mxArray *plhs[]){
  Double *ptr;
  mxArray *return_model, **rhs;
  int out_id = 0;
  rhs = (mxArray **)mxMalloc(sizeof(mxArray *)*FIKSVM_APPROX_NRF);
  
  //feat_dim
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = feat_dim;
  out_id++;

  //num_sv
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = num_sv;
  out_id++;

  //num_bins;
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = num_bins;
  out_id++;

  //rho
  rhs[out_id] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  ptr[0] = rho;
  out_id++;

  //min_sv
  rhs[out_id] = mxCreateNumericMatrix(1,feat_dim,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  for(int i = 0; i < feat_dim; i++)
    ptr[i] =  min_sv[i];
  out_id++;
      
  //max_sv
  rhs[out_id] = mxCreateNumericMatrix(1,feat_dim,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  for(int i = 0; i < feat_dim; i++)
    ptr[i] =  max_sv[i];
  out_id++;

  //a
  rhs[out_id] = mxCreateNumericMatrix(1,feat_dim,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  for(int i = 0; i < feat_dim; i++)
    ptr[i] =  a[i];
  out_id++;

  //b
  rhs[out_id] = mxCreateNumericMatrix(1,feat_dim,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  for(int i = 0; i < feat_dim; i++)
    ptr[i] =  b[i];
  out_id++;

  //h
  rhs[out_id] = mxCreateNumericMatrix(feat_dim,num_bins+1,mxDOUBLE_CLASS_ID,mxREAL);
  ptr = (Double *)mxGetPr(rhs[out_id]);
  int rhs_indx = 0;
  //matlab arrays are row first indices
  for(int j = 0; j < num_bins+1 ; j++)
    for(int i = 0; i < feat_dim; i++)
      ptr[rhs_indx++] =  h[i][j];
  out_id++;
  
  /* Create a struct matrix contains fields */
  return_model = mxCreateStructMatrix(1, 1, FIKSVM_APPROX_NRF, fiksvm_approx_field_names);
  
  /* Fill struct matrix with input arguments */
  for(int i = 0; i < FIKSVM_APPROX_NRF; i++)
    mxSetField(return_model,0,fiksvm_approx_field_names[i],mxDuplicateArray(rhs[i]));
  /* return */
  plhs[0] = return_model;
  mxFree(rhs);
  return NULL;
}
void fiksvm_approx_classifier::compute_bin_indx(int &indx, Double dim_value, int dim_indx){
  Double f_indx;
  f_indx = a[dim_indx]*dim_value + b[dim_indx];
  indx   = (int)(f_indx+0.5); //round to nearest integer
  indx   = min(max(0,indx),num_bins);
}
void fiksvm_approx_classifier::compute_lrbin_indx(int &l_indx, int &r_indx, Double &alpha, Double dim_value, int dim_indx){
  Double f_indx;
  f_indx = a[dim_indx]*dim_value + b[dim_indx];
  l_indx = (int)f_indx;
  r_indx = l_indx+1;
  alpha  = f_indx-l_indx;
  l_indx = min(max(0,l_indx),num_bins);
  r_indx = min(max(0,r_indx),num_bins);
}
Double fiksvm_approx_classifier::pwl_predict(const Double* x){
  Double dec_value = -rho, alpha;
  int l_indx, r_indx;
  for(int i = 0; i < feat_dim; i++){
    compute_lrbin_indx(l_indx,r_indx,alpha,x[i],i);
    dec_value += h[i][l_indx]*(1-alpha) + h[i][r_indx]*alpha;
  }
  return dec_value;
}
Double fiksvm_approx_classifier::pwc_predict(const Double *x){
  Double dec_value = -rho;
  int indx;
  for(int i=0;i<feat_dim;i++){
    compute_bin_indx(indx,x[i],i);
    dec_value += h[i][indx];
  }
  return dec_value;
}

fiksvm_approx_classifier::~fiksvm_approx_classifier(){
  free(a);
  free(b);
  free_2d_array(h,feat_dim,num_bins+1);
  free(min_sv);
  free(max_sv);
}
