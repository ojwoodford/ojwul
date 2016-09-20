#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "svm.h"
#include "fiksvm.h"
#include "mex.h"
#include "svm_model_matlab.h"

#if _MATLAB_ < 703
typedef int mwIndex;
#endif 

#define CMD_LEN 2048

using namespace std;

void exit_with_help()
{
	mexPrintf(
	"Usage: model = precomp_model(model,'-m 0/1 -n num_bins')\n"
	"\t -m 0/1      : 0 - exact model [default]\n"
	"\t             : 1 - approx model\n"
	"\t -n num_bins : number of bins for approximation\n"
	"\t               default=100; num_bins>=1\n"
	);
}
static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] ){
  struct svm_model *model;
  /*default values for the model computation */
  int model_type = 0;
  int num_bins = 100;

  if(nlhs > 1 || nrhs < 1){
    exit_with_help();
    fake_answer(plhs);
    return;
  }
  if(mxIsStruct(prhs[0])){
    const char *error_msg;
    model = (struct svm_model *) malloc(sizeof(struct svm_model));
    error_msg = matlab_matrix_to_model(model, prhs[0]);
    if(error_msg){
      mexPrintf("Error: can't read model: %s\n", error_msg);
      fake_answer(plhs);
    }else{ /*parse the command line arguments */
      if(nrhs==2){
	int i, argc = 1;
	char cmd[CMD_LEN], *argv[CMD_LEN/2];
	
	/* put options in argv[] */
	mxGetString(prhs[1], cmd,  mxGetN(prhs[1]) + 1);
	if((argv[argc] = strtok(cmd, " ")) != NULL)
	  while((argv[++argc] = strtok(NULL, " ")) != NULL)
	    ;
	
	for(i=1;i<argc;i++){
	  if(argv[i][0] != '-') break;
	  if(++i>=argc){
	    exit_with_help();
	    fake_answer(plhs);
	    return;
	  }
	  switch(argv[i-1][1]){
	  case 'm':
	    model_type = atoi(argv[i]);
	    break;
	  case 'n':
	    num_bins = atoi(argv[i]);
	    break;
	  default:
	    mexPrintf("unknown option\n");
	    svm_destroy_model(model);
	    exit_with_help();
	    fake_answer(plhs);
	    return;
	  }
	}
      }
      /*check validity of input parameters */
      if(num_bins < 1 || (model_type != 0 && model_type !=1)){
	svm_destroy_model(model);
	exit_with_help();
	fake_answer(plhs);
	return;
      }
      /*check if the model type is compatible */
      if(model->nr_class != 2){
	mexPrintf("Error : only binary classifier supported..\n");
	svm_destroy_model(model);
	fake_answer(plhs);
	return;
      }
      int svm_type=svm_get_svm_type(model);
      if(svm_type != C_SVC){
	mexPrintf("Error : only svc classifier supported..\n");
	svm_destroy_model(model);
	fake_answer(plhs);
	return;
      }
      /*write the output to the rhs */
      const char *error_msg_fiksvm;
      if(model_type == 0){
	fiksvm_exact_classifier *exact_model = new fiksvm_exact_classifier(model);
	error_msg_fiksvm = exact_model->model_to_matlab(plhs); 
	delete exact_model;
      }else{
	fiksvm_approx_classifier *approx_model = new fiksvm_approx_classifier(model,num_bins);
	error_msg_fiksvm = approx_model->model_to_matlab(plhs); 
	delete approx_model;
      }
      if(error_msg_fiksvm){
	mexPrintf("Error: can't precompute model: %s\n", error_msg_fiksvm);
	fake_answer(plhs);
      }
    }
    svm_destroy_model(model);
    return;
  }else{ /*arg 0 is not model file */
    mexPrintf("model file should be a struct array\n");
    fake_answer(plhs);
  }
  return;
}

