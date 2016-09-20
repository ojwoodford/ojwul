#include <time.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "svm.h"
#include "fiksvm.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif 

#define CMD_LEN 2048

void exit_with_help()
{
	mexPrintf(
	" Usage:" 
        " dec_values = fiksvm_predict(testing_label_vector, testing_instance_matrix, model,'options')\n"
	" \n"
	" Output:\n"
	"           dec_values      : predictions using model\n"
	" \n"
	" Options:\n"
	"   -e exact flag           : whether input model is exact or approx, 0 or 1 (default 0);\n"
	"   -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0);\n"
	"   -v verbose flag         : 0 or 1 (default 0);\n"
	"   -a approximation type   : 0 - piecewise const or 1 - piecewise linear (default 1);\n"
	"\n"
	);
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void predict(int nlhs, mxArray *plhs[], const mxArray *prhs[], 
	     const int predict_probability, const int verbose_flag, const int exact_flag, const int approx_type){
  //read the data and get the parameters
  int label_vector_row_num, label_vector_col_num;
  int feat_dim, num_feats;
  
  Double *ptr_instance; 
  Double *ptr_dec_values;

  fiksvm_exact_classifier *exact_classifier = NULL;
  fiksvm_approx_classifier *approx_classifier = NULL;

  feat_dim  = mxGetN(prhs[1]);
  num_feats = mxGetM(prhs[1]);
  label_vector_row_num = mxGetM(prhs[0]);
  label_vector_col_num = mxGetN(prhs[0]);

  if(label_vector_row_num!=num_feats){
      mexPrintf("# of labels (# of column in 1st argument) does not match # of instances (# of rows in 2nd argument).\n");
      fake_answer(plhs);
      return;
  }
  if(label_vector_col_num!=1){
      mexPrintf("label (1st argument) should be a vector (# of column is 1).\n");
      fake_answer(plhs);
      return;
  }
  plhs[0] = mxCreateNumericMatrix(num_feats, 1, mxDOUBLE_CLASS_ID, mxREAL);
  ptr_dec_values = (Double *)mxGetPr(plhs[0]);
  ptr_instance = (Double *)mxGetPr(prhs[1]);

  if(exact_flag)
    exact_classifier = new fiksvm_exact_classifier(prhs[2]);
  else
    approx_classifier = new fiksvm_approx_classifier(prhs[2]);

  

  if(mxIsSparse(prhs[1])){
    mexPrintf("Error : sparse data format not supported..\n");
    fake_answer(plhs);
    return;
  }
  Double *feat_data = (Double*)malloc(feat_dim*sizeof(Double));
  for(int i=0; i<num_feats; i++){
    for(int j=0; j<feat_dim ; j++){
      //matlab arrays go from down first and then sideways
      feat_data[j] = (Double)ptr_instance[num_feats*j + i];
    }
    if(exact_flag)
      ptr_dec_values[i] = exact_classifier->predict(feat_data);
    else{
      if(approx_type == 0)
	ptr_dec_values[i] = approx_classifier->pwc_predict(feat_data);
      else
	ptr_dec_values[i] = approx_classifier->pwl_predict(feat_data);
    }
  }
  if(exact_flag)
    delete exact_classifier;
  else
    delete approx_classifier;

  free(feat_data);
}
//outputs the predictions
void mexFunction( int nlhs, mxArray *plhs[],
  		 int nrhs, const mxArray *prhs[] )
{
	int prob_estimate_flag = 0;
	int verbose_flag = 0; 
	int approx_type  = 1;
	int exact_flag   = 0;
	if(nrhs > 4 || nrhs < 3)
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}
	if(mxIsStruct(prhs[2]))
	{
		if(nrhs==4)
		{
			int i, argc = 1;
			char cmd[CMD_LEN], *argv[CMD_LEN/2];

			// put options in argv[]
			mxGetString(prhs[3], cmd,  mxGetN(prhs[3]) + 1);
			if((argv[argc] = strtok(cmd, " ")) != NULL)
				while((argv[++argc] = strtok(NULL, " ")) != NULL)
					;

			for(i=1;i<argc;i++)
			{
				if(argv[i][0] != '-') break;
				if(++i>=argc)
				{
					exit_with_help();
					fake_answer(plhs);
					return;
				}
				switch(argv[i-1][1])
				{
					case 'b':
						prob_estimate_flag = atoi(argv[i]);
						break;
					case 'v':
						verbose_flag = atoi(argv[i]);
						break;
					case 'a':
						approx_type = atoi(argv[i]);
						break;
					case 'e':
						exact_flag = atoi(argv[i]);
						break;
					default:
						mexPrintf("unknown option\n");
						exit_with_help();
						fake_answer(plhs);
						return;
				}
			}
		}
		predict(nlhs, plhs, prhs,prob_estimate_flag,verbose_flag,exact_flag,approx_type);
	}
	else
	{
		mexPrintf("model file should be a struct array\n");
		fake_answer(plhs);
	}

	return;
}
