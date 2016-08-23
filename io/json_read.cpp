#include "private/json.hpp"
#include "mex.h"
#include <fstream>
#include <vector>

using json = nlohmann::json;

mxArray* recurse_json(const json &obj)
{
    mxArray* out = NULL;
    // Create a MATLAB structure from the the name,value pairs
    switch (obj.type())
    {
        case json::value_t::object:
        {
            std::vector<const char*> fieldname_ptrs(obj.size());
            std::vector<std::string> fieldnames(obj.size());
            int a = 0;
            for (json::const_iterator it = obj.begin(); it != obj.end(); ++it, ++a) {
                fieldnames[a] = it.key();
                fieldname_ptrs[a] = fieldnames[a].c_str();
            }
            out = mxCreateStructMatrix(1, 1, obj.size(), &fieldname_ptrs[0]);
            a = 0;
            for (json::const_iterator it = obj.begin(); it != obj.end(); ++it, ++a)
                mxSetFieldByNumber(out, 0, a, recurse_json(it.value()));
            break;
        }
        case json::value_t::array:
        {
            out = mxCreateCellMatrix(1, obj.size());
            int a = 0;
            for (json::const_iterator it = obj.begin(); it != obj.end(); ++it, ++a)
                mxSetCell(out, a, recurse_json(it.value()));
            break;
        }
        case json::value_t::string:
            out = mxCreateString(obj.get<std::string>().c_str());
            break;
        case json::value_t::boolean:
            out = mxCreateLogicalScalar(mxLogical(obj.get<bool>()));
            break;
        case json::value_t::number_integer:
            out = mxCreateDoubleScalar(double(obj.get<long long>()));
            break;
        case json::value_t::number_unsigned:
            out = mxCreateDoubleScalar(double(obj.get<unsigned long long>()));
            break;
        case json::value_t::number_float:
            out = mxCreateDoubleScalar(double(obj.get<double>()));
            break;
        default:
            out = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            break;
    }
    return out;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check number of arguments
	if (nlhs != 1 || nrhs != 1)
		mexErrMsgTxt("Unexpected number of arguments.");
    const char* fname = mxArrayToString(prhs[0]);
    if (fname == NULL )
        mexErrMsgTxt("Input argument expected to be a filename.");
    
    // Open the file
    std::ifstream fs(fname);
    if (!fs.is_open())
        mexErrMsgTxt("Failed to open file.");
    
    // Parse the file
    json j;
    try {
        fs >> j;
        fs.close();
    } catch (...) {
        fs.close();
        mexErrMsgTxt("Failed to parse file.");
    }
    
    // Create the MATLAB structure
    plhs[0] = recurse_json(j);
}