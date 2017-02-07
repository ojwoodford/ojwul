#include "private/json.hpp"
#include "mex.h"
#include <fstream>
#include <vector>

using json = nlohmann::json;

json recurse_object(const mxArray* in)
{
    json obj;      
    // Go over all types
    if (mxIsCell(in)) {
        int N = mxGetNumberOfElements(in);
        for (int a = 0; a < N; ++a)
            obj.push_back(recurse_object(mxGetCell(in, a)));
    } else if (mxIsStruct(in)) {
        int M = mxGetNumberOfFields(in);
        std::vector<const char*> fnames(M);
        for (int b = 0; b < M; ++b)
            fnames[b] = mxGetFieldNameByNumber(in, b);
        int N = mxGetNumberOfElements(in);
        for (int a = 0; a < N; ++a) {
            json obj_;
            for (int b = 0; b < M; ++b)
                obj_[fnames[b]] = recurse_object(mxGetFieldByNumber(in, a, b));
            if (N == 1)
                obj = obj_;
            else
                obj.push_back(obj_);
        }
    } else if (mxIsChar(in)) {
        obj = mxArrayToString(in);
    } else if (mxIsLogical(in)) {
        int N = mxGetNumberOfElements(in);
        const mxLogical* data = (const mxLogical*)mxGetData(in);
        if (N == 1) {
            obj = data[0] != 0;
        } else {
            for (int a = 0; a < N; ++a)
                obj.push_back(data[a] != 0);
        }
    } else if (mxIsNumeric(in)) {
        int N = mxGetNumberOfElements(in);
        const void* data = (const void*)mxGetData(in);
        switch (mxGetClassID(in)) {
            case mxDOUBLE_CLASS:
                if (N == 1) { obj = ((const double*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const double*)data)[a]); }
                break;
            case mxSINGLE_CLASS:
                if (N == 1) { obj = ((const float*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const float*)data)[a]); }
                break;
            case mxUINT8_CLASS:
                if (N == 1) { obj = ((const uint8_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const uint8_t*)data)[a]); }
                break;
            case mxINT8_CLASS:
                if (N == 1) { obj = ((const int8_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const int8_t*)data)[a]); }
                break;
            case mxUINT16_CLASS:
                if (N == 1) { obj = ((const uint16_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const uint16_t*)data)[a]); }
                break;
            case mxINT16_CLASS:
                if (N == 1) { obj = ((const int16_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const int16_t*)data)[a]); }
                break;
            case mxUINT32_CLASS:
                if (N == 1) { obj = ((const uint32_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const uint32_t*)data)[a]); }
                break;
            case mxINT32_CLASS:
                if (N == 1) { obj = ((const int32_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const int32_t*)data)[a]); }
                break;
            case mxUINT64_CLASS:
                if (N == 1) { obj = ((const uint64_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const uint64_t*)data)[a]); }
                break;
            case mxINT64_CLASS:
                if (N == 1) { obj = ((const int64_t*)data)[0]; } else {
                for (int a = 0; a < N; ++a) obj.push_back(((const int64_t*)data)[a]); }
                break;
            default:
                mexErrMsgTxt("Unsupported type.");
                break;
        }
    } else {
		mexErrMsgTxt("Unrecognized type.");
    }
    return obj;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check number of arguments
	if (nlhs != 0 || nrhs != 2)
		mexErrMsgTxt("Unexpected number of arguments.");
    const char* fname = mxArrayToString(prhs[1]);
    if (fname == NULL )
        mexErrMsgTxt("Second input argument expected to be a filename.");
        
    // Parse the object
    const json j = recurse_object(prhs[0]);
    
    // Write out the file
    std::ofstream fs(fname);
    if (!fs.is_open())
        mexErrMsgTxt("Failed to open file for writing.");
    std::stringstream jsonStream;
    jsonStream << std::setw(2) << std::setprecision(9) << j << std::endl;
    fs << jsonStream.str();
    fs.close();
}