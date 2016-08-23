#include "private/json.hpp"
#include "mex.h"
#include <fstream>
#include <vector>

using json = nlohmann::json;

mxArray* recurse_json(const json &obj);

mxArray* hetergeneous_array(const json &obj)
{
    mxArray* out = mxCreateCellMatrix(1, obj.size());
    int a = 0;
    for (json::const_iterator it = obj.begin(); it != obj.end(); ++it, ++a)
        mxSetCell(out, a, recurse_json(it.value()));
    return out;
}

mxArray* recurse_json(const json &obj)
{
    mxArray* out;
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
            out = mxCreateStructMatrix(1, 1, fieldname_ptrs.size(), &fieldname_ptrs[0]);
            a = 0;
            for (json::const_iterator it = obj.begin(); it != obj.end(); ++it, ++a)
                mxSetFieldByNumber(out, 0, a, recurse_json(it.value()));
            break;
        }
        case json::value_t::array:
        {
            // Check for special cases
            if (obj.size() == 0) {
                out = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
                break;
            }
            json::const_iterator it = obj.begin();
            if (obj.size() == 1) {
                out = recurse_json(it.value());
                break;
            }
            json::value_t type = it->type();
            if (type == json::value_t::array || type == json::value_t::string) {
                // Array of arrays or strings - assume heterogeneous
                out = hetergeneous_array(obj);
                break;
            }
            // Check if homogeneous or heterogeneous
            bool homogeneous = true;
            for (++it; it != obj.end() && homogeneous; ++it)
                homogeneous = type == it->type();
            if (!homogeneous) {
                out = hetergeneous_array(obj);
                break;
            }
            if (type == json::value_t::object) {
                // Check the fieldnames are the same
                const json& obj_ = obj.begin().value();
                std::vector<const char*> fieldname_ptrs(obj_.size());
                std::vector<std::string> fieldnames(obj_.size());
                int a = 0;
                for (it = obj_.begin(); it != obj_.end(); ++it, ++a) {
                    fieldnames[a] = it.key();
                    fieldname_ptrs[a] = fieldnames[a].c_str();
                }
                for (it = obj.begin()+1; (it != obj.end()) && homogeneous; ++it) {
                    if (it->size() != fieldnames.size()) {
                        homogeneous = false;
                        break;
                    }
                    for (std::vector<const char*>::const_iterator str = fieldname_ptrs.begin(); str != fieldname_ptrs.end(); ++str) {
                        if (it->count(*str) == 0) {
                            homogeneous = false;
                            break;
                        }
                    }
                }
                if (!homogeneous) {
                    out = hetergeneous_array(obj);
                    break;
                }
                // Homogeneous, so make an array of structs
                out = mxCreateStructMatrix(1, obj.size(), fieldname_ptrs.size(), &fieldname_ptrs[0]);
                a = 0;
                for (it = obj.begin(); it != obj.end(); ++it, ++a) {
                    int b = 0;
                    for (std::vector<const char*>::const_iterator str = fieldname_ptrs.begin(); str != fieldname_ptrs.end(); ++str, ++b)
                        mxSetFieldByNumber(out, a, b, recurse_json(it->find(*str).value()));
                }
            }
            switch (type) {
                case json::value_t::boolean:
                {
                    out = mxCreateLogicalMatrix(1, obj.size());
                    mxLogical* data = (mxLogical*)mxGetData(out);
                    for (it = obj.begin(); it != obj.end(); ++it)
                        *data++ = mxLogical(it->get<bool>());
                    break;
                }
                case json::value_t::number_unsigned:
                case json::value_t::number_integer:
                {
                    out = mxCreateNumericMatrix(1, obj.size(), mxINT64_CLASS, mxREAL);
                    int64_t* data = (int64_t*)mxGetData(out);
                    for (it = obj.begin(); it != obj.end(); ++it)
                        *data++ = it->get<int64_t>();
                    break;
                }
                case json::value_t::number_float:
                {
                    out = mxCreateNumericMatrix(1, obj.size(), mxDOUBLE_CLASS, mxREAL);
                    double* data = (double*)mxGetData(out);
                    for (it = obj.begin(); it != obj.end(); ++it)
                        *data++ = it->get<double>();
                    break;
                }
                default:
                    break;
            }
            break;
        }
        case json::value_t::string:
            out = mxCreateString(obj.get<std::string>().c_str());
            break;
        case json::value_t::boolean:
            out = mxCreateLogicalScalar(mxLogical(obj.get<bool>()));
            break;
        case json::value_t::number_integer:
            out = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
            *((int64_t*)mxGetData(out)) = obj.get<int64_t>();
            break;
        case json::value_t::number_unsigned:
            out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
            *((uint64_t*)mxGetData(out)) = obj.get<uint64_t>();
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