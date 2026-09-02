#pragma once
#include "ast.hpp"
#include <string>
using namespace std;

namespace calc {

struct LaplaceResult {
    bool ok;
    string text;
    string latex;
    string method;
    string message;
};

LaplaceResult laplaceTransform(const ExprPtr& expr, const string& t_var,
                               const string& s_var);

LaplaceResult inverseLaplaceTransform(const ExprPtr& expr, const string& s_var,
                                      const string& t_var);

}
