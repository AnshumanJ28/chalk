#pragma once
#include "ast.hpp"
#include <string>
#include <vector>
using namespace std;

namespace calc {

struct VectorCalcResult {
    bool ok;
    vector<string> components;
    string scalar;
    vector<string> latex_components;
    string latex_scalar;
    string message;
};

VectorCalcResult gradient(const string& expr, const vector<string>& vars);

VectorCalcResult divergence(const vector<string>& componentExprs,
                            const vector<string>& vars);

VectorCalcResult curl(const vector<string>& componentExprs,
                      const vector<string>& vars);

VectorCalcResult laplacian(const string& expr, const vector<string>& vars);

double line_integral(const vector<string>& F_exprs, const vector<string>& path_exprs, const string& t_var, double t_start, double t_end, int steps=1000);
double surface_integral(const vector<string>& F_exprs, const vector<string>& surface_exprs, const vector<string>& uv_vars, double u_start, double u_end, double v_start, double v_end, int steps=100);

}
