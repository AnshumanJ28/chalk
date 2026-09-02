#pragma once
#include "ast.hpp"
#include <string>
#include <complex>

namespace calc {

using namespace std;

complex<double> evaluate_complex(ExprPtr expr, const string& var, complex<double> val);
complex<double> compute_residue(const string& expr_str, const string& var, complex<double> z0);
complex<double> contour_integral_circle(const string& expr_str, const string& var, complex<double> center, double radius, int steps = 1000);

}
