#pragma once
#include "ast.hpp"
#include <string>

namespace calc {

using namespace std;

// Compute the continuous Fourier transform F(w) = \int f(t) e^{-iwt} dt
// Returns a string representation (since it involves delta functions which aren't in our standard AST)
string fourier_transform(const string& expr_str, const string& t_var = "t", const string& w_var = "w");

}
