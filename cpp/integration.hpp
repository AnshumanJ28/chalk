#pragma once
#include "ast.hpp"
#include <string>
#include <vector>
using namespace std;

namespace calc {

struct CalcResult {
    bool ok;
    string text;
    string latex;
};

CalcResult integrate(const ExprPtr& e, const string& wrt);

struct DefiniteResult {
    bool ok;
    double value;
    string text;
    string latex;
    bool is_numeric;
};

DefiniteResult definiteIntegral(const ExprPtr& e, const string& wrt,
                                double lower, double upper);

struct IntegralSpec {
    string var;
    double lower;
    double upper;
};

struct MultiIntegralResult {
    bool ok;
    string text;
    string latex;
    bool is_numeric;
};

MultiIntegralResult multiIntegral(const ExprPtr& e, const vector<IntegralSpec>& specs);

}
