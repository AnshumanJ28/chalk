#pragma once
#include "ast.hpp"
#include <string>
#include <vector>
#include <map>
using namespace std;

namespace calc {

struct SolveResult {
    bool ok;
    vector<string> roots;
    vector<string> latex_roots;
    string message;
};

SolveResult solve(const ExprPtr& e, const string& wrt);

struct SystemSolveResult {
    bool ok;
    map<string, double> solution;
    string message;
};

SystemSolveResult solveSystem(const vector<string>& equations,
                              const vector<string>& variables);

}
