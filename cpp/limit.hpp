#pragma once
#include "ast.hpp"
#include <string>
using namespace std;

namespace calc {

struct LimitResult {
    bool ok;
    string text;
    string latex;
};

LimitResult limit(const ExprPtr& e, const string& v,
                  double point, bool approachInfPos, bool approachInfNeg);

}
