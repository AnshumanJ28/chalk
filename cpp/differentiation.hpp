#pragma once
#include "ast.hpp"
#include <string>
#include <vector>
using namespace std;

namespace calc {

ExprPtr differentiate(const ExprPtr& e, const string& v);

ExprPtr differentiateOrder(const ExprPtr& e, const string& v, int order);

ExprPtr mixedPartial(const ExprPtr& e, const vector<string>& vars);

}
