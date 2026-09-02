#pragma once
#include "ast.hpp"
#include <string>
using namespace std;

namespace calc {

string toString(const ExprPtr& e);

string toLatex(const ExprPtr& e);

}
