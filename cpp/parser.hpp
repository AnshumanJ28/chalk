#pragma once
#include "ast.hpp"
#include <string>
using namespace std;

namespace calc {

ExprPtr parse(const string& text);

}
