#pragma once
#include "ast.hpp"
#include <string>
using namespace std;

namespace calc {

struct SeriesResult {
    bool ok;
    string text;
    string latex;
    string message;
};

SeriesResult taylorSeries(const ExprPtr& expr, const string& var,
                          double a, int terms);

}
