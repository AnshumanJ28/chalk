#pragma once
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include <stdexcept>
using namespace std;

namespace calc {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

enum class Op {
    NUM, VAR, ADD, SUB, MUL, DIV, POW, NEG,
    SIN, COS, TAN, ASIN, ACOS, ATAN,
    LN, LOG10, EXP, SQRT, ABS,
    FRESNELS, FRESNELC, ERF
};

struct Expr;
using ExprPtr = shared_ptr<Expr>;

struct Expr {
    Op op;
    double num = 0.0;
    string var;
    vector<ExprPtr> args;
    Expr(Op o) : op(o) {}
};

ExprPtr num(double v);
ExprPtr var(const string& name);
ExprPtr add(ExprPtr a, ExprPtr b);
ExprPtr sub(ExprPtr a, ExprPtr b);
ExprPtr mul(ExprPtr a, ExprPtr b);
ExprPtr divi(ExprPtr a, ExprPtr b);
ExprPtr powr(ExprPtr a, ExprPtr b);
ExprPtr neg(ExprPtr a);
ExprPtr fn(Op op, ExprPtr a);

bool isNum(const ExprPtr& e, double v = NAN);
bool contains(const ExprPtr& e, const string& v);
bool structEq(const ExprPtr& a, const ExprPtr& b);
bool linearCoeffs(const ExprPtr& e, const string& v, double& a, double& b);
string numToStr(double v);

ExprPtr simplify(const ExprPtr& e);
double evaluate(const ExprPtr& e, const map<string, double>& vals);

}
