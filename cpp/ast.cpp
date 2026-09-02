#include "ast.hpp"
#include <sstream>
#include <cmath>
#include <algorithm>
using namespace std;

namespace calc {

ExprPtr num(double v) { auto e = make_shared<Expr>(Op::NUM); e->num = v; return e; }
ExprPtr var(const string& n) { auto e = make_shared<Expr>(Op::VAR); e->var = n; return e; }

static ExprPtr bin(Op op, ExprPtr a, ExprPtr b) {
    auto e = make_shared<Expr>(op);
    e->args = {a, b};
    return e;
}

ExprPtr add(ExprPtr a, ExprPtr b) { return bin(Op::ADD, a, b); }
ExprPtr sub(ExprPtr a, ExprPtr b) { return bin(Op::SUB, a, b); }
ExprPtr mul(ExprPtr a, ExprPtr b) { return bin(Op::MUL, a, b); }
ExprPtr divi(ExprPtr a, ExprPtr b) { return bin(Op::DIV, a, b); }
ExprPtr powr(ExprPtr a, ExprPtr b) { return bin(Op::POW, a, b); }
ExprPtr neg(ExprPtr a) { auto e = make_shared<Expr>(Op::NEG); e->args = {a}; return e; }
ExprPtr fn(Op op, ExprPtr a) { auto e = make_shared<Expr>(op); e->args = {a}; return e; }

bool isNum(const ExprPtr& e, double v) {
    if (e->op != Op::NUM) return false;
    if (isnan(v)) return true;
    return abs(e->num - v) < 1e-12;
}

bool contains(const ExprPtr& e, const string& v) {
    if (e->op == Op::VAR) return e->var == v;
    for (auto& a : e->args) if (contains(a, v)) return true;
    return false;
}

bool structEq(const ExprPtr& a, const ExprPtr& b) {
    if (a->op != b->op) return false;
    if (a->op == Op::NUM) return abs(a->num - b->num) < 1e-12;
    if (a->op == Op::VAR) return a->var == b->var;
    if (a->args.size() != b->args.size()) return false;
    for (size_t i = 0; i < a->args.size(); i++)
        if (!structEq(a->args[i], b->args[i])) return false;
    return true;
}

bool linearCoeffs(const ExprPtr& e, const string& v, double& a, double& b) {
    switch (e->op) {
        case Op::NUM: a = 0; b = e->num; return true;
        case Op::VAR:
            if (e->var == v) { a = 1; b = 0; return true; }
            return false;
        case Op::NEG: {
            double a1, b1;
            if (!linearCoeffs(e->args[0], v, a1, b1)) return false;
            a = -a1; b = -b1; return true;
        }
        case Op::ADD: {
            double a1, b1, a2, b2;
            if (!linearCoeffs(e->args[0], v, a1, b1) ||
                !linearCoeffs(e->args[1], v, a2, b2)) return false;
            a = a1 + a2; b = b1 + b2; return true;
        }
        case Op::SUB: {
            double a1, b1, a2, b2;
            if (!linearCoeffs(e->args[0], v, a1, b1) ||
                !linearCoeffs(e->args[1], v, a2, b2)) return false;
            a = a1 - a2; b = b1 - b2; return true;
        }
        case Op::MUL: {
            double a1, b1, a2, b2;
            bool l1 = linearCoeffs(e->args[0], v, a1, b1);
            bool l2 = linearCoeffs(e->args[1], v, a2, b2);
            if (l1 && a1 == 0 && l2) { a = b1 * a2; b = b1 * b2; return true; }
            if (l2 && a2 == 0 && l1) { a = a1 * b2; b = b1 * b2; return true; }
            return false;
        }
        case Op::DIV: {
            double a1, b1, a2, b2;
            if (!linearCoeffs(e->args[0], v, a1, b1)) return false;
            if (!linearCoeffs(e->args[1], v, a2, b2) || a2 != 0) return false;
            if (b2 == 0) return false;
            a = a1 / b2; b = b1 / b2; return true;
        }
        default: return false;
    }
}

string numToStr(double v) {
    if (abs(v) < 1e-6) v = 0.0;
    if (abs(v - round(v)) < 1e-9) {
        long long r = (long long)llround(v);
        return to_string(r);
    }
    ostringstream ss;
    ss.precision(10);
    ss << v;
    return ss.str();
}

ExprPtr simplify(const ExprPtr& e) {
    if (e->op == Op::NUM || e->op == Op::VAR) return e;

    vector<ExprPtr> sargs;
    for (auto& a : e->args) sargs.push_back(simplify(a));

    switch (e->op) {
        case Op::NEG: {
            auto a = sargs[0];
            if (a->op == Op::NUM) return num(-a->num);
            if (a->op == Op::NEG) return a->args[0];
            return neg(a);
        }
        case Op::ADD: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(a, 0)) return b;
            if (isNum(b, 0)) return a;
            if (a->op == Op::NUM && b->op == Op::NUM) return num(a->num + b->num);
            if (structEq(a, b)) return simplify(mul(num(2), a));
            return add(a, b);
        }
        case Op::SUB: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(b, 0)) return a;
            if (a->op == Op::NUM && b->op == Op::NUM) return num(a->num - b->num);
            if (structEq(a, b)) return num(0);
            return sub(a, b);
        }
        case Op::MUL: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(a, 0) || isNum(b, 0)) return num(0);
            if (isNum(a, 1)) return b;
            if (isNum(b, 1)) return a;
            if (isNum(a, -1)) return simplify(neg(b));
            if (isNum(b, -1)) return simplify(neg(a));
            if (a->op == Op::NUM && b->op == Op::NUM) return num(a->num * b->num);
            return mul(a, b);
        }
        case Op::DIV: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(a, 0) && !isNum(b, 0)) return num(0);
            if (isNum(b, 1)) return a;
            if (a->op == Op::NUM && b->op == Op::NUM && b->num != 0)
                return num(a->num / b->num);
            if (structEq(a, b)) return num(1);
            return divi(a, b);
        }
        case Op::POW: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(b, 0)) return num(1);
            if (isNum(b, 1)) return a;
            if (isNum(a, 0)) return num(0);
            if (isNum(a, 1)) return num(1);
            if (a->op == Op::NUM && b->op == Op::NUM)
                return num(pow(a->num, b->num));
            return powr(a, b);
        }
        default: {

            auto a = sargs[0];
            if (a->op == Op::NUM) {
                double x = a->num, r = NAN;
                switch (e->op) {
                    case Op::SIN: r = sin(x); break;
                    case Op::COS: r = cos(x); break;
                    case Op::TAN: r = tan(x); break;
                    case Op::ASIN: r = asin(x); break;
                    case Op::ACOS: r = acos(x); break;
                    case Op::ATAN: r = atan(x); break;
                    case Op::LN: r = log(x); break;
                    case Op::LOG10: r = log10(x); break;
                    case Op::EXP: r = exp(x); break;
                    case Op::SQRT: r = sqrt(x); break;
                    case Op::ABS: r = abs(x); break;
                    default: break;
                }
                if (!isnan(r)) return num(r);
            }
            return fn(e->op, a);
        }
    }
}

double evaluate(const ExprPtr& e, const map<string, double>& vals) {
    switch (e->op) {
        case Op::NUM: return e->num;
        case Op::VAR: {
            auto it = vals.find(e->var);
            if (it == vals.end()) throw runtime_error("Unbound variable: " + e->var);
            return it->second;
        }
        case Op::ADD: return evaluate(e->args[0], vals) + evaluate(e->args[1], vals);
        case Op::SUB: return evaluate(e->args[0], vals) - evaluate(e->args[1], vals);
        case Op::MUL: return evaluate(e->args[0], vals) * evaluate(e->args[1], vals);
        case Op::DIV: return evaluate(e->args[0], vals) / evaluate(e->args[1], vals);
        case Op::POW: return pow(evaluate(e->args[0], vals), evaluate(e->args[1], vals));
        case Op::NEG: return -evaluate(e->args[0], vals);
        case Op::SIN: return sin(evaluate(e->args[0], vals));
        case Op::COS: return cos(evaluate(e->args[0], vals));
        case Op::TAN: return tan(evaluate(e->args[0], vals));
        case Op::ASIN: return asin(evaluate(e->args[0], vals));
        case Op::ACOS: return acos(evaluate(e->args[0], vals));
        case Op::ATAN: return atan(evaluate(e->args[0], vals));
        case Op::LN: return log(evaluate(e->args[0], vals));
        case Op::LOG10: return log10(evaluate(e->args[0], vals));
        case Op::EXP: return exp(evaluate(e->args[0], vals));
        case Op::SQRT: return sqrt(evaluate(e->args[0], vals));
        case Op::ABS: return abs(evaluate(e->args[0], vals));
        case Op::FRESNELS:
        case Op::FRESNELC:
        case Op::ERF:
            return NAN;
    }
    return NAN;
}

}
