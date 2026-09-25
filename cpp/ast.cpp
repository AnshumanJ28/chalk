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

bool quadraticCoeffs(const ExprPtr& e, const string& v, double& a, double& b, double& c) {
    // Try to express e as a*v^2 + b*v + c
    // Strategy: evaluate at v=0, v=1, v=-1 numerically on the simplified expression
    // and check if the expression is actually quadratic by also checking v=2.
    // This is a heuristic approach that works for polynomial expressions.
    
    // First, make sure the expression only contains v as a variable
    // and is polynomial (no trig, exp, ln, etc. involving v)
    
    // Try symbolic approach: collect terms
    switch (e->op) {
        case Op::NUM: a = 0; b = 0; c = e->num; return true;
        case Op::VAR:
            if (e->var == v) { a = 0; b = 1; c = 0; return true; }
            return false; // other variable
        case Op::NEG: {
            double a1, b1, c1;
            if (!quadraticCoeffs(e->args[0], v, a1, b1, c1)) return false;
            a = -a1; b = -b1; c = -c1; return true;
        }
        case Op::ADD: {
            double a1, b1, c1, a2, b2, c2;
            if (!quadraticCoeffs(e->args[0], v, a1, b1, c1) ||
                !quadraticCoeffs(e->args[1], v, a2, b2, c2)) return false;
            a = a1 + a2; b = b1 + b2; c = c1 + c2; return true;
        }
        case Op::SUB: {
            double a1, b1, c1, a2, b2, c2;
            if (!quadraticCoeffs(e->args[0], v, a1, b1, c1) ||
                !quadraticCoeffs(e->args[1], v, a2, b2, c2)) return false;
            a = a1 - a2; b = b1 - b2; c = c1 - c2; return true;
        }
        case Op::MUL: {
            double a1, b1, c1, a2, b2, c2;
            bool l1 = quadraticCoeffs(e->args[0], v, a1, b1, c1);
            bool l2 = quadraticCoeffs(e->args[1], v, a2, b2, c2);
            if (!l1 || !l2) return false;
            // (a1*v^2 + b1*v + c1) * (a2*v^2 + b2*v + c2)
            // Result is degree <= 2 only if both are degree <= 1 or one is constant
            if (a1 == 0 && a2 == 0) {
                // linear * linear = quadratic
                a = b1 * b2;
                b = b1 * c2 + c1 * b2;
                c = c1 * c2;
                return true;
            }
            if (a1 == 0 && b1 == 0) {
                // constant * quadratic
                a = c1 * a2; b = c1 * b2; c = c1 * c2;
                return true;
            }
            if (a2 == 0 && b2 == 0) {
                // quadratic * constant
                a = a1 * c2; b = b1 * c2; c = c1 * c2;
                return true;
            }
            return false; // degree > 2
        }
        case Op::POW: {
            if (e->args[1]->op == Op::NUM && abs(e->args[1]->num - 2.0) < 1e-12) {
                // base^2 where base is linear in v
                double la, lb;
                if (linearCoeffs(e->args[0], v, la, lb)) {
                    a = la * la;
                    b = 2 * la * lb;
                    c = lb * lb;
                    return true;
                }
            }
            if (e->args[1]->op == Op::NUM && abs(e->args[1]->num - 1.0) < 1e-12) {
                return quadraticCoeffs(e->args[0], v, a, b, c);
            }
            if (!contains(e, v)) { a = 0; b = 0; /* evaluate later */ return false; }
            return false;
        }
        case Op::DIV: {
            double a1, b1, c1;
            if (!quadraticCoeffs(e->args[0], v, a1, b1, c1)) return false;
            // denominator must be constant
            double la2, lb2;
            if (!linearCoeffs(e->args[1], v, la2, lb2) || la2 != 0) return false;
            if (lb2 == 0) return false;
            a = a1 / lb2; b = b1 / lb2; c = c1 / lb2;
            return true;
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
            // (a/b) * c = (a*c)/b if simplification helps
            if (a->op == Op::DIV) {
                return simplify(divi(simplify(mul(a->args[0], b)), a->args[1]));
            }
            // a * (b/c) = (a*b)/c
            if (b->op == Op::DIV) {
                return simplify(divi(simplify(mul(a, b->args[0])), b->args[1]));
            }
            return mul(a, b);
        }
        case Op::DIV: {
            auto a = sargs[0], b = sargs[1];
            if (isNum(a, 0) && !isNum(b, 0)) return num(0);
            if (isNum(b, 1)) return a;
            if (a->op == Op::NUM && b->op == Op::NUM && b->num != 0)
                return num(a->num / b->num);
            if (structEq(a, b)) return num(1);
            
            // Negation cancellation
            if (b->op == Op::NEG && structEq(a, b->args[0])) return num(-1);
            if (a->op == Op::NEG && structEq(a->args[0], b)) return num(-1);
            if (a->op == Op::NEG && b->op == Op::NEG) return simplify(divi(a->args[0], b->args[0]));

            // Nested division
            if (a->op == Op::DIV && b->op == Op::DIV) {
                return simplify(divi(mul(a->args[0], b->args[1]), mul(a->args[1], b->args[0])));
            }
            if (a->op == Op::DIV) {
                return simplify(divi(a->args[0], mul(a->args[1], b)));
            }
            if (b->op == Op::DIV) {
                return simplify(divi(mul(a, b->args[1]), b->args[0]));
            }

            // Cancel common factors: f / (c*f) = 1/c, (c*f) / f = c
            // a / (c*a) -> 1/c
            if (b->op == Op::MUL) {
                if (structEq(a, b->args[0])) return simplify(divi(num(1), b->args[1]));
                if (structEq(a, b->args[1])) return simplify(divi(num(1), b->args[0]));
            }
            // (c*a) / b where b matches one factor -> other factor
            if (a->op == Op::MUL) {
                if (structEq(a->args[0], b)) return a->args[1];
                if (structEq(a->args[1], b)) return a->args[0];
            }
            // (c1*f) / (c2*f) -> c1/c2
            if (a->op == Op::MUL && b->op == Op::MUL) {
                if (structEq(a->args[0], b->args[0]))
                    return simplify(divi(a->args[1], b->args[1]));
                if (structEq(a->args[1], b->args[1]))
                    return simplify(divi(a->args[0], b->args[0]));
                if (structEq(a->args[0], b->args[1]))
                    return simplify(divi(a->args[1], b->args[0]));
                if (structEq(a->args[1], b->args[0]))
                    return simplify(divi(a->args[0], b->args[1]));
            }
            // (c*f) / (c*g) -> f/g (cancel numeric coefficient)
            if (a->op == Op::MUL && b->op == Op::MUL &&
                a->args[0]->op == Op::NUM && b->args[0]->op == Op::NUM) {
                return simplify(divi(mul(num(a->args[0]->num / b->args[0]->num), a->args[1]), b->args[1]));
            }
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
                    case Op::SINH: r = sinh(x); break;
                    case Op::COSH: r = cosh(x); break;
                    case Op::TANH: r = tanh(x); break;
                    case Op::ATANH: r = atanh(x); break;
                    case Op::SEC: r = 1.0/cos(x); break;
                    case Op::CSC: r = 1.0/sin(x); break;
                    case Op::COT: r = cos(x)/sin(x); break;
                    default: break;
                }
                if (!isnan(r) && isfinite(r)) return num(r);
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
        case Op::SINH: return sinh(evaluate(e->args[0], vals));
        case Op::COSH: return cosh(evaluate(e->args[0], vals));
        case Op::TANH: return tanh(evaluate(e->args[0], vals));
        case Op::ATANH: return atanh(evaluate(e->args[0], vals));
        case Op::SEC: return 1.0/cos(evaluate(e->args[0], vals));
        case Op::CSC: return 1.0/sin(evaluate(e->args[0], vals));
        case Op::COT: { double x = evaluate(e->args[0], vals); return cos(x)/sin(x); }
        case Op::FRESNELS:
        case Op::FRESNELC:
        case Op::ERF:
            return NAN;
    }
    return NAN;
}

ExprPtr substitute(const ExprPtr& e, const string& v, const ExprPtr& replacement) {
    if (e->op == Op::NUM) return e;
    if (e->op == Op::VAR) {
        if (e->var == v) return replacement;
        return e;
    }
    vector<ExprPtr> newArgs;
    for (auto& a : e->args) newArgs.push_back(substitute(a, v, replacement));
    if (e->args.size() == 1) {
        if (e->op == Op::NEG) return neg(newArgs[0]);
        return fn(e->op, newArgs[0]);
    }
    switch (e->op) {
        case Op::ADD: return add(newArgs[0], newArgs[1]);
        case Op::SUB: return sub(newArgs[0], newArgs[1]);
        case Op::MUL: return mul(newArgs[0], newArgs[1]);
        case Op::DIV: return divi(newArgs[0], newArgs[1]);
        case Op::POW: return powr(newArgs[0], newArgs[1]);
        default: return e;
    }
}

static string opNameForSig(Op op) {
    switch (op) {
        case Op::NUM: return "NUM";
        case Op::VAR: return "VAR";
        case Op::ADD: return "ADD"; case Op::SUB: return "SUB";
        case Op::MUL: return "MUL"; case Op::DIV: return "DIV";
        case Op::POW: return "POW"; case Op::NEG: return "NEG";
        case Op::SIN: return "SIN"; case Op::COS: return "COS"; case Op::TAN: return "TAN";
        case Op::ASIN: return "ASIN"; case Op::ACOS: return "ACOS"; case Op::ATAN: return "ATAN";
        case Op::SINH: return "SINH"; case Op::COSH: return "COSH"; case Op::TANH: return "TANH";
        case Op::ATANH: return "ATANH"; case Op::SEC: return "SEC"; case Op::CSC: return "CSC";
        case Op::COT: return "COT"; case Op::LN: return "LN"; case Op::LOG10: return "LOG10";
        case Op::EXP: return "EXP"; case Op::SQRT: return "SQRT"; case Op::ABS: return "ABS";
        case Op::FRESNELS: return "FRESNELS"; case Op::FRESNELC: return "FRESNELC"; case Op::ERF: return "ERF";
        default: return "?";
    }
}

string getSignature(const ExprPtr& e, const string& v) {
    if (e->op == Op::NUM) return "NUM";
    if (e->op == Op::VAR) return (e->var == v) ? "VAR" : "NUM";

    if (e->args.size() == 1) {
        return opNameForSig(e->op) + "(" + getSignature(e->args[0], v) + ")";
    }

    if (e->args.size() == 2) {
        string s1 = getSignature(e->args[0], v);
        string s2 = getSignature(e->args[1], v);
        if (e->op == Op::ADD || e->op == Op::MUL) {
            if (s1 > s2) swap(s1, s2);
        }
        return opNameForSig(e->op) + "(" + s1 + "," + s2 + ")";
    }
    return "UNKNOWN";
}

}
