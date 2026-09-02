#include "integration.hpp"
#include "differentiation.hpp"
#include "printing.hpp"
#include <cmath>
#include <functional>
using namespace std;

namespace calc {

static ExprPtr integrateLinearForm(const ExprPtr& e, const string& v, bool& ok) {
    double a, b;
    ok = false;
    if (e->op == Op::SIN || e->op == Op::COS || e->op == Op::EXP ||
        e->op == Op::LN || e->op == Op::TAN) {
        if (!linearCoeffs(e->args[0], v, a, b)) return nullptr;
        if (a == 0) return nullptr;
        ExprPtr u = e->args[0];
        ok = true;
        switch (e->op) {
            case Op::SIN: return divi(neg(fn(Op::COS, u)), num(a));
            case Op::COS: return divi(fn(Op::SIN, u), num(a));
            case Op::EXP: return divi(fn(Op::EXP, u), num(a));
            case Op::TAN: return divi(neg(fn(Op::LN, fn(Op::ABS, fn(Op::COS, u)))), num(a));
            case Op::LN:  return divi(sub(mul(u, fn(Op::LN, u)), u), num(a));
            default: break;
        }
    }
    return nullptr;
}

static ExprPtr tryIntegrate(const ExprPtr& e, const string& v, bool& ok) {
    ok = true;
    if (!contains(e, v)) return mul(e, var(v));

    switch (e->op) {
        case Op::VAR:
            if (e->var == v) return divi(powr(var(v), num(2)), num(2));
            ok = false; return nullptr;

        case Op::ADD: {
            bool ok1, ok2;
            ExprPtr ra = tryIntegrate(e->args[0], v, ok1);
            ExprPtr rb = tryIntegrate(e->args[1], v, ok2);
            if (!ok1 || !ok2) { ok = false; return nullptr; }
            return add(ra, rb);
        }

        case Op::SUB: {
            bool ok1, ok2;
            ExprPtr ra = tryIntegrate(e->args[0], v, ok1);
            ExprPtr rb = tryIntegrate(e->args[1], v, ok2);
            if (!ok1 || !ok2) { ok = false; return nullptr; }
            return sub(ra, rb);
        }

        case Op::NEG: {
            bool ok1;
            ExprPtr ra = tryIntegrate(e->args[0], v, ok1);
            if (!ok1) { ok = false; return nullptr; }
            return neg(ra);
        }

        case Op::MUL: {
            auto& a = e->args[0]; auto& b = e->args[1];
            if (!contains(a, v)) {
                bool ok1; ExprPtr rb = tryIntegrate(b, v, ok1);
                if (!ok1) { ok = false; return nullptr; }
                return mul(a, rb);
            }
            if (!contains(b, v)) {
                bool ok1; ExprPtr ra = tryIntegrate(a, v, ok1);
                if (!ok1) { ok = false; return nullptr; }
                return mul(b, ra);
            }
            ok = false; return nullptr;
        }

        case Op::DIV: {
            auto& a = e->args[0]; auto& b = e->args[1];
            if (!contains(b, v)) {
                bool ok1; ExprPtr ra = tryIntegrate(a, v, ok1);
                if (!ok1) { ok = false; return nullptr; }
                return divi(ra, b);
            }

            double la, lb;
            if (isNum(a, 1) && linearCoeffs(b, v, la, lb) && la != 0) {
                return divi(fn(Op::LN, fn(Op::ABS, b)), num(la));
            }
            ok = false; return nullptr;
        }

        case Op::POW: {
            auto& base = e->args[0]; auto& expo = e->args[1];
            double la, lb;
            if (expo->op == Op::NUM && linearCoeffs(base, v, la, lb) && la != 0) {
                double n = expo->num;
                if (abs(n + 1) < 1e-12) {
                    return divi(fn(Op::LN, fn(Op::ABS, base)), num(la));
                }
                return divi(powr(base, num(n + 1)), num(la * (n + 1)));
            }
            ok = false; return nullptr;
        }

        case Op::SIN: case Op::COS: case Op::EXP: case Op::LN: case Op::TAN: {
            bool sok; ExprPtr r = integrateLinearForm(e, v, sok);
            if (sok) return r;

            if (e->op == Op::SIN || e->op == Op::COS) {
                double coeff = 0.0;
                if (e->args[0]->op == Op::POW && e->args[0]->args[0]->op == Op::VAR && e->args[0]->args[0]->var == v && isNum(e->args[0]->args[1], 2)) {
                    coeff = 1.0;
                } else if (e->args[0]->op == Op::MUL) {
                    auto left = e->args[0]->args[0];
                    auto right = e->args[0]->args[1];
                    ExprPtr varPow = nullptr;
                    if (isNum(left)) { coeff = left->num; varPow = right; }
                    else if (isNum(right)) { coeff = right->num; varPow = left; }
                    if (varPow && varPow->op == Op::POW && varPow->args[0]->op == Op::VAR && varPow->args[0]->var == v && isNum(varPow->args[1], 2)) {
                    } else {
                        coeff = 0.0;
                    }
                }
                if (coeff > 0) {
                    ExprPtr coeffTerm = num(sqrt(M_PI / (2 * coeff)));
                    ExprPtr inner = mul(num(sqrt((2 * coeff) / M_PI)), var(v));
                    if (e->op == Op::SIN) {
                        return mul(coeffTerm, fn(Op::FRESNELS, inner));
                    } else {
                        return mul(coeffTerm, fn(Op::FRESNELC, inner));
                    }
                }
            }

            if (e->op == Op::EXP) {
                double coeff = 0.0;
                ExprPtr inner = e->args[0];
                if (inner->op == Op::NEG && inner->args[0]->op == Op::POW && inner->args[0]->args[0]->op == Op::VAR && inner->args[0]->args[0]->var == v && isNum(inner->args[0]->args[1], 2)) {
                    coeff = -1.0;
                } else if (inner->op == Op::MUL) {
                    auto left = inner->args[0];
                    auto right = inner->args[1];
                    ExprPtr varPow = nullptr;
                    if (isNum(left)) { coeff = left->num; varPow = right; }
                    else if (isNum(right)) { coeff = right->num; varPow = left; }
                    if (varPow && varPow->op == Op::POW && varPow->args[0]->op == Op::VAR && varPow->args[0]->var == v && isNum(varPow->args[1], 2)) {
                    } else {
                        coeff = 0.0;
                    }
                }
                if (coeff < 0) {
                    double a = -coeff;
                    ExprPtr coeffTerm = divi(num(sqrt(M_PI)), mul(num(2), num(sqrt(a))));
                    ExprPtr erfInner = mul(num(sqrt(a)), var(v));
                    return mul(coeffTerm, fn(Op::ERF, erfInner));
                }
            }

            ok = false; return nullptr;
        }

        default:
            ok = false; return nullptr;
    }
}

CalcResult integrate(const ExprPtr& e, const string& wrt) {
    bool ok;
    ExprPtr r = tryIntegrate(simplify(e), wrt, ok);
    if (!ok || !r) {
        return {false,
            "Cannot find a closed-form antiderivative (pattern not recognized). "
            "Try a numeric definite integral instead, or simplify the expression.",
            ""};
    }
    ExprPtr s = simplify(r);
    return {true, toString(s) + " + C", toLatex(s) + " + C"};
}

static double adaptiveSimpsonRec(
    const function<double(double)>& f,
    double a, double b,
    double fa, double fm, double fb,
    double whole, double eps, int depth)
{
    double m = (a + b) / 2.0;
    double h = (b - a) / 2.0;

    double m1 = (a + m) / 2.0;
    double m2 = (m + b) / 2.0;
    double fm1 = f(m1);
    double fm2 = f(m2);

    double left  = (h / 6.0) * (fa + 4.0 * fm1 + fm);
    double right = (h / 6.0) * (fm + 4.0 * fm2 + fb);
    double total = left + right;

    if (depth >= 20 || abs(total - whole) <= 15.0 * eps) {
        return total + (total - whole) / 15.0;
    }

    return adaptiveSimpsonRec(f, a, m, fa, fm1, fm, left, eps / 2.0, depth + 1) +
           adaptiveSimpsonRec(f, m, b, fm, fm2, fb, right, eps / 2.0, depth + 1);
}

static double adaptiveSimpson(const function<double(double)>& f,
                              double a, double b, double eps = 1e-10) {
    double fa = f(a);
    double fb = f(b);
    double fm = f((a + b) / 2.0);
    double whole = ((b - a) / 6.0) * (fa + 4.0 * fm + fb);
    return adaptiveSimpsonRec(f, a, b, fa, fm, fb, whole, eps, 0);
}

DefiniteResult definiteIntegral(const ExprPtr& e, const string& wrt,
                                double lower, double upper) {
    DefiniteResult res;
    res.ok = false;
    res.is_numeric = false;

    ExprPtr se = simplify(e);

    bool symOk;
    ExprPtr antideriv = tryIntegrate(se, wrt, symOk);
    if (symOk && antideriv) {
        ExprPtr sAntideriv = simplify(antideriv);
        map<string, double> upperVals, lowerVals;
        upperVals[wrt] = upper;
        lowerVals[wrt] = lower;
        try {
            double fUpper = evaluate(sAntideriv, upperVals);
            double fLower = evaluate(sAntideriv, lowerVals);
            double result = fUpper - fLower;
            if (isfinite(result)) {
                res.ok = true;
                res.value = result;
                res.is_numeric = false;
                res.text = numToStr(result);
                res.latex = numToStr(result);
                return res;
            }
        } catch (...) {

        }
    }

    try {
        auto f = [&](double x) -> double {
            map<string, double> vals;
            vals[wrt] = x;
            return evaluate(se, vals);
        };

        double result = adaptiveSimpson(f, lower, upper);
        if (isfinite(result)) {
            res.ok = true;
            res.value = result;
            res.is_numeric = true;
            res.text = numToStr(result) + " (numeric)";
            res.latex = numToStr(result);
            return res;
        }
    } catch (...) {}

    res.text = "Could not compute the definite integral.";
    return res;
}

MultiIntegralResult multiIntegral(const ExprPtr& e, const vector<IntegralSpec>& specs) {
    MultiIntegralResult res;
    res.ok = false;
    res.is_numeric = false;

    if (specs.empty()) {
        res.ok = true;
        res.text = toString(simplify(e));
        res.latex = toLatex(simplify(e));
        return res;
    }

    ExprPtr current = simplify(e);

    for (size_t i = 0; i < specs.size(); i++) {
        const auto& spec = specs[i];
        bool isDefinite = isfinite(spec.lower) && isfinite(spec.upper);

        if (isDefinite) {

            bool symOk;
            ExprPtr antideriv = tryIntegrate(current, spec.var, symOk);
            if (symOk && antideriv) {
                ExprPtr sAntideriv = simplify(antideriv);

                map<string, double> upperVals, lowerVals;
                upperVals[spec.var] = spec.upper;
                lowerVals[spec.var] = spec.lower;

                bool hasOtherVars = false;
                for (size_t j = i + 1; j < specs.size(); j++) {
                    if (contains(sAntideriv, specs[j].var)) {
                        hasOtherVars = true;
                        break;
                    }
                }

                if (!hasOtherVars) {

                    try {
                        double fUpper = evaluate(sAntideriv, upperVals);
                        double fLower = evaluate(sAntideriv, lowerVals);
                        double val = fUpper - fLower;
                        if (isfinite(val)) {
                            current = num(val);
                            continue;
                        }
                    } catch (...) {}
                }

            }

            auto f = [&](double x) -> double {
                map<string, double> vals;
                vals[spec.var] = x;
                return evaluate(current, vals);
            };

            try {
                double val = adaptiveSimpson(f, spec.lower, spec.upper);
                if (isfinite(val)) {
                    current = num(val);
                    res.is_numeric = true;
                    continue;
                }
            } catch (...) {}

            res.text = "Could not compute the iterated integral at dimension " +
                       to_string(i + 1) + " (variable: " + spec.var + ").";
            return res;

        } else {

            bool symOk;
            ExprPtr antideriv = tryIntegrate(current, spec.var, symOk);
            if (!symOk || !antideriv) {
                res.text = "Cannot find a closed-form antiderivative for variable " +
                           spec.var + " at integral " + to_string(i + 1) + ".";
                return res;
            }
            current = simplify(antideriv);
        }
    }

    current = simplify(current);
    res.ok = true;
    res.text = toString(current);
    res.latex = toLatex(current);

    bool allIndefinite = true;
    for (const auto& spec : specs) {
        if (isfinite(spec.lower) && isfinite(spec.upper)) {
            allIndefinite = false;
            break;
        }
    }
    if (allIndefinite) {
        res.text += " + C";
        res.latex += " + C";
    }

    return res;
}

}
