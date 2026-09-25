#include "laplace.hpp"
#include "differentiation.hpp"
#include "printing.hpp"
#include <cmath>
using namespace std;

namespace calc {

static bool extractLinearCoeff(const ExprPtr& e, const string& v,
                               double& a, double& b) {
    return linearCoeffs(e, v, a, b);
}

static ExprPtr forwardLaplace(const ExprPtr& e, const string& t,
                              const string& s, bool& ok, string& rule) {
    ok = true;

    if (!contains(e, t)) {

        rule = "L{c} = c/s";
        return divi(e, var(s));
    }

    switch (e->op) {
        case Op::VAR: {
            if (e->var == t) {

                rule = "L{t} = 1/s^2";
                return divi(num(1), powr(var(s), num(2)));
            }
            break;
        }
        case Op::ADD: {
            bool ok1, ok2;
            string r1, r2;
            ExprPtr L1 = forwardLaplace(e->args[0], t, s, ok1, r1);
            ExprPtr L2 = forwardLaplace(e->args[1], t, s, ok2, r2);
            if (ok1 && ok2) {
                rule = "Linearity";
                return add(L1, L2);
            }
            break;
        }
        case Op::SUB: {
            bool ok1, ok2;
            string r1, r2;
            ExprPtr L1 = forwardLaplace(e->args[0], t, s, ok1, r1);
            ExprPtr L2 = forwardLaplace(e->args[1], t, s, ok2, r2);
            if (ok1 && ok2) {
                rule = "Linearity";
                return sub(L1, L2);
            }
            break;
        }
        case Op::NEG: {
            bool ok1;
            string r1;
            ExprPtr L1 = forwardLaplace(e->args[0], t, s, ok1, r1);
            if (ok1) {
                rule = "Linearity";
                return neg(L1);
            }
            break;
        }
        case Op::MUL: {
            auto& a = e->args[0];
            auto& b = e->args[1];
            if (!contains(a, t)) {
                bool ok1;
                string r1;
                ExprPtr L1 = forwardLaplace(b, t, s, ok1, r1);
                if (ok1) {
                    rule = "Linearity";
                    return mul(a, L1);
                }
            } else if (!contains(b, t)) {
                bool ok1;
                string r1;
                ExprPtr L1 = forwardLaplace(a, t, s, ok1, r1);
                if (ok1) {
                    rule = "Linearity";
                    return mul(b, L1);
                }
            }

            auto checkShifting = [&](const ExprPtr& expTerm, const ExprPtr& fTerm) -> ExprPtr {
                if (expTerm->op == Op::EXP) {
                    double ca, cb;
                    if (extractLinearCoeff(expTerm->args[0], t, ca, cb)) {
                        bool fok;
                        string frule;
                        ExprPtr Fs = forwardLaplace(fTerm, t, s, fok, frule);
                        if (fok) {

                            auto replaceVar = [&](auto& self, const ExprPtr& node, const string& target, const ExprPtr& repl) -> ExprPtr {
                                if (node->op == Op::VAR && node->var == target) return repl;
                                if (node->op == Op::NUM || node->op == Op::VAR) return node;
                                auto newNode = make_shared<Expr>(node->op);
                                for (auto& arg : node->args) {
                                    newNode->args.push_back(self(self, arg, target, repl));
                                }
                                return newNode;
                            };
                            rule = "First Shifting Theorem: L{e^{at} f(t)} = F(s-a)";
                            ExprPtr shift = sub(var(s), num(ca));
                            return replaceVar(replaceVar, Fs, s, shift);
                        }
                    }
                }
                return nullptr;
            };

            ExprPtr res1 = checkShifting(a, b);
            if (res1) return res1;
            ExprPtr res2 = checkShifting(b, a);
            if (res2) return res2;

            auto checkTDerivative = [&](const ExprPtr& tTerm, const ExprPtr& fTerm) -> ExprPtr {
                if (tTerm->op == Op::VAR && tTerm->var == t) {
                    bool fok;
                    string frule;
                    ExprPtr Fs = forwardLaplace(fTerm, t, s, fok, frule);
                    if (fok) {
                        rule = "L{t f(t)} = -F'(s)";
                        return neg(differentiate(Fs, s));
                    }
                } else if (tTerm->op == Op::POW && tTerm->args[0]->op == Op::VAR &&
                           tTerm->args[0]->var == t && tTerm->args[1]->op == Op::NUM) {
                    double n = tTerm->args[1]->num;
                    if (abs(n - round(n)) < 1e-9 && n > 0 && n < 5) {
                        bool fok;
                        string frule;
                        ExprPtr Fs = forwardLaplace(fTerm, t, s, fok, frule);
                        if (fok) {
                            rule = "L{t^n f(t)} = (-1)^n F^(n)(s)";
                            ExprPtr deriv = differentiateOrder(Fs, s, (int)round(n));
                            if ((int)round(n) % 2 != 0) deriv = neg(deriv);
                            return deriv;
                        }
                    }
                }
                return nullptr;
            };

            res1 = checkTDerivative(a, b);
            if (res1) return res1;
            res2 = checkTDerivative(b, a);
            if (res2) return res2;

            break;
        }
        case Op::POW: {
            if (e->args[0]->op == Op::VAR && e->args[0]->var == t && e->args[1]->op == Op::NUM) {
                double n = e->args[1]->num;
                if (abs(n - round(n)) < 1e-9 && n >= 0) {

                    rule = "L{t^n} = n!/s^{n+1}";
                    double fact = 1;
                    for(int i = 2; i <= (int)round(n); i++) fact *= i;
                    return divi(num(fact), powr(var(s), num(n + 1)));
                }
            }
            break;
        }
        case Op::EXP: {
            double ca, cb;
            if (extractLinearCoeff(e->args[0], t, ca, cb)) {

                rule = "L{e^{at}} = 1/(s-a)";
                ExprPtr res = divi(num(1), sub(var(s), num(ca)));
                if (abs(cb) > 1e-9) {
                    res = mul(fn(Op::EXP, num(cb)), res);
                }
                return res;
            }
            break;
        }
        case Op::SIN: {
            double ca, cb;
            if (extractLinearCoeff(e->args[0], t, ca, cb) && abs(cb) < 1e-9) {

                rule = "L{sin(at)} = a/(s^2+a^2)";
                return divi(num(ca), add(powr(var(s), num(2)), num(ca * ca)));
            }
            break;
        }
        case Op::COS: {
            double ca, cb;
            if (extractLinearCoeff(e->args[0], t, ca, cb) && abs(cb) < 1e-9) {

                rule = "L{cos(at)} = s/(s^2+a^2)";
                return divi(var(s), add(powr(var(s), num(2)), num(ca * ca)));
            }
            break;
        }
        default:
            break;
    }

    ok = false;
    return nullptr;
}

static ExprPtr inverseLaplace(const ExprPtr& e, const string& s,
                              const string& t, bool& ok, string& rule) {
    ok = true;

    if (!contains(e, s)) {
        ok = false;
        return nullptr;
    }

    if (e->op == Op::ADD || e->op == Op::SUB) {
        bool ok1, ok2;
        string r1, r2;
        ExprPtr L1 = inverseLaplace(e->args[0], s, t, ok1, r1);
        ExprPtr L2 = inverseLaplace(e->args[1], s, t, ok2, r2);
        if (ok1 && ok2) {
            rule = "Linearity";
            return (e->op == Op::ADD) ? add(L1, L2) : sub(L1, L2);
        }
    }
    if (e->op == Op::NEG) {
        bool ok1;
        string r1;
        ExprPtr L1 = inverseLaplace(e->args[0], s, t, ok1, r1);
        if (ok1) {
            rule = "Linearity";
            return neg(L1);
        }
    }
    if (e->op == Op::MUL) {
        auto& a = e->args[0];
        auto& b = e->args[1];
        if (!contains(a, s)) {
            bool ok1;
            string r1;
            ExprPtr L1 = inverseLaplace(b, s, t, ok1, r1);
            if (ok1) {
                rule = "Linearity";
                return mul(a, L1);
            }
        } else if (!contains(b, s)) {
            bool ok1;
            string r1;
            ExprPtr L1 = inverseLaplace(a, s, t, ok1, r1);
            if (ok1) {
                rule = "Linearity";
                return mul(b, L1);
            }
        }
    }

    if (e->op == Op::DIV) {
        auto& num_expr = e->args[0];
        auto& den = e->args[1];

        if (!contains(num_expr, s)) {

            if (den->op == Op::VAR && den->var == s) {
                rule = "L^{-1}{1/s} = 1";
                return num_expr;
            }

            double da, db;
            if (extractLinearCoeff(den, s, da, db) && abs(da) > 1e-9) {
                rule = "L^{-1}{1/(s-a)} = e^{at}";
                double a = -db / da;
                ExprPtr expTerm = fn(Op::EXP, mul(num(a), var(t)));
                return mul(divi(num_expr, num(da)), expTerm);
            }

            if (den->op == Op::POW && den->args[0]->op == Op::VAR &&
                den->args[0]->var == s && den->args[1]->op == Op::NUM) {
                double n = den->args[1]->num;
                if (abs(n - round(n)) < 1e-9 && n > 1) {
                    rule = "L^{-1}{1/s^n} = t^{n-1}/(n-1)!";
                    double fact = 1;
                    for (int i = 2; i < (int)round(n); i++) fact *= i;
                    return divi(mul(num_expr, powr(var(t), num(n - 1))), num(fact));
                }
            }

            if (den->op == Op::ADD) {

                 auto isSquareOf = [&](const ExprPtr& node, const string& target_var) -> bool {
                     if (node->op == Op::POW && node->args[0]->op == Op::VAR &&
                         node->args[0]->var == target_var &&
                         node->args[1]->op == Op::NUM && abs(node->args[1]->num - 2.0) < 1e-9) {
                         return true;
                     }
                     return false;
                 };

                 bool s_squared_first = isSquareOf(den->args[0], s) && !contains(den->args[1], s);
                 bool s_squared_second = isSquareOf(den->args[1], s) && !contains(den->args[0], s);

                 if (s_squared_first || s_squared_second) {
                     ExprPtr const_term = s_squared_first ? den->args[1] : den->args[0];
                     if (const_term->op == Op::NUM && const_term->num > 0) {
                         double a = sqrt(const_term->num);
                         rule = "L^{-1}{1/(s^2+a^2)} = (1/a)sin(at)";
                         return mul(divi(num_expr, num(a)), fn(Op::SIN, mul(num(a), var(t))));
                     }
                 }
            }
        } else if (num_expr->op == Op::VAR && num_expr->var == s) {

             if (den->op == Op::ADD) {
                 auto isSquareOf = [&](const ExprPtr& node, const string& target_var) -> bool {
                     if (node->op == Op::POW && node->args[0]->op == Op::VAR &&
                         node->args[0]->var == target_var &&
                         node->args[1]->op == Op::NUM && abs(node->args[1]->num - 2.0) < 1e-9) {
                         return true;
                     }
                     return false;
                 };

                 bool s_squared_first = isSquareOf(den->args[0], s) && !contains(den->args[1], s);
                 bool s_squared_second = isSquareOf(den->args[1], s) && !contains(den->args[0], s);

                 if (s_squared_first || s_squared_second) {
                     ExprPtr const_term = s_squared_first ? den->args[1] : den->args[0];
                     if (const_term->op == Op::NUM && const_term->num > 0) {
                         double a = sqrt(const_term->num);
                         rule = "L^{-1}{s/(s^2+a^2)} = cos(at)";
                         return fn(Op::COS, mul(num(a), var(t)));
                     }
                 }
            }
        }
    }

    ok = false;
    return nullptr;
}

LaplaceResult laplaceTransform(const ExprPtr& expr, const string& t_var,
                               const string& s_var) {
    LaplaceResult res;
    res.ok = false;

    try {
        ExprPtr se = simplify(expr);
        bool ok;
        string rule;
        ExprPtr Fs = forwardLaplace(se, t_var, s_var, ok, rule);

        if (ok && Fs) {
            Fs = simplify(Fs);
            res.ok = true;
            res.text = toString(Fs);
            res.latex = toLatex(Fs);
            res.method = rule;
        } else {
            res.message = "Could not find a Laplace transform using the built-in rules. "
                          "It supports basic combinations of polynomials, exponentials, and sines/cosines.";
        }
    } catch (const exception& ex) {
        res.message = string("Error computing Laplace transform: ") + ex.what();
    }

    return res;
}

LaplaceResult inverseLaplaceTransform(const ExprPtr& expr, const string& s_var,
                                      const string& t_var) {
    LaplaceResult res;
    res.ok = false;

    try {
        ExprPtr se = simplify(expr);
        bool ok;
        string rule;
        ExprPtr ft = inverseLaplace(se, s_var, t_var, ok, rule);

        if (ok && ft) {
            ft = simplify(ft);
            res.ok = true;
            res.text = toString(ft);
            res.latex = toLatex(ft);
            res.method = rule;
        } else {
            res.message = "Could not find an Inverse Laplace transform using the built-in rules. "
                          "It currently supports simple rational functions and basic patterns.";
        }
    } catch (const exception& ex) {
        res.message = string("Error computing Inverse Laplace transform: ") + ex.what();
    }

    return res;
}

}
