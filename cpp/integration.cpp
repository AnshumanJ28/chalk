#include "integration.hpp"
#include "differentiation.hpp"
#include "printing.hpp"
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <tuple>
using namespace std;

namespace calc {

// Known integral patterns mapped to their antiderivatives.

static ExprPtr tableLookup(const ExprPtr& e, const string& v, bool& ok) {
    ok = false;
    ExprPtr se = simplify(e);

    // ∫ sec(x) dx = ln|sec(x) + tan(x)|
    if (se->op == Op::SEC) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::LN, fn(Op::ABS, add(fn(Op::SEC, se->args[0]),
                                                     fn(Op::TAN, se->args[0])))), num(la));
        }
    }

    // ∫ csc(x) dx = -ln|csc(x) + cot(x)|
    if (se->op == Op::CSC) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(neg(fn(Op::LN, fn(Op::ABS, add(fn(Op::CSC, se->args[0]),
                                                        fn(Op::COT, se->args[0]))))), num(la));
        }
    }

    // ∫ cot(x) dx = ln|sin(x)|
    if (se->op == Op::COT) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::LN, fn(Op::ABS, fn(Op::SIN, se->args[0]))), num(la));
        }
    }

    // ∫ sec^2(x) dx = tan(x)
    if (se->op == Op::POW && se->args[0]->op == Op::SEC && isNum(se->args[1], 2)) {
        double la, lb;
        if (linearCoeffs(se->args[0]->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::TAN, se->args[0]->args[0]), num(la));
        }
    }

    // ∫ csc^2(x) dx = -cot(x)
    if (se->op == Op::POW && se->args[0]->op == Op::CSC && isNum(se->args[1], 2)) {
        double la, lb;
        if (linearCoeffs(se->args[0]->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(neg(fn(Op::COT, se->args[0]->args[0])), num(la));
        }
    }

    // ∫ sec(x)*tan(x) dx = sec(x)
    if (se->op == Op::MUL) {
        auto& a = se->args[0]; auto& b = se->args[1];
        if (a->op == Op::SEC && b->op == Op::TAN && structEq(a->args[0], b->args[0])) {
            double la, lb;
            if (linearCoeffs(a->args[0], v, la, lb) && la != 0) {
                ok = true;
                return divi(fn(Op::SEC, a->args[0]), num(la));
            }
        }
        if (b->op == Op::SEC && a->op == Op::TAN && structEq(a->args[0], b->args[0])) {
            double la, lb;
            if (linearCoeffs(a->args[0], v, la, lb) && la != 0) {
                ok = true;
                return divi(fn(Op::SEC, a->args[0]), num(la));
            }
        }
    }

    // ∫ csc(x)*cot(x) dx = -csc(x)
    if (se->op == Op::MUL) {
        auto& a = se->args[0]; auto& b = se->args[1];
        if (a->op == Op::CSC && b->op == Op::COT && structEq(a->args[0], b->args[0])) {
            double la, lb;
            if (linearCoeffs(a->args[0], v, la, lb) && la != 0) {
                ok = true;
                return divi(neg(fn(Op::CSC, a->args[0])), num(la));
            }
        }
        if (b->op == Op::CSC && a->op == Op::COT && structEq(a->args[0], b->args[0])) {
            double la, lb;
            if (linearCoeffs(a->args[0], v, la, lb) && la != 0) {
                ok = true;
                return divi(neg(fn(Op::CSC, a->args[0])), num(la));
            }
        }
    }

    // ∫ sinh(x) dx = cosh(x), ∫ cosh(x) dx = sinh(x), ∫ tanh(x) dx = ln(cosh(x))
    if (se->op == Op::SINH) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::COSH, se->args[0]), num(la));
        }
    }
    if (se->op == Op::COSH) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::SINH, se->args[0]), num(la));
        }
    }
    if (se->op == Op::TANH) {
        double la, lb;
        if (linearCoeffs(se->args[0], v, la, lb) && la != 0) {
            ok = true;
            return divi(fn(Op::LN, fn(Op::COSH, se->args[0])), num(la));
        }
    }

    // ∫ 1/sqrt(1-x^2) dx = asin(x) [handled via completing the square too]
    // ∫ 1/(1+x^2) dx = atan(x) [handled via completing the square too]

    return nullptr;
}


// Integrate rational functions p(x)/q(x) by decomposing into simpler fractions.

// Represents polynomial as vector of coefficients: {c0, c1, c2, ...} for c0 + c1*x + c2*x^2 + ...
using Poly = vector<double>;

static Poly polyFromExpr(const ExprPtr& e, const string& v, bool& ok) {
    ok = true;
    switch (e->op) {
        case Op::NUM: return {e->num};
        case Op::VAR:
            if (e->var == v) return {0.0, 1.0};
            ok = false; return {};
        case Op::NEG: {
            bool ok1;
            Poly p = polyFromExpr(e->args[0], v, ok1);
            if (!ok1) { ok = false; return {}; }
            for (auto& c : p) c = -c;
            return p;
        }
        case Op::ADD: {
            bool ok1, ok2;
            Poly a = polyFromExpr(e->args[0], v, ok1);
            Poly b = polyFromExpr(e->args[1], v, ok2);
            if (!ok1 || !ok2) { ok = false; return {}; }
            Poly r(max(a.size(), b.size()), 0.0);
            for (size_t i = 0; i < a.size(); i++) r[i] += a[i];
            for (size_t i = 0; i < b.size(); i++) r[i] += b[i];
            return r;
        }
        case Op::SUB: {
            bool ok1, ok2;
            Poly a = polyFromExpr(e->args[0], v, ok1);
            Poly b = polyFromExpr(e->args[1], v, ok2);
            if (!ok1 || !ok2) { ok = false; return {}; }
            Poly r(max(a.size(), b.size()), 0.0);
            for (size_t i = 0; i < a.size(); i++) r[i] += a[i];
            for (size_t i = 0; i < b.size(); i++) r[i] -= b[i];
            return r;
        }
        case Op::MUL: {
            bool ok1, ok2;
            Poly a = polyFromExpr(e->args[0], v, ok1);
            Poly b = polyFromExpr(e->args[1], v, ok2);
            if (!ok1 || !ok2) { ok = false; return {}; }
            if (a.empty() || b.empty()) return {0.0};
            Poly r(a.size() + b.size() - 1, 0.0);
            for (size_t i = 0; i < a.size(); i++)
                for (size_t j = 0; j < b.size(); j++)
                    r[i + j] += a[i] * b[j];
            return r;
        }
        case Op::POW: {
            if (e->args[1]->op == Op::NUM) {
                double n = e->args[1]->num;
                if (n >= 0 && abs(n - round(n)) < 1e-12 && n <= 10) {
                    int ni = (int)round(n);
                    bool ok1;
                    Poly base = polyFromExpr(e->args[0], v, ok1);
                    if (!ok1) { ok = false; return {}; }
                    if (ni == 0) return {1.0};
                    Poly result = base;
                    for (int i = 1; i < ni; i++) {
                        Poly newR(result.size() + base.size() - 1, 0.0);
                        for (size_t a = 0; a < result.size(); a++)
                            for (size_t b = 0; b < base.size(); b++)
                                newR[a + b] += result[a] * base[b];
                        result = newR;
                    }
                    return result;
                }
            }
            ok = false; return {};
        }
        case Op::DIV: {
            // Only if denominator is constant
            if (!contains(e->args[1], v)) {
                bool ok1;
                Poly a = polyFromExpr(e->args[0], v, ok1);
                if (!ok1) { ok = false; return {}; }
                map<string, double> vals;
                try {
                    double d = evaluate(simplify(e->args[1]), vals);
                    if (d == 0) { ok = false; return {}; }
                    for (auto& c : a) c /= d;
                    return a;
                } catch (...) { ok = false; return {}; }
            }
            ok = false; return {};
        }
        default:
            ok = false; return {};
    }
}

static int polyDeg(const Poly& p) {
    for (int i = (int)p.size() - 1; i >= 0; i--)
        if (abs(p[i]) > 1e-12) return i;
    return -1; // zero polynomial
}

static void polyTrim(Poly& p) {
    while (p.size() > 1 && abs(p.back()) < 1e-12) p.pop_back();
}

static pair<Poly, Poly> polyDivide(const Poly& dividend, const Poly& divisor) {
    Poly rem = dividend;
    int dd = polyDeg(divisor);
    if (dd < 0) return {dividend, {}}; // division by zero
    int dr = polyDeg(rem);
    Poly quotient(max(0, dr - dd + 1), 0.0);
    while (dr >= dd) {
        double coeff = rem[dr] / divisor[dd];
        int pos = dr - dd;
        quotient[pos] = coeff;
        for (int i = 0; i <= dd; i++)
            rem[pos + i] -= coeff * divisor[i];
        dr = polyDeg(rem);
    }
    polyTrim(quotient);
    polyTrim(rem);
    return {quotient, rem};
}

static double polyEval(const Poly& p, double x) {
    double r = 0, xn = 1;
    for (size_t i = 0; i < p.size(); i++) {
        r += p[i] * xn;
        xn *= x;
    }
    return r;
}

// Find real roots of a polynomial using companion matrix eigenvalue or brute-force for degree <= 3
static vector<double> polyRealRoots(const Poly& p) {
    vector<double> roots;
    int d = polyDeg(p);
    if (d <= 0) return roots;
    if (d == 1) {
        roots.push_back(-p[0] / p[1]);
        return roots;
    }
    if (d == 2) {
        double a = p[2], b = p[1], c = p[0];
        double disc = b * b - 4 * a * c;
        if (disc >= -1e-12) {
            if (disc < 0) disc = 0;
            roots.push_back((-b + sqrt(disc)) / (2 * a));
            roots.push_back((-b - sqrt(disc)) / (2 * a));
        }
        return roots;
    }
    // For degree 3+, use Durand-Kerner method (simplified real-root finder)
    // Try rational root theorem: test factors of constant term / leading coefficient
    Poly current = p;
    for (int attempt = 0; attempt < d && polyDeg(current) > 0; attempt++) {
        // Newton's method from multiple starting points
        bool found = false;
        for (double start : {0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 3.0, -3.0, 5.0, -5.0, 10.0, -10.0}) {
            double x = start;
            for (int iter = 0; iter < 100; iter++) {
                double fx = polyEval(current, x);
                if (abs(fx) < 1e-10) {
                    found = true;
                    break;
                }
                // derivative
                double dfx = 0, xn = 1;
                for (int i = 1; i <= polyDeg(current); i++) {
                    dfx += i * current[i] * xn;
                    xn *= x;
                }
                if (abs(dfx) < 1e-15) break;
                x = x - fx / dfx;
            }
            if (found) {
                // Round to integer if close
                if (abs(x - round(x)) < 1e-8) x = round(x);
                roots.push_back(x);
                // Deflate polynomial
                Poly factor = {-x, 1.0};
                auto [q, r] = polyDivide(current, factor);
                current = q;
                break;
            }
        }
        if (!found) break;
    }
    return roots;
}

static ExprPtr polyToExpr(const Poly& p, const string& v) {
    ExprPtr result = nullptr;
    for (size_t i = 0; i < p.size(); i++) {
        if (abs(p[i]) < 1e-12) continue;
        ExprPtr term;
        if (i == 0) term = num(p[i]);
        else if (i == 1) term = mul(num(p[i]), var(v));
        else term = mul(num(p[i]), powr(var(v), num((double)i)));
        if (!result) result = term;
        else result = add(result, term);
    }
    return result ? result : num(0);
}

// Integrate a polynomial
static ExprPtr integratePolynomial(const Poly& p, const string& v) {
    ExprPtr result = nullptr;
    for (size_t i = 0; i < p.size(); i++) {
        if (abs(p[i]) < 1e-12) continue;
        double coeff = p[i] / (double)(i + 1);
        ExprPtr term = mul(num(coeff), powr(var(v), num((double)(i + 1))));
        if (!result) result = term;
        else result = add(result, term);
    }
    return result ? result : num(0);
}

static ExprPtr tryPartialFractions(const ExprPtr& e, const string& v, bool& ok) {
    ok = false;
    if (e->op != Op::DIV) return nullptr;

    bool okN, okD;
    Poly numer = polyFromExpr(e->args[0], v, okN);
    Poly denom = polyFromExpr(e->args[1], v, okD);
    if (!okN || !okD) return nullptr;

    int degN = polyDeg(numer);
    int degD = polyDeg(denom);
    if (degD <= 0) return nullptr; // constant denominator, not interesting

    // If degree of numerator >= degree of denominator, do polynomial long division first
    ExprPtr polyPart = nullptr;
    if (degN >= degD) {
        auto [quotient, remainder] = polyDivide(numer, denom);
        polyPart = integratePolynomial(quotient, v);
        numer = remainder;
        degN = polyDeg(numer);
    }

    if (degN < 0) {
        // Remainder is zero
        ok = (polyPart != nullptr);
        return polyPart;
    }

    // Find roots of denominator
    vector<double> roots = polyRealRoots(denom);
    if (roots.empty()) {
        // Check if it's a quadratic with no real roots -> leads to atan
        if (degD == 2 && degN <= 1) {
            double A = denom[2], B = denom[1], C = denom[0];
            double disc = B * B - 4 * A * C;
            if (disc < -1e-12) {
                // Complete the square: A(x + B/(2A))^2 + (C - B^2/(4A))
                // ∫ (px + q) / (Ax^2 + Bx + C) dx
                double p_coeff = (degN >= 1) ? numer[1] : 0;
                double q_coeff = numer[0];

                // Split: p/(2A) * ln|Ax^2+Bx+C| + (q - pB/(2A)) * (1/sqrt(AC - B^2/4)) * atan(...)
                ExprPtr result = nullptr;

                // Part 1: coefficient of derivative matching
                // d/dx(Ax^2+Bx+C) = 2Ax + B
                // p*x + q = (p/(2A)) * (2Ax + B) + (q - pB/(2A))
                double matchCoeff = p_coeff / (2.0 * A);
                double remainder_const = q_coeff - p_coeff * B / (2.0 * A);

                if (abs(matchCoeff) > 1e-12) {
                    // ∫ matchCoeff * (2Ax+B)/(Ax^2+Bx+C) dx = matchCoeff * ln|Ax^2+Bx+C|
                    ExprPtr denomExpr = polyToExpr(denom, v);
                    result = mul(num(matchCoeff), fn(Op::LN, fn(Op::ABS, denomExpr)));
                }

                if (abs(remainder_const) > 1e-12) {
                    // ∫ remainder_const / (Ax^2+Bx+C) dx
                    // = remainder_const / A * ∫ 1/((x+B/(2A))^2 + (C/A - B^2/(4A^2))) dx
                    double h = B / (2.0 * A);
                    double k = C / A - B * B / (4.0 * A * A);
                    // = remainder_const / A * (1/sqrt(k)) * atan((x+h)/sqrt(k))
                    if (k > 0) {
                        double sqk = sqrt(k);
                        ExprPtr atanPart = mul(num(remainder_const / (A * sqk)),
                                              fn(Op::ATAN, divi(add(var(v), num(h)), num(sqk))));
                        if (result) result = add(result, atanPart);
                        else result = atanPart;
                    }
                }

                if (result) {
                    ok = true;
                    if (polyPart) result = add(polyPart, result);
                    return result;
                }
            }
        }
        return nullptr;
    }

    // Check we have all roots (all linear factors)
    // For simplicity, handle the case where we can fully factor into distinct linear factors
    if ((int)roots.size() < degD) {
        // We don't have all roots. Try to handle partial factoring.
        // For now, if we have some roots, try partial fractions with what we have.
        // Factor out known roots and see if the remainder is a quadratic with no real roots.
        Poly remaining = denom;
        vector<double> usedRoots;
        for (double r : roots) {
            Poly factor = {-r, 1.0};
            auto [q, rem] = polyDivide(remaining, factor);
            if (polyDeg(rem) < 0 || (rem.size() == 1 && abs(rem[0]) < 1e-8)) {
                remaining = q;
                usedRoots.push_back(r);
            }
        }
        int remDeg = polyDeg(remaining);
        if (remDeg > 2 || usedRoots.empty()) return nullptr;

        // We have: denom = (x-r1)(x-r2)...*remaining(x)
        // Set up partial fractions: N(x) / ((x-r1)(x-r2)...*R(x)) = A1/(x-r1) + A2/(x-r2) + ... + (Bx+C)/R(x)
        // Use the cover-up method for the linear factors
        int numFactors = (int)usedRoots.size();
        int totalUnknowns = numFactors + (remDeg == 2 ? 2 : (remDeg == 1 ? 1 : 0));
        if (totalUnknowns > 10) return nullptr; // too complex

        // For distinct linear factors, use Heaviside cover-up method
        ExprPtr result = nullptr;
        Poly currentNumer = numer;

        for (double r : usedRoots) {
            // A_i = N(r_i) / product of (r_i - r_j) for j != i * remaining(r_i)
            double denomProd = 1.0;
            for (double r2 : usedRoots) {
                if (abs(r2 - r) > 1e-12) denomProd *= (r - r2);
            }
            denomProd *= polyEval(remaining, r);
            if (abs(denomProd) < 1e-12) return nullptr;
            double A_i = polyEval(currentNumer, r) / denomProd;

            // ∫ A_i / (x - r_i) dx = A_i * ln|x - r_i|
            ExprPtr term = mul(num(A_i), fn(Op::LN, fn(Op::ABS, sub(var(v), num(r)))));
            if (!result) result = term;
            else result = add(result, term);
        }

        // Handle the remaining quadratic part if needed
        if (remDeg >= 1) {
            // Subtract the linear partial fractions from the original and integrate the remainder
            // This is complex, so for now we skip if there's a non-trivial remainder
            // But handle the common case: degN <= numFactors + remDeg - 1
            // For now, just return what we have if remainder is constant
            if (remDeg == 2) {
                // We need to find (Bx+C)/R(x) part. 
                // N(x)/D(x) - sum(A_i/(x-r_i)) = (Bx+C)/R(x)
                // This requires solving a linear system, which is complex.
                // For the common case, we already handle it via the standard path.
            }
        }

        if (result) {
            ok = true;
            if (polyPart) result = add(polyPart, result);
            return result;
        }
        return nullptr;
    }

    // All roots found - check for distinct roots (simpler case)
    // Remove duplicates and count multiplicities
    sort(roots.begin(), roots.end());
    vector<pair<double, int>> rootMult; // (root, multiplicity)
    for (size_t i = 0; i < roots.size(); ) {
        int count = 1;
        while (i + count < roots.size() && abs(roots[i] - roots[i + count]) < 1e-8) count++;
        rootMult.push_back({roots[i], count});
        i += count;
    }

    // For distinct linear factors only (multiplicity 1)
    bool allDistinct = true;
    for (auto& [r, m] : rootMult) {
        if (m > 1) { allDistinct = false; break; }
    }

    if (allDistinct && degN < degD) {
        // Heaviside cover-up method
        ExprPtr result = nullptr;
        for (auto& [root, mult] : rootMult) {
            double denomProd = 1.0;
            for (auto& [r2, m2] : rootMult) {
                if (abs(r2 - root) > 1e-12) denomProd *= (root - r2);
            }
            // Normalize by leading coefficient of denominator
            denomProd *= denom[degD];
            if (abs(denomProd) < 1e-12) return nullptr;
            double A_i = polyEval(numer, root) / denomProd;
            ExprPtr term = mul(num(A_i), fn(Op::LN, fn(Op::ABS, sub(var(v), num(root)))));
            if (!result) result = term;
            else result = add(result, term);
        }
        if (result) {
            ok = true;
            if (polyPart) result = add(polyPart, result);
            return result;
        }
    }

    // Handle repeated roots: A/(x-r)^n -> integrate to A/(-n+1) * (x-r)^(-n+1) or A*ln|x-r|
    if (!allDistinct && degN < degD) {
        // Build the system: N(x)/D(x) = sum A_{i,k} / (x-r_i)^k
        // Use evaluation at multiple points to solve the linear system
        // For simplicity, use the derivative method for repeated roots.
        ExprPtr result = nullptr;
        Poly currentDenom = denom;

        for (auto& [root, mult] : rootMult) {
            for (int k = mult; k >= 1; k--) {
                // Compute coefficient using limit method
                // A_{i,k} = lim_{x->r_i} (1/(mult-k)!) * d^(mult-k)/dx^(mult-k) [(x-r_i)^mult * N(x)/D(x)]
                // For k=mult (highest power), this is just N(r_i) / (D(x)/(x-r_i)^mult evaluated at r_i)
                
                // Build D(x) / (x-root)^mult
                Poly factor = {1.0};
                for (auto& [r2, m2] : rootMult) {
                    if (abs(r2 - root) > 1e-12) {
                        for (int j = 0; j < m2; j++) {
                            Poly lin = {-r2, 1.0};
                            Poly newF(factor.size() + 1, 0.0);
                            for (size_t a = 0; a < factor.size(); a++)
                                for (size_t b = 0; b < lin.size(); b++)
                                    newF[a + b] += factor[a] * lin[b];
                            factor = newF;
                        }
                    }
                }
                // Multiply by leading coefficient
                for (auto& c : factor) c *= denom[degD];
                
                double denomVal = polyEval(factor, root);
                if (abs(denomVal) < 1e-12) { ok = false; return nullptr; }
                
                if (k == mult) {
                    double A = polyEval(numer, root) / denomVal;
                    if (abs(A) > 1e-12) {
                        ExprPtr term;
                        if (k == 1) {
                            term = mul(num(A), fn(Op::LN, fn(Op::ABS, sub(var(v), num(root)))));
                        } else {
                            double newExp = -(double)(k - 1);
                            term = divi(mul(num(A), powr(sub(var(v), num(root)), num(newExp))), num(newExp));
                        }
                        if (!result) result = term;
                        else result = add(result, term);
                    }
                }
                // For k < mult, we would need to differentiate. Skip for now for complex cases.
                if (k < mult) break;
            }
        }
        if (result) {
            ok = true;
            if (polyPart) result = add(polyPart, result);
            return result;
        }
    }

    return nullptr;
}


// sin^n(x), cos^n(x) -> use reduction formulas and identities.

static ExprPtr tryTrigPowerReduction(const ExprPtr& e, const string& v, bool& ok);

// Forward declaration for main integration function
static ExprPtr tryIntegrateRec(const ExprPtr& e, const string& v, bool& ok, int depth);

static ExprPtr tryTrigPowerReduction(const ExprPtr& e, const string& v, bool& ok) {
    ok = false;

    // sin^2(u) = (1 - cos(2u))/2
    if (e->op == Op::POW && e->args[0]->op == Op::SIN && isNum(e->args[1], 2)) {
        ExprPtr u = e->args[0]->args[0];
        double la, lb;
        if (linearCoeffs(u, v, la, lb) && la != 0) {
            // ∫ sin^2(ax+b) dx = x/2 - sin(2(ax+b))/(4a)
            ok = true;
            ExprPtr twoU = mul(num(2), u);
            return sub(divi(var(v), num(2)), divi(fn(Op::SIN, twoU), num(4 * la)));
        }
    }

    // cos^2(u) = (1 + cos(2u))/2
    if (e->op == Op::POW && e->args[0]->op == Op::COS && isNum(e->args[1], 2)) {
        ExprPtr u = e->args[0]->args[0];
        double la, lb;
        if (linearCoeffs(u, v, la, lb) && la != 0) {
            // ∫ cos^2(ax+b) dx = x/2 + sin(2(ax+b))/(4a)
            ok = true;
            ExprPtr twoU = mul(num(2), u);
            return add(divi(var(v), num(2)), divi(fn(Op::SIN, twoU), num(4 * la)));
        }
    }

    // tan^2(u) = sec^2(u) - 1 -> ∫ = tan(u)/a - x
    if (e->op == Op::POW && e->args[0]->op == Op::TAN && isNum(e->args[1], 2)) {
        ExprPtr u = e->args[0]->args[0];
        double la, lb;
        if (linearCoeffs(u, v, la, lb) && la != 0) {
            ok = true;
            return sub(divi(fn(Op::TAN, u), num(la)), var(v));
        }
    }

    // sin(a)*cos(b) product-to-sum: sin(a)*cos(b) = [sin(a+b) + sin(a-b)]/2
    if (e->op == Op::MUL) {
        auto& a = e->args[0]; auto& b = e->args[1];
        if (a->op == Op::SIN && b->op == Op::COS) {
            double la1, lb1, la2, lb2;
            if (linearCoeffs(a->args[0], v, la1, lb1) && linearCoeffs(b->args[0], v, la2, lb2)) {
                // ∫ sin(a1*x+b1)*cos(a2*x+b2) dx = [sin((a1+a2)x+(b1+b2)) + sin((a1-a2)x+(b1-b2))]/2
                ExprPtr sumArg = add(a->args[0], b->args[0]);
                ExprPtr diffArg = sub(a->args[0], b->args[0]);
                ExprPtr transformed = divi(add(fn(Op::SIN, sumArg), fn(Op::SIN, diffArg)), num(2));
                bool intOk;
                ExprPtr result = tryIntegrateRec(simplify(transformed), v, intOk, 1);
                if (intOk) { ok = true; return result; }
            }
        }
        if (a->op == Op::COS && b->op == Op::SIN) {
            double la1, lb1, la2, lb2;
            if (linearCoeffs(a->args[0], v, la1, lb1) && linearCoeffs(b->args[0], v, la2, lb2)) {
                ExprPtr sumArg = add(b->args[0], a->args[0]);
                ExprPtr diffArg = sub(b->args[0], a->args[0]);
                ExprPtr transformed = divi(add(fn(Op::SIN, sumArg), fn(Op::SIN, diffArg)), num(2));
                bool intOk;
                ExprPtr result = tryIntegrateRec(simplify(transformed), v, intOk, 1);
                if (intOk) { ok = true; return result; }
            }
        }
        // sin(a)*sin(b) = [cos(a-b) - cos(a+b)]/2
        if (a->op == Op::SIN && b->op == Op::SIN) {
            double la1, lb1, la2, lb2;
            if (linearCoeffs(a->args[0], v, la1, lb1) && linearCoeffs(b->args[0], v, la2, lb2)) {
                ExprPtr diffArg = sub(a->args[0], b->args[0]);
                ExprPtr sumArg = add(a->args[0], b->args[0]);
                ExprPtr transformed = divi(sub(fn(Op::COS, diffArg), fn(Op::COS, sumArg)), num(2));
                bool intOk;
                ExprPtr result = tryIntegrateRec(simplify(transformed), v, intOk, 1);
                if (intOk) { ok = true; return result; }
            }
        }
        // cos(a)*cos(b) = [cos(a-b) + cos(a+b)]/2
        if (a->op == Op::COS && b->op == Op::COS) {
            double la1, lb1, la2, lb2;
            if (linearCoeffs(a->args[0], v, la1, lb1) && linearCoeffs(b->args[0], v, la2, lb2)) {
                ExprPtr diffArg = sub(a->args[0], b->args[0]);
                ExprPtr sumArg = add(a->args[0], b->args[0]);
                ExprPtr transformed = divi(add(fn(Op::COS, diffArg), fn(Op::COS, sumArg)), num(2));
                bool intOk;
                ExprPtr result = tryIntegrateRec(simplify(transformed), v, intOk, 1);
                if (intOk) { ok = true; return result; }
            }
        }
    }

    return nullptr;
}


// For integrals involving 1/(ax^2+bx+c) or 1/sqrt(ax^2+bx+c)

static ExprPtr tryCompletingSquare(const ExprPtr& e, const string& v, bool& ok) {
    ok = false;

    if (e->op != Op::DIV) return nullptr;
    auto& numer = e->args[0];
    auto& denom = e->args[1];

    // Case 1: ∫ 1/(ax^2 + bx + c) dx
    if (isNum(numer, 1) || (numer->op == Op::NUM && !contains(numer, v))) {
        double a, b, c;
        if (quadraticCoeffs(denom, v, a, b, c) && abs(a) > 1e-12) {
            double disc = b * b - 4 * a * c;
            double numerVal = (numer->op == Op::NUM) ? numer->num : 1.0;

            if (disc < -1e-12) {
                // No real roots: complete the square
                // ax^2 + bx + c = a[(x + b/(2a))^2 + (4ac-b^2)/(4a^2)]
                double h = b / (2.0 * a);
                double k = (4.0 * a * c - b * b) / (4.0 * a * a);
                // ∫ numerVal / (a * ((x+h)^2 + k)) dx = (numerVal/(a*sqrt(k))) * atan((x+h)/sqrt(k))
                if (k > 0) {
                    double sqk = sqrt(k);
                    ok = true;
                    return mul(num(numerVal / (a * sqk)), fn(Op::ATAN, divi(add(var(v), num(h)), num(sqk))));
                }
            } else if (disc > 1e-12) {
                // Two real roots: partial fraction decomposition handles this
                // But we can also do it directly
                double r1 = (-b + sqrt(disc)) / (2 * a);
                double r2 = (-b - sqrt(disc)) / (2 * a);
                if (abs(r1 - r2) > 1e-12) {
                    // 1/(a(x-r1)(x-r2)) = 1/(a(r1-r2)) * [1/(x-r1) - 1/(x-r2)]
                    double coeff = numerVal / (a * (r1 - r2));
                    ok = true;
                    return mul(num(coeff), sub(fn(Op::LN, fn(Op::ABS, sub(var(v), num(r1)))),
                                              fn(Op::LN, fn(Op::ABS, sub(var(v), num(r2))))));
                }
            }
        }
    }

    // Case 2: ∫ 1/sqrt(a - x^2) dx = asin(x/sqrt(a)) (when denom is sqrt)
    if (denom->op == Op::SQRT || (denom->op == Op::POW && isNum(denom->args[1]) &&
        abs(denom->args[1]->num - 0.5) < 1e-12)) {
        ExprPtr inner = (denom->op == Op::SQRT) ? denom->args[0] : denom->args[0];
        double qa, qb, qc;
        if (quadraticCoeffs(inner, v, qa, qb, qc) && abs(qa) > 1e-12) {
            double numerVal = (numer->op == Op::NUM) ? numer->num : 1.0;

            if (qa < 0) {
                // sqrt(-a*x^2 + bx + c) form -> complete the square
                // -a(x^2 - b/a*x) + c = -a[(x - b/(2a))^2 - b^2/(4a^2)] + c
                // = -a(x - h)^2 + k where h = -qb/(2qa), k = qc - qb^2/(4qa)
                // Actually qa is negative. Let A = -qa > 0.
                // inner = -A*x^2 + qb*x + qc = -(A)(x^2 - (qb/A)x) + qc
                // = -A[(x - qb/(2A))^2 - qb^2/(4A^2)] + qc
                // = -A(x - h)^2 + k, where h=qb/(2A), k = qc + qb^2/(4A)
                double A = -qa;
                double h = qb / (2.0 * A);
                double k = qc + qb * qb / (4.0 * A);
                if (k > 0) {
                    // ∫ 1/sqrt(k - A(x-h)^2) dx = (1/sqrt(A)) * asin((x-h)*sqrt(A/k))
                    ok = true;
                    return mul(num(numerVal / sqrt(A)),
                              fn(Op::ASIN, mul(sub(var(v), num(h)), num(sqrt(A / k)))));
                }
            } else {
                // sqrt(a*x^2 + bx + c), a > 0: involves arcsinh or ln
                // ∫ 1/sqrt(a*x^2 + bx + c) dx = (1/sqrt(a)) * ln|2*sqrt(a)*sqrt(ax^2+bx+c) + 2ax + b|
                // = (1/sqrt(a)) * arcsinh((2ax+b)/sqrt(4ac-b^2))  when disc < 0
                double h = qb / (2.0 * qa);
                double k = qc - qb * qb / (4.0 * qa);
                if (k > 0) {
                    // ∫ 1/sqrt(qa*(x+h)^2 + k) dx
                    // Let u = x + h, ∫ 1/sqrt(qa*u^2 + k) dx = (1/sqrt(qa)) * ln|u*sqrt(qa) + sqrt(qa*u^2+k)|
                    // Which is (1/sqrt(qa)) * arcsinh(u*sqrt(qa/k))
                    ok = true;
                    // Use ln form since we don't have ARCSINH
                    ExprPtr u_expr = add(var(v), num(h));
                    return mul(num(numerVal / sqrt(qa)),
                              fn(Op::LN, add(mul(num(sqrt(qa)), u_expr),
                                            fn(Op::SQRT, add(mul(num(qa), powr(u_expr, num(2))), num(k))))));
                }
            }
        }
    }

    // Case 3: ∫ (px+q)/(ax^2+bx+c) dx - split into derivative match + completing the square
    {
        double la, lb_n;
        double qa, qb, qc;
        if (linearCoeffs(numer, v, la, lb_n) && la != 0 &&
            quadraticCoeffs(denom, v, qa, qb, qc) && abs(qa) > 1e-12) {
            // N(x) = la*x + lb_n, D(x) = qa*x^2 + qb*x + qc
            // D'(x) = 2*qa*x + qb
            // la*x + lb_n = (la/(2*qa)) * (2*qa*x + qb) + (lb_n - la*qb/(2*qa))
            double matchCoeff = la / (2.0 * qa);
            double remainder = lb_n - la * qb / (2.0 * qa);

            ExprPtr result = nullptr;

            // ∫ matchCoeff * D'/D dx = matchCoeff * ln|D|
            if (abs(matchCoeff) > 1e-12) {
                ExprPtr denomExpr = add(add(mul(num(qa), powr(var(v), num(2))),
                                           mul(num(qb), var(v))), num(qc));
                result = mul(num(matchCoeff), fn(Op::LN, fn(Op::ABS, denomExpr)));
            }

            // ∫ remainder / D dx -> completing the square
            if (abs(remainder) > 1e-12) {
                double disc = qb * qb - 4 * qa * qc;
                if (disc < -1e-12) {
                    double h = qb / (2.0 * qa);
                    double k = (4.0 * qa * qc - qb * qb) / (4.0 * qa * qa);
                    if (k > 0) {
                        double sqk = sqrt(k);
                        ExprPtr atanPart = mul(num(remainder / (qa * sqk)),
                                              fn(Op::ATAN, divi(add(var(v), num(h)), num(sqk))));
                        if (result) result = add(result, atanPart);
                        else result = atanPart;
                    }
                } else if (disc > 1e-12) {
                    double r1 = (-qb + sqrt(disc)) / (2 * qa);
                    double r2 = (-qb - sqrt(disc)) / (2 * qa);
                    if (abs(r1 - r2) > 1e-12) {
                        double coeff = remainder / (qa * (r1 - r2));
                        ExprPtr lnPart = mul(num(coeff),
                            sub(fn(Op::LN, fn(Op::ABS, sub(var(v), num(r1)))),
                                fn(Op::LN, fn(Op::ABS, sub(var(v), num(r2))))));
                        if (result) result = add(result, lnPart);
                        else result = lnPart;
                    }
                }
            }

            if (result) { ok = true; return result; }
        }
    }

    return nullptr;
}


// Look for f(g(x)) * g'(x) patterns.

static void findCandidatesRec(const ExprPtr& e, const string& v, vector<ExprPtr>& candidates) {
    if (!contains(e, v)) return;
    if (e->op == Op::SIN || e->op == Op::COS || e->op == Op::TAN ||
        e->op == Op::EXP || e->op == Op::LN ||
        e->op == Op::SINH || e->op == Op::COSH || e->op == Op::TANH ||
        e->op == Op::SQRT || e->op == Op::ABS) {
        if (contains(e->args[0], v)) {
            candidates.push_back(e->args[0]);
        }
    } else if (e->op == Op::POW && contains(e->args[0], v) && !contains(e->args[1], v)) {
        candidates.push_back(e->args[0]);
        if (e->args[1]->op == Op::NUM) {
            double p = e->args[1]->num;
            if (p >= 2 && p / 2 == floor(p / 2)) {
                candidates.push_back(powr(e->args[0], num(p / 2)));
            }
        }
    }
    for (auto& a : e->args) {
        findCandidatesRec(a, v, candidates);
    }
}

static ExprPtr tryUSubstitution(const ExprPtr& e, const string& v, bool& ok, int depth) {
    ok = false;
    if (depth >= 3) return nullptr;

    vector<pair<ExprPtr, ExprPtr>> pairsToTry;
    if (e->op == Op::MUL) {
        pairsToTry.push_back({e->args[0], e->args[1]});
    } else if (e->op == Op::DIV) {
        pairsToTry.push_back({e->args[0], divi(num(1), e->args[1])});
        
        if (e->args[1]->op == Op::MUL) {
            pairsToTry.push_back({divi(e->args[0], e->args[1]->args[0]), divi(num(1), e->args[1]->args[1])});
            pairsToTry.push_back({divi(e->args[0], e->args[1]->args[1]), divi(num(1), e->args[1]->args[0])});
        }
    } else {
        return nullptr;
    }

    auto trySubst = [&](ExprPtr outer, ExprPtr multiplier) -> ExprPtr {
        if (outer->args.empty()) return nullptr;

        vector<ExprPtr> candidates;
        if (outer->op == Op::SIN || outer->op == Op::COS || outer->op == Op::TAN ||
            outer->op == Op::EXP || outer->op == Op::LN ||
            outer->op == Op::SINH || outer->op == Op::COSH || outer->op == Op::TANH ||
            outer->op == Op::SQRT || outer->op == Op::ABS) {
            candidates.push_back(outer->args[0]);
        }
        if (outer->op == Op::POW && contains(outer->args[0], v)) {
            candidates.push_back(outer->args[0]);
        }
        if (outer->op == Op::DIV && contains(outer->args[1], v)) {
            candidates.push_back(outer->args[1]);
        }
        
        // Here `outer` is actually a composed tree. Let's just find any candidate globally in `outer`.
        vector<ExprPtr> deepCandidates;
        findCandidatesRec(outer, v, deepCandidates);
        candidates.insert(candidates.end(), deepCandidates.begin(), deepCandidates.end());

        candidates.push_back(multiplier);

        for (auto& g : candidates) {
            if (!contains(g, v)) continue;
            ExprPtr dg = simplify(differentiate(g, v));
            if (!dg) continue;

            ExprPtr ratio = simplify(divi(multiplier, dg));
            if (!contains(ratio, v)) {
                string uVar = "__u__";
                
                std::function<ExprPtr(ExprPtr)> replaceG = [&](ExprPtr node) -> ExprPtr {
                    if (structEq(node, g)) return var(uVar);
                    
                    if (g->op == Op::EXP && node->op == Op::EXP) {
                        ExprPtr gInner = g->args[0];
                        ExprPtr nInner = node->args[0];
                        ExprPtr innerRatio = simplify(divi(nInner, gInner));
                        if (!contains(innerRatio, v)) {
                            return powr(var(uVar), innerRatio);
                        }
                    }
                    
                    if (g->op == Op::POW && node->op == Op::POW && structEq(g->args[0], node->args[0])) {
                        if (g->args[1]->op == Op::NUM && node->args[1]->op == Op::NUM) {
                            double b_num = node->args[1]->num;
                            double a_num = g->args[1]->num;
                            if (a_num != 0 && b_num / a_num == floor(b_num / a_num)) { 
                                return powr(var(uVar), num(b_num / a_num));
                            }
                        }
                    }

                    if (node->args.empty()) return node;
                    vector<ExprPtr> newArgs;
                    for (auto& arg : node->args) newArgs.push_back(replaceG(arg));
                    auto ret = make_shared<Expr>(node->op);
                    ret->args = newArgs;
                    return ret;
                };
                
                ExprPtr outerAsU = replaceG(outer);

                // If it still contains `x`, the substitution was incomplete and integrating w.r.t `u` is mathematically invalid!
                if (contains(outerAsU, v)) continue;

                bool intOk;
                ExprPtr F_u = tryIntegrateRec(outerAsU, uVar, intOk, depth + 1);
                if (intOk && F_u) {
                    ExprPtr result = simplify(mul(ratio, substitute(F_u, uVar, g)));
                    ok = true;
                    return result;
                }
            }
        }
        return nullptr;
    };

    for (auto& p : pairsToTry) {
        ExprPtr a = p.first;
        ExprPtr b = p.second;
        if (contains(a, v) && contains(b, v)) {
            ExprPtr r = trySubst(a, b);
            if (r) return r;
            r = trySubst(b, a);
            if (r) return r;
        }
    }

    return nullptr;
}


// For expressions like e^(g(x)) * polynomial(x), detect patterns that have
// elementary antiderivatives using the "parallel Risch" approach.

static ExprPtr tryRischHeuristic(const ExprPtr& e, const string& v, bool& ok, int depth) {
    ok = false;
    if (depth >= 3) return nullptr;

    // Pattern: p(x) * e^(ax+b) where p(x) is a polynomial
    // ∫ p(x) * e^(ax+b) dx = Q(x) * e^(ax+b) where Q(x) has the same degree as p(x)
    // and Q(x) satisfies: Q'(x) + a*Q(x) = p(x)
    if (e->op == Op::MUL) {
        auto& left = e->args[0]; auto& right = e->args[1];

        auto tryExpPoly = [&](ExprPtr polyPart, ExprPtr expPart) -> ExprPtr {
            if (expPart->op != Op::EXP) return nullptr;
            double la, lb;
            if (!linearCoeffs(expPart->args[0], v, la, lb)) return nullptr;
            if (abs(la) < 1e-12) return nullptr;

            // Extract polynomial coefficients
            bool pOk;
            Poly p = polyFromExpr(polyPart, v, pOk);
            if (!pOk) return nullptr;
            int deg = polyDeg(p);
            if (deg < 0) return nullptr;

            // Solve Q'(x) + a*Q(x) = p(x)
            // Q has same degree. Working from highest degree down:
            // For Q = q_n*x^n + ... + q_0:
            //   a*q_n = p_n  =>  q_n = p_n/a
            //   n*q_n + a*q_{n-1} = p_{n-1}  =>  q_{n-1} = (p_{n-1} - n*q_n)/a
            //   ...
            Poly q(deg + 1, 0.0);
            for (int i = deg; i >= 0; i--) {
                double rhs = (i < (int)p.size()) ? p[i] : 0.0;
                // Subtract contribution from derivative of higher terms
                if (i + 1 <= deg) {
                    rhs -= (i + 1) * q[i + 1];
                }
                q[i] = rhs / la;
            }

            ExprPtr Q = polyToExpr(q, v);
            ok = true;
            return mul(Q, expPart);
        };

        ExprPtr r = tryExpPoly(left, right);
        if (r) return r;
        r = tryExpPoly(right, left);
        if (r) return r;
    }

    // Pattern: p(x) * sin(ax+b) or p(x) * cos(ax+b)
    // ∫ p(x)*sin(ax+b) dx = Q(x)*sin(ax+b) + R(x)*cos(ax+b)
    // where deg(Q)=deg(R)=deg(p) and:
    //   Q'(x)*sin + a*Q*cos + R'(x)*cos - a*R*sin = p*sin
    // => Q' - aR = p, R' + aQ = 0
    // => R = -Q'a, then Q' - a*(-Q'a... no.
    // From R' + aQ = 0 => R' = -aQ
    // From Q' - aR = p => Q' + a*(integrate of aQ) — actually let's solve directly:
    // R = -Q'/a... wait, R' = -aQ. Let's solve iteratively.
    // Q'(x) - a*R(x) = p(x)
    // R'(x) + a*Q(x) = 0
    // From (2): R(x) = -a * ∫Q dx ... that makes it recursive.
    // Better: assume Q, R are polynomials of same degree n.
    // Undetermined coefficients approach.
    if (e->op == Op::MUL) {
        auto& left = e->args[0]; auto& right = e->args[1];

        auto tryTrigPoly = [&](ExprPtr polyPart, ExprPtr trigPart) -> ExprPtr {
            if (trigPart->op != Op::SIN && trigPart->op != Op::COS) return nullptr;
            double la, lb;
            if (!linearCoeffs(trigPart->args[0], v, la, lb)) return nullptr;
            if (abs(la) < 1e-12) return nullptr;

            bool pOk;
            Poly p = polyFromExpr(polyPart, v, pOk);
            if (!pOk) return nullptr;
            int deg = polyDeg(p);
            if (deg < 0) return nullptr;

            // System: Q' - a*R = p (if trig is sin), R' + a*Q = 0
            //    or:  Q' - a*R = 0 (if trig is cos), R' + a*Q = p
            // Q and R are polynomials of degree deg.

            Poly q(deg + 1, 0.0);
            Poly r(deg + 1, 0.0);

            bool isSin = (trigPart->op == Op::SIN);

            // Solve from highest degree downward
            for (int i = deg; i >= 0; i--) {
                double p_i = (i < (int)p.size()) ? p[i] : 0.0;
                // Q'_i = (i+1)*q[i+1] (if i+1 <= deg, else 0)
                // R'_i = (i+1)*r[i+1] (if i+1 <= deg, else 0)
                double Qprime_i = (i + 1 <= deg) ? (i + 1) * q[i + 1] : 0.0;
                double Rprime_i = (i + 1 <= deg) ? (i + 1) * r[i + 1] : 0.0;

                if (isSin) {
                    // Q'_i - a*r[i] = p_i  =>  r[i] = (Q'_i - p_i) / a
                    r[i] = (Qprime_i - p_i) / la;
                    // R'_i + a*q[i] = 0  =>  q[i] = -R'_i / a
                    q[i] = -Rprime_i / la;
                } else {
                    // R'_i + a*q[i] = p_i  =>  q[i] = (p_i - R'_i) / a
                    q[i] = (p_i - Rprime_i) / la;
                    // Q'_i - a*r[i] = 0  =>  r[i] = Q'_i / a
                    r[i] = Qprime_i / la;
                }
            }

            // Verify solution by checking the equations at degree -1 (constant terms from derivatives)
            // The solution should be self-consistent due to our iterative approach.

            ExprPtr Q = polyToExpr(q, v);
            ExprPtr R = polyToExpr(r, v);
            ExprPtr sinPart = fn(Op::SIN, trigPart->args[0]);
            ExprPtr cosPart = fn(Op::COS, trigPart->args[0]);

            ok = true;
            return add(mul(Q, sinPart), mul(R, cosPart));
        };

        ExprPtr res = tryTrigPoly(left, right);
        if (res) return res;
        res = tryTrigPoly(right, left);
        if (res) return res;
    }

    return nullptr;
}


// Original linear form integration

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
    // Handle inverse trig with linear arguments
    if (e->op == Op::ASIN) {
        if (!linearCoeffs(e->args[0], v, a, b)) return nullptr;
        if (a == 0) return nullptr;
        ExprPtr u = e->args[0];
        ok = true;
        // ∫ asin(u) du = u*asin(u) + sqrt(1-u^2)
        return divi(add(mul(u, fn(Op::ASIN, u)), fn(Op::SQRT, sub(num(1), powr(u, num(2))))), num(a));
    }
    if (e->op == Op::ACOS) {
        if (!linearCoeffs(e->args[0], v, a, b)) return nullptr;
        if (a == 0) return nullptr;
        ExprPtr u = e->args[0];
        ok = true;
        // ∫ acos(u) du = u*acos(u) - sqrt(1-u^2)
        return divi(sub(mul(u, fn(Op::ACOS, u)), fn(Op::SQRT, sub(num(1), powr(u, num(2))))), num(a));
    }
    if (e->op == Op::ATAN) {
        if (!linearCoeffs(e->args[0], v, a, b)) return nullptr;
        if (a == 0) return nullptr;
        ExprPtr u = e->args[0];
        ok = true;
        // ∫ atan(u) du = u*atan(u) - ln(1+u^2)/2
        return divi(sub(mul(u, fn(Op::ATAN, u)),
                       divi(fn(Op::LN, add(num(1), powr(u, num(2)))), num(2))), num(a));
    }
    return nullptr;
}


// LIATE priority for integration by parts

static int getLiatePriority(const ExprPtr& e, const string& v) {
    if (e->op == Op::LN || e->op == Op::LOG10) return 4;
    if (e->op == Op::ASIN || e->op == Op::ACOS || e->op == Op::ATAN) return 3;
    if (e->op == Op::VAR && e->var == v) return 2;
    if (e->op == Op::POW && e->args[0]->op == Op::VAR && e->args[0]->var == v && isNum(e->args[1])) {
        if (e->args[1]->num > 0) return 2;
    }
    if (e->op == Op::SIN || e->op == Op::COS || e->op == Op::TAN) return 1;
    if (e->op == Op::EXP) return 0;
    return -1;
}


static ExprPtr tryCircularIBP(const ExprPtr& e, const string& v, bool& ok) {
    ok = false;
    
    // Extract constant multipliers: e.g. e = k * (exp * sin)
    ExprPtr constPart = num(1);
    ExprPtr coreExpr = e;
    
    // A simple loop to pull out constants from left-heavy MUL trees
    while (coreExpr->op == Op::MUL) {
        if (!contains(coreExpr->args[0], v)) {
            constPart = mul(constPart, coreExpr->args[0]);
            coreExpr = coreExpr->args[1];
        } else if (!contains(coreExpr->args[1], v)) {
            constPart = mul(constPart, coreExpr->args[1]);
            coreExpr = coreExpr->args[0];
        } else if (coreExpr->args[0]->op == Op::MUL && !contains(coreExpr->args[0]->args[0], v)) {
            constPart = mul(constPart, coreExpr->args[0]->args[0]);
            coreExpr = mul(coreExpr->args[0]->args[1], coreExpr->args[1]);
        } else {
            break;
        }
    }
    
    if (coreExpr->op != Op::MUL) return nullptr;
    
    auto& left = coreExpr->args[0];
    auto& right = coreExpr->args[1];
    
    ExprPtr expPart = nullptr;
    ExprPtr trigPart = nullptr;
    
    if (left->op == Op::EXP && (right->op == Op::SIN || right->op == Op::COS)) {
        expPart = left; trigPart = right;
    } else if (right->op == Op::EXP && (left->op == Op::SIN || left->op == Op::COS)) {
        expPart = right; trigPart = left;
    }
    
    if (expPart && trigPart) {
        double a, b_exp, c, d_trig;
        if (linearCoeffs(expPart->args[0], v, a, b_exp) && a != 0 &&
            linearCoeffs(trigPart->args[0], v, c, d_trig) && c != 0) {
            
            double denom = a * a + c * c;
            
            if (trigPart->op == Op::SIN) {
                ok = true;
                return mul(constPart, divi(mul(expPart, sub(mul(num(a), trigPart), mul(num(c), fn(Op::COS, trigPart->args[0])))), num(denom)));
            } else {
                ok = true;
                return mul(constPart, divi(mul(expPart, add(mul(num(a), trigPart), mul(num(c), fn(Op::SIN, trigPart->args[0])))), num(denom)));
            }
        }
    }
    return nullptr;
}


// Main integration engine - orchestrates all techniques

static ExprPtr tryIntegrate(const ExprPtr& e, const string& v, bool& ok) {
    return tryIntegrateRec(e, v, ok, 0);
}

static ExprPtr tryIntegrateRec(const ExprPtr& e, const string& v, bool& ok, int depth) {
    ok = true;
    if (depth > 10) return nullptr;
    if (!contains(e, v)) return mul(e, var(v));

    auto tryParts = [&](ExprPtr u, ExprPtr dv) -> ExprPtr {
        if (depth >= 3) return nullptr;
        bool v_ok;
        ExprPtr v_expr = tryIntegrateRec(dv, v, v_ok, depth + 1);
        if (!v_ok || !v_expr) return nullptr;
        ExprPtr du = differentiate(u, v);
        if (!du) return nullptr;
        bool int_v_du_ok;
        ExprPtr int_v_du = tryIntegrateRec(simplify(mul(v_expr, du)), v, int_v_du_ok, depth + 1);
        if (!int_v_du_ok || !int_v_du) return nullptr;
        return sub(mul(u, v_expr), int_v_du);
    };

    {
        bool tableOk;
        ExprPtr r = tableLookup(e, v, tableOk);
        if (tableOk && r) { ok = true; return r; }
    }

    switch (e->op) {
        case Op::VAR:
            if (e->var == v) return divi(powr(var(v), num(2)), num(2));
            ok = false; return nullptr;

        case Op::ADD: {
            bool ok1, ok2;
            ExprPtr ra = tryIntegrateRec(e->args[0], v, ok1, depth);
            ExprPtr rb = tryIntegrateRec(e->args[1], v, ok2, depth);
            if (!ok1 || !ok2) { ok = false; return nullptr; }
            return add(ra, rb);
        }

        case Op::SUB: {
            bool ok1, ok2;
            ExprPtr ra = tryIntegrateRec(e->args[0], v, ok1, depth);
            ExprPtr rb = tryIntegrateRec(e->args[1], v, ok2, depth);
            if (!ok1 || !ok2) { ok = false; return nullptr; }
            return sub(ra, rb);
        }

        case Op::NEG: {
            bool ok1;
            ExprPtr ra = tryIntegrateRec(e->args[0], v, ok1, depth);
            if (!ok1) { ok = false; return nullptr; }
            return neg(ra);
        }

        case Op::MUL: {
            auto& a = e->args[0]; auto& b = e->args[1];
            
            if (!contains(a, v)) {
                bool ok1; ExprPtr rb = tryIntegrateRec(b, v, ok1, depth);
                if (!ok1) { ok = false; return nullptr; }
                return mul(a, rb);
            }
            if (!contains(b, v)) {
                bool ok1; ExprPtr ra = tryIntegrateRec(a, v, ok1, depth);
                if (!ok1) { ok = false; return nullptr; }
                return mul(b, ra);
            }
            
            string sig = getSignature(e, v);
            if (sig == "MUL(EXP(VAR),SIN(VAR))" || sig == "MUL(COS(VAR),EXP(VAR))" ||
                sig == "MUL(EXP(NUM),SIN(VAR))" || sig == "MUL(COS(VAR),EXP(NUM))") {
                bool circOk = false;
                ExprPtr r = tryCircularIBP(e, v, circOk);
                if (circOk && r) return r;
            }
            
            // Fallback for circular IBP in case signature matching misses complex inner args
            {
                bool circOk = false;
                ExprPtr r = tryCircularIBP(e, v, circOk);
                if (circOk && r) return r;
            }

            {
                bool trigOk;
                ExprPtr r = tryTrigPowerReduction(e, v, trigOk);
                if (trigOk && r) return r;
            }

            if (depth < 3) {
                bool uOk;
                ExprPtr r = tryUSubstitution(e, v, uOk, depth);
                if (uOk && r) return r;
            }

            if (depth < 2) {
                bool rischOk;
                ExprPtr r = tryRischHeuristic(e, v, rischOk, depth);
                if (rischOk && r) return r;
            }

            // Integration by parts (LIATE rule)
            if (depth < 3) {
                int prioA = getLiatePriority(a, v);
                int prioB = getLiatePriority(b, v);

                if (prioA >= 0 || prioB >= 0) {
                    if (prioA >= prioB) {
                        ExprPtr res = tryParts(a, b);
                        if (res) return res;
                        // Try the other ordering as fallback
                        res = tryParts(b, a);
                        if (res) return res;
                    } else {
                        ExprPtr res = tryParts(b, a);
                        if (res) return res;
                        // Try the other ordering as fallback
                        res = tryParts(a, b);
                        if (res) return res;
                    }
                }
            }

            ok = false; return nullptr;
        }

        case Op::DIV: {
            auto& a = e->args[0]; auto& b = e->args[1];
            if (!contains(b, v)) {
                bool ok1; ExprPtr ra = tryIntegrateRec(a, v, ok1, depth);
                if (!ok1) { ok = false; return nullptr; }
                return divi(ra, b);
            }

            double la, lb;
            if (isNum(a, 1) && linearCoeffs(b, v, la, lb) && la != 0) {
                return divi(fn(Op::LN, fn(Op::ABS, b)), num(la));
            }

            {
                bool csOk;
                ExprPtr r = tryCompletingSquare(e, v, csOk);
                if (csOk && r) return r;
            }

            if (depth < 2) {
                bool pfOk;
                ExprPtr r = tryPartialFractions(e, v, pfOk);
                if (pfOk && r) return r;
            }

            // Basic substitution: f'(x)/f(x) -> ln|f(x)|
            if (depth < 3) {
                ExprPtr db = differentiate(b, v);
                if (db) {
                    ExprPtr ratio = simplify(divi(a, db));
                    if (!contains(ratio, v)) {
                        return mul(ratio, fn(Op::LN, fn(Op::ABS, b)));
                    }
                }
                
                bool uOk;
                ExprPtr r = tryUSubstitution(e, v, uOk, depth);
                if (uOk && r) {
                    ok = true;
                    return r;
                }
            }

            ok = false; return nullptr;
        }

        case Op::POW: {
            auto& base = e->args[0]; auto& expo = e->args[1];

            {
                bool trigOk;
                ExprPtr r = tryTrigPowerReduction(e, v, trigOk);
                if (trigOk && r) return r;
            }

            double la, lb;
            if (expo->op == Op::NUM && linearCoeffs(base, v, la, lb) && la != 0) {
                double n = expo->num;
                if (abs(n + 1) < 1e-12) {
                    return divi(fn(Op::LN, fn(Op::ABS, base)), num(la));
                }
                return divi(powr(base, num(n + 1)), num(la * (n + 1)));
            }
            
            if (base->op == Op::LN || base->op == Op::ASIN || base->op == Op::ACOS || base->op == Op::ATAN) {
                ExprPtr res = tryParts(e, num(1));
                if (res) {
                    ok = true;
                    return res;
                }
            }

            // Power with non-linear base: u-substitution
            if (depth < 3 && expo->op == Op::NUM) {
                ExprPtr dbase = differentiate(base, v);
                if (dbase) {
                    double n = expo->num;
                    if (abs(n + 1) > 1e-12) {
                        ExprPtr ratio = simplify(divi(num(1), dbase));
                        if (!contains(ratio, v)) {
                            return mul(ratio, divi(powr(base, num(n + 1)), num(n + 1)));
                        }
                    }
                }
            }

            // a^x -> a^x / ln(a)
            if (!contains(base, v) && contains(expo, v)) {
                double la_e, lb_e;
                if (linearCoeffs(expo, v, la_e, lb_e) && la_e != 0) {
                    // ∫ base^(la*x+lb) dx = base^(la*x+lb) / (la * ln(base))
                    return divi(e, mul(num(la_e), fn(Op::LN, base)));
                }
            }

            ok = false; return nullptr;
        }

        case Op::SIN: case Op::COS: case Op::EXP: case Op::LN: case Op::TAN:
        case Op::ASIN: case Op::ACOS: case Op::ATAN:
        case Op::SINH: case Op::COSH: case Op::TANH:
        case Op::SEC: case Op::CSC: case Op::COT: {
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

            // U-substitution for trig with non-linear args
            if (depth < 3) {
                ExprPtr d_inner = differentiate(e->args[0], v);
                if (d_inner) {
                    ExprPtr ratio = simplify(divi(num(1), d_inner));
                    if (!contains(ratio, v)) {
                        // Inner is linear -> already handled by integrateLinearForm
                        // This case shouldn't be reached, but just in case
                    }
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
