#include "series.hpp"
#include "differentiation.hpp"
#include "printing.hpp"
#include <cmath>
#include <map>
#include <sstream>
using namespace std;

namespace calc {

static double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

SeriesResult taylorSeries(const ExprPtr& expr, const string& var,
                          double a, int terms) {
    SeriesResult res;
    res.ok = false;

    if (terms <= 0 || terms > 30) {
        res.message = "Number of terms must be between 1 and 30.";
        return res;
    }

    try {
        ExprPtr se = simplify(expr);

        ExprPtr poly = num(0);
        ExprPtr currentDeriv = se;

        for (int k = 0; k < terms; k++) {

            ExprPtr simpDeriv = simplify(currentDeriv);
            map<string, double> evalPoint;
            evalPoint[var] = a;

            double coeff;
            try {
                coeff = evaluate(simpDeriv, evalPoint);
            } catch (...) {
                res.message = "Could not evaluate derivative at the center point a = " +
                              numToStr(a) + " (term " + to_string(k) + ").";
                return res;
            }

            if (!isfinite(coeff)) {
                res.message = "Derivative is infinite or undefined at a = " +
                              numToStr(a) + " (term " + to_string(k) + ").";
                return res;
            }

            double termCoeff = coeff / factorial(k);

            if (abs(termCoeff) > 1e-15) {
                ExprPtr termExpr;
                if (k == 0) {
                    termExpr = num(termCoeff);
                } else if (abs(a) < 1e-12) {

                    if (k == 1) {
                        termExpr = mul(num(termCoeff), calc::var(var));
                    } else {
                        termExpr = mul(num(termCoeff), powr(calc::var(var), num(k)));
                    }
                } else {

                    ExprPtr shift = sub(calc::var(var), num(a));
                    if (k == 1) {
                        termExpr = mul(num(termCoeff), shift);
                    } else {
                        termExpr = mul(num(termCoeff), powr(shift, num(k)));
                    }
                }
                poly = add(poly, termExpr);
            }

            if (k < terms - 1) {
                currentDeriv = differentiate(currentDeriv, var);
            }
        }

        poly = simplify(poly);
        res.ok = true;
        res.text = toString(poly);
        res.latex = toLatex(poly);
    } catch (const exception& ex) {
        res.message = string("Error computing Taylor series: ") + ex.what();
    }

    return res;
}

}
