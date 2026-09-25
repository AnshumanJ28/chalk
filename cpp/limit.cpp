#include "limit.hpp"
#include "differentiation.hpp"
#include "printing.hpp"
#include <cmath>
#include <functional>
#include <map>
using namespace std;

namespace calc {

static bool finiteVal(double x) { return isfinite(x); }

LimitResult limit(const ExprPtr& e, const string& v,
                  double point, bool posInf, bool negInf) {
    ExprPtr se = simplify(e);

    auto tryEval = [&](double x) -> double {
        map<string, double> m{{v, x}};
        try { return evaluate(se, m); } catch (...) { return NAN; }
    };

    if (posInf || negInf) {
        double sign = posInf ? 1.0 : -1.0;
        double v1 = tryEval(sign * 1e4);
        double v2 = tryEval(sign * 1e6);
        double v3 = tryEval(sign * 1e8);

        if (finiteVal(v1) && finiteVal(v2) && finiteVal(v3) &&
            abs(v2 - v3) < 1e-4 * max(1.0, abs(v3))) {
            string txt = numToStr(v3);
            return {true, txt, txt};
        }
        if (finiteVal(v3) && abs(v3) > 1e8) {
            string txt = v3 > 0 ? "\\infty" : "-\\infty";
            string plain = v3 > 0 ? "infinity" : "-infinity";
            return {true, plain, txt};
        }
        if (finiteVal(v1) && finiteVal(v3) &&
            abs(v3) > abs(v1) * 10) {
            string txt = v3 > 0 ? "\\infty" : "-\\infty";
            string plain = v3 > 0 ? "infinity" : "-infinity";
            return {true, plain, txt};
        }
        return {false,
            "Could not determine the limit numerically as it does not appear to converge.",
            ""};
    }

    double direct = tryEval(point);
    if (finiteVal(direct)) {
        string txt = numToStr(direct);
        return {true, txt, txt};
    }

    function<LimitResult(const ExprPtr&, int)> lhopital =
        [&](const ExprPtr& expr, int depth) -> LimitResult {
        if (depth > 6) return {false, "", ""};
        if (expr->op != Op::DIV) return {false, "", ""};

        map<string, double> m{{v, point}};
        double n0, d0;
        try { n0 = evaluate(expr->args[0], m); } catch (...) { n0 = NAN; }
        try { d0 = evaluate(expr->args[1], m); } catch (...) { d0 = NAN; }

        bool zeroZero = (isnan(n0) || abs(n0) < 1e-9) &&
                        (isnan(d0) || abs(d0) < 1e-9);
        bool infInf = (!finiteVal(n0) || abs(n0) > 1e8) &&
                      (!finiteVal(d0) || abs(d0) > 1e8);

        if (!zeroZero && !infInf) return {false, "", ""};

        ExprPtr dn = simplify(differentiate(expr->args[0], v));
        ExprPtr dd = simplify(differentiate(expr->args[1], v));
        ExprPtr ratio = divi(dn, dd);

        map<string, double> mm{{v, point}};
        double val;
        try { val = evaluate(simplify(ratio), mm); } catch (...) { val = NAN; }
        if (finiteVal(val)) {
            string txt = numToStr(val);
            return {true, txt, txt};
        }

        return lhopital(simplify(ratio), depth + 1);
    };

    LimitResult lh = lhopital(se, 0);
    if (lh.ok) return lh;

    double hl  = tryEval(point - 1e-5);
    double hl2 = tryEval(point - 1e-7);
    double hr  = tryEval(point + 1e-5);
    double hr2 = tryEval(point + 1e-7);

    bool leftOk  = finiteVal(hl) && finiteVal(hl2) && abs(hl - hl2) < 1e-3;
    bool rightOk = finiteVal(hr) && finiteVal(hr2) && abs(hr - hr2) < 1e-3;

    if (leftOk && rightOk && abs(hl2 - hr2) < 1e-3) {
        string txt = numToStr((hl2 + hr2) / 2.0);
        return {true, txt, txt};
    }
    if (leftOk && rightOk) {
        return {false,
            "The left-hand limit (" + numToStr(hl2) + ") and right-hand limit (" +
            numToStr(hr2) + ") differ, so the limit does not exist.",
            ""};
    }
    if ((finiteVal(hl2) && abs(hl2) > 1e6) ||
        (finiteVal(hr2) && abs(hr2) > 1e6)) {
        return {true,
            "infinity (or -infinity) -- function diverges near this point",
            "\\pm\\infty"};
    }

    return {false,
        "Could not determine this limit with the available methods.",
        ""};
}

}
