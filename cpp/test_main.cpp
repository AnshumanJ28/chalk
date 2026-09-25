#include "ast.hpp"
#include "parser.hpp"
#include "printing.hpp"
#include "differentiation.hpp"
#include "integration.hpp"
#include "vector_calc.hpp"
#include "ode_solver.hpp"
#include "limit.hpp"
#include "solver.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

using namespace calc;

static int passed = 0, failed = 0;

static void check(const string& name, bool cond) {
    if (cond) { cout << "  [PASS] " << name << "\n"; passed++; }
    else      { cout << "  [FAIL] " << name << "\n"; failed++; }
}

int main() {
    cout << "=== Parser Tests ===\n";
    {

        ExprPtr e = parse("x^2 + 3*x - 5");
        check("parse basic polynomial", e != nullptr);

        ExprPtr e1 = parse("asin(x)");
        check("parse asin(x)", e1 != nullptr && e1->op == Op::ASIN);

        ExprPtr e2 = parse("arcsin(x)");
        check("parse arcsin(x)", e2 != nullptr && e2->op == Op::ASIN);

        ExprPtr e3 = parse("arccos(x)");
        check("parse arccos(x)", e3 != nullptr && e3->op == Op::ACOS);

        ExprPtr e4 = parse("arctan(x)");
        check("parse arctan(x)", e4 != nullptr && e4->op == Op::ATAN);

        ExprPtr e5 = parse("sin^-1(x)");
        check("parse sin^-1(x)", e5 != nullptr && e5->op == Op::ASIN);

        ExprPtr e6 = parse("cos^(-1)(x)");
        check("parse cos^(-1)(x)", e6 != nullptr && e6->op == Op::ACOS);

        ExprPtr e7 = parse("tan^-1(x)");
        check("parse tan^-1(x)", e7 != nullptr && e7->op == Op::ATAN);
    }

    cout << "\n=== Printing Tests ===\n";
    {
        ExprPtr e = parse("x^2 + 3*x - 5");
        string s = toString(simplify(e));
        check("toString polynomial", !s.empty());
        cout << "    toString: " << s << "\n";

        string l = toLatex(simplify(e));
        check("toLatex polynomial", !l.empty());
        cout << "    toLatex: " << l << "\n";

        ExprPtr frac = parse("(x+1)/(x-2)");
        string fl = toLatex(frac);
        check("toLatex fraction has \\frac", fl.find("\\frac") != string::npos);
        cout << "    toLatex fraction: " << fl << "\n";

        ExprPtr invTrig = parse("asin(x)");
        string itl = toLatex(invTrig);
        check("toLatex asin -> sin^{-1}", itl.find("\\sin^{-1}") != string::npos);
        cout << "    toLatex invTrig: " << itl << "\n";
    }

    cout << "\n=== Differentiation Tests ===\n";
    {

        ExprPtr e = parse("x^3");
        ExprPtr d = simplify(differentiate(e, "x"));
        map<string, double> m{{"x", 2.0}};
        double val = evaluate(d, m);
        check("d/dx x^3 at x=2 -> 12", abs(val - 12.0) < 1e-6);

        ExprPtr e2 = parse("sin(x)");
        ExprPtr d2 = simplify(differentiate(e2, "x"));
        map<string, double> m2{{"x", 0.0}};
        double val2 = evaluate(d2, m2);
        check("d/dx sin(x) at x=0 -> 1", abs(val2 - 1.0) < 1e-6);

        ExprPtr e3 = parse("x^4");
        ExprPtr d3 = differentiateOrder(e3, "x", 2);
        map<string, double> m3{{"x", 1.0}};
        double val3 = evaluate(d3, m3);
        check("d²/dx² x^4 at x=1 -> 12", abs(val3 - 12.0) < 1e-6);

        ExprPtr e4 = parse("sin(x)");
        ExprPtr d4 = differentiateOrder(e4, "x", 3);
        map<string, double> m4{{"x", 0.0}};
        double val4 = evaluate(d4, m4);
        check("d³/dx³ sin(x) at x=0 -> -1", abs(val4 + 1.0) < 1e-3);

        ExprPtr e5 = parse("x^2*y^3");
        ExprPtr d5 = mixedPartial(e5, {"x", "y"});
        map<string, double> m5{{"x", 1.0}, {"y", 2.0}};
        double val5 = evaluate(d5, m5);
        check("∂²/∂x∂y x²y³ at (1,2) -> 24", abs(val5 - 24.0) < 1e-6);
    }

    cout << "\n=== Integration Tests ===\n";
    {

        CalcResult r1 = integrate(parse("x^2"), "x");
        check("∫ x² dx found", r1.ok);
        cout << "    Result: " << r1.text << "\n";

        DefiniteResult dr = definiteIntegral(parse("x^2"), "x", 0, 1);
        check("∫₀¹ x² dx ≈ 1/3", dr.ok && abs(dr.value - 1.0/3.0) < 1e-6);
        cout << "    ∫₀¹ x² dx = " << dr.value << "\n";

        DefiniteResult dr2 = definiteIntegral(parse("sin(x)"), "x", 0, M_PI);
        check("∫₀^π sin(x) dx ≈ 2", dr2.ok && abs(dr2.value - 2.0) < 1e-6);
        cout << "    ∫₀^π sin(x) dx = " << dr2.value << "\n";
    }

    cout << "\n=== Vector Calculus Tests ===\n";
    {

        VectorCalcResult gr = gradient("x^2 + y^2 + z^2", {"x", "y", "z"});
        check("gradient computed", gr.ok && gr.components.size() == 3);
        if (gr.ok) cout << "    ∇f = [" << gr.components[0] << ", " << gr.components[1] << ", " << gr.components[2] << "]\n";

        VectorCalcResult dv = divergence({"x^2", "y^2", "z^2"}, {"x", "y", "z"});
        check("divergence computed", dv.ok);
        if (dv.ok) cout << "    ∇·F = " << dv.scalar << "\n";

        VectorCalcResult cr = curl({"-y", "x", "0"}, {"x", "y", "z"});
        check("curl computed", cr.ok && cr.components.size() == 3);
        if (cr.ok) cout << "    ∇×F = [" << cr.components[0] << ", " << cr.components[1] << ", " << cr.components[2] << "]\n";

        VectorCalcResult lp = laplacian("x^3 + y^3 + z^3", {"x", "y", "z"});
        check("laplacian computed", lp.ok);
        if (lp.ok) cout << "    ∇²f = " << lp.scalar << "\n";
    }

    cout << "\n=== Limit Tests ===\n";
    {

        LimitResult lr = limit(parse("sin(x)/x"), "x", 0, false, false);
        check("lim sin(x)/x -> 0 = 1", lr.ok && lr.text == "1");
    }

    cout << "\n=== ODE Tests ===\n";
    {

        ODEResult or1 = solveODE("x");
        check("ODE dy/dx=x symbolic", or1.ok);
        if (or1.ok) cout << "    Solution: " << or1.symbolic_solution << "\n";

        ODETrajectory tr = computeTrajectory("x + y", 0, 1, -3, 3, 200);
        check("RK4 trajectory computed", tr.ok && tr.points.size() > 10);

        SlopeFieldData sf = computeSlopeField("x + y", -3, 3, -3, 3, 10, 10);
        check("slope field computed", sf.ok && sf.points.size() > 0);
    }

    cout << "\n=== Solver Tests ===\n";
    {

        SolveResult sr = solve(parse("x^2 - 5*x + 6"), "x");
        check("quadratic x²-5x+6 solved", sr.ok && sr.roots.size() == 2);
    }

    cout << "\n==================\n";
    cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    return failed > 0 ? 1 : 0;
}
