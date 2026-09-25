#pragma once
#include "ode_solver.hpp"
#include "parser.hpp"
#include "printing.hpp"
#include "differentiation.hpp"
#include "integration.hpp"
#include <cmath>
#include <map>
#include <sstream>
using namespace std;

namespace calc {

ODEResult solveODE(const string& dydx_expr) {
    ODEResult res;
    res.ok = false;

    try {
        ExprPtr expr = parse(dydx_expr);
        ExprPtr se = simplify(expr);

        if (!contains(se, "y")) {
            bool intOk;
            ExprPtr antideriv = simplify(se);
            CalcResult intRes = integrate(antideriv, "x");
            if (intRes.ok) {
                res.ok = true;
                res.symbolic_solution = "y = " + intRes.text;
                res.latex_solution = "y = " + intRes.latex;
                res.method = "direct integration";
                return res;
            }
        }

        if (se->op == Op::MUL && se->args.size() == 2) {
            auto& a = se->args[0];
            auto& b = se->args[1];
            bool aHasX = contains(a, "x"), aHasY = contains(a, "y");
            bool bHasX = contains(b, "x"), bHasY = contains(b, "y");

            ExprPtr gx = nullptr, hy = nullptr;
            if (aHasX && !aHasY && bHasY && !bHasX) { gx = a; hy = b; }
            else if (bHasX && !bHasY && aHasY && !aHasX) { gx = b; hy = a; }

            if (gx && hy) {
                ExprPtr invH = simplify(divi(num(1), hy));
                CalcResult lhsInt = integrate(invH, "y");
                CalcResult rhsInt = integrate(gx, "x");

                if (lhsInt.ok && rhsInt.ok) {
                    res.ok = true;
                    res.symbolic_solution = lhsInt.text + " = " + rhsInt.text;
                    res.latex_solution = lhsInt.latex + " = " + rhsInt.latex;
                    res.method = "separable";
                    return res;
                }
            }
        }

        res.message = "Could not find a symbolic closed-form solution. "
                      "Use the numerical solver (RK4) with an initial condition for a trajectory plot.";
    } catch (const exception& ex) {
        res.message = string("Error: ") + ex.what();
    }

    return res;
}

ODEResult solveODE2(double a, double b, double c) {
    ODEResult res;
    res.ok = false;

    if (abs(a) < 1e-12) {
        res.message = "Not a 2nd order ODE (a = 0).";
        return res;
    }

    double discriminant = b * b - 4 * a * c;

    ostringstream text, latex;

    if (discriminant > 1e-12) {

        double r1 = (-b + sqrt(discriminant)) / (2 * a);
        double r2 = (-b - sqrt(discriminant)) / (2 * a);
        text << "y(x) = C1 * e^(" << numToStr(r1) << "*x) + C2 * e^(" << numToStr(r2) << "*x)";
        latex << "y(x) = C_1 e^{" << numToStr(r1) << "x} + C_2 e^{" << numToStr(r2) << "x}";
        res.method = "2nd order linear homogeneous (distinct real roots)";
    } else if (abs(discriminant) <= 1e-12) {

        double r = -b / (2 * a);
        text << "y(x) = C1 * e^(" << numToStr(r) << "*x) + C2 * x * e^(" << numToStr(r) << "*x)";
        latex << "y(x) = C_1 e^{" << numToStr(r) << "x} + C_2 x e^{" << numToStr(r) << "x}";
        res.method = "2nd order linear homogeneous (repeated real root)";
    } else {

        double alpha = -b / (2 * a);
        double beta = sqrt(-discriminant) / (2 * a);
        if (abs(alpha) < 1e-12) {
            text << "y(x) = C1 * cos(" << numToStr(beta) << "*x) + C2 * sin(" << numToStr(beta) << "*x)";
            latex << "y(x) = C_1 \\cos(" << numToStr(beta) << "x) + C_2 \\sin(" << numToStr(beta) << "x)";
        } else {
            text << "y(x) = e^(" << numToStr(alpha) << "*x) * (C1 * cos(" << numToStr(beta) << "*x) + C2 * sin(" << numToStr(beta) << "*x))";
            latex << "y(x) = e^{" << numToStr(alpha) << "x} \\left( C_1 \\cos(" << numToStr(beta) << "x) + C_2 \\sin(" << numToStr(beta) << "x) \\right)";
        }
        res.method = "2nd order linear homogeneous (complex roots)";
    }

    res.symbolic_solution = text.str();
    res.latex_solution = latex.str();
    res.ok = true;
    return res;
}

ODETrajectory computeTrajectory(const string& dydx_expr,
                                double x0, double y0,
                                double xmin, double xmax,
                                int steps) {
    ODETrajectory res;
    res.ok = false;

    try {
        ExprPtr expr = parse(dydx_expr);

        auto f = [&](double x, double y) -> double {
            map<string, double> vals;
            vals["x"] = x;
            vals["y"] = y;
            return evaluate(expr, vals);
        };

        vector<pair<double, double>> forwardPts;
        {
            double h = (xmax - x0) / max(steps / 2, 1);
            double x = x0, y = y0;
            forwardPts.push_back({x, y});

            for (int i = 0; i < steps / 2 && x < xmax; i++) {
                try {
                    double k1 = h * f(x, y);
                    double k2 = h * f(x + h / 2, y + k1 / 2);
                    double k3 = h * f(x + h / 2, y + k2 / 2);
                    double k4 = h * f(x + h, y + k3);
                    y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
                    x += h;
                    if (!isfinite(y) || abs(y) > 1e12) break;
                    forwardPts.push_back({x, y});
                } catch (...) { break; }
            }
        }

        vector<pair<double, double>> backwardPts;
        {
            double h = (x0 - xmin) / max(steps / 2, 1);
            double x = x0, y = y0;

            for (int i = 0; i < steps / 2 && x > xmin; i++) {
                try {
                    double nh = -h;
                    double k1 = nh * f(x, y);
                    double k2 = nh * f(x + nh / 2, y + k1 / 2);
                    double k3 = nh * f(x + nh / 2, y + k2 / 2);
                    double k4 = nh * f(x + nh, y + k3);
                    y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
                    x -= h;
                    if (!isfinite(y) || abs(y) > 1e12) break;
                    backwardPts.push_back({x, y});
                } catch (...) { break; }
            }
        }

        res.points.reserve(backwardPts.size() + forwardPts.size());
        for (auto it = backwardPts.rbegin(); it != backwardPts.rend(); ++it)
            res.points.push_back(*it);
        for (auto& pt : forwardPts)
            res.points.push_back(pt);

        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing trajectory: ") + ex.what();
    }

    return res;
}

ODETrajectory computeTrajectory2(double a, double b, double c,
                                 const string& d2ydx2_expr,
                                 double x0, double y0, double v0,
                                 double xmin, double xmax,
                                 int steps) {
    ODETrajectory res;
    res.ok = false;

    try {
        ExprPtr expr = parse(d2ydx2_expr);

        auto fv = [&](double x, double y, double v) -> double {
            map<string, double> vals;
            vals["x"] = x;
            vals["y"] = y;
            vals["v"] = v;
            double fx = evaluate(expr, vals);
            return (fx - b * v - c * y) / a;
        };

        vector<pair<double, double>> forwardPts;
        {
            double h = (xmax - x0) / max(steps / 2, 1);
            double x = x0, y = y0, v = v0;
            forwardPts.push_back({x, y});

            for (int i = 0; i < steps / 2 && x < xmax; i++) {
                try {
                    double ky1 = h * v;
                    double kv1 = h * fv(x, y, v);

                    double ky2 = h * (v + kv1 / 2);
                    double kv2 = h * fv(x + h / 2, y + ky1 / 2, v + kv1 / 2);

                    double ky3 = h * (v + kv2 / 2);
                    double kv3 = h * fv(x + h / 2, y + ky2 / 2, v + kv2 / 2);

                    double ky4 = h * (v + kv3);
                    double kv4 = h * fv(x + h, y + ky3, v + kv3);

                    y += (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6.0;
                    v += (kv1 + 2 * kv2 + 2 * kv3 + kv4) / 6.0;
                    x += h;

                    if (!isfinite(y) || abs(y) > 1e12 || !isfinite(v)) break;
                    forwardPts.push_back({x, y});
                } catch (...) { break; }
            }
        }

        vector<pair<double, double>> backwardPts;
        {
            double h = (x0 - xmin) / max(steps / 2, 1);
            double x = x0, y = y0, v = v0;

            for (int i = 0; i < steps / 2 && x > xmin; i++) {
                try {
                    double nh = -h;
                    double ky1 = nh * v;
                    double kv1 = nh * fv(x, y, v);

                    double ky2 = nh * (v + kv1 / 2);
                    double kv2 = nh * fv(x + nh / 2, y + ky1 / 2, v + kv1 / 2);

                    double ky3 = nh * (v + kv2 / 2);
                    double kv3 = nh * fv(x + nh / 2, y + ky2 / 2, v + kv2 / 2);

                    double ky4 = nh * (v + kv3);
                    double kv4 = nh * fv(x + nh, y + ky3, v + kv3);

                    y += (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6.0;
                    v += (kv1 + 2 * kv2 + 2 * kv3 + kv4) / 6.0;
                    x -= h;

                    if (!isfinite(y) || abs(y) > 1e12 || !isfinite(v)) break;
                    backwardPts.push_back({x, y});
                } catch (...) { break; }
            }
        }

        res.points.reserve(backwardPts.size() + forwardPts.size());
        for (auto it = backwardPts.rbegin(); it != backwardPts.rend(); ++it)
            res.points.push_back(*it);
        for (auto& pt : forwardPts)
            res.points.push_back(pt);

        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing 2nd order trajectory: ") + ex.what();
    }

    return res;
}

SlopeFieldData computeSlopeField(const string& dydx_expr,
                                 double xmin, double xmax,
                                 double ymin, double ymax,
                                 int grid_nx, int grid_ny) {
    SlopeFieldData res;
    res.ok = false;

    try {
        ExprPtr expr = parse(dydx_expr);

        double dx = (xmax - xmin) / max(grid_nx - 1, 1);
        double dy = (ymax - ymin) / max(grid_ny - 1, 1);

        res.points.reserve(grid_nx * grid_ny);

        for (int i = 0; i < grid_nx; i++) {
            for (int j = 0; j < grid_ny; j++) {
                double x = xmin + i * dx;
                double y = ymin + j * dy;

                SlopeFieldPoint pt;
                pt.x = x;
                pt.y = y;
                pt.valid = true;

                try {
                    map<string, double> vals;
                    vals["x"] = x;
                    vals["y"] = y;
                    pt.slope = evaluate(expr, vals);
                    if (!isfinite(pt.slope)) pt.valid = false;
                } catch (...) {
                    pt.slope = 0;
                    pt.valid = false;
                }

                res.points.push_back(pt);
            }
        }

        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing slope field: ") + ex.what();
    }

    return res;
}

}
