#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ast.hpp"
#include "parser.hpp"
#include "printing.hpp"
#include "differentiation.hpp"
#include "integration.hpp"
#include "vector_calc.hpp"
#include "ode_solver.hpp"
#include "limit.hpp"
#include "solver.hpp"
#include "series.hpp"
#include "laplace.hpp"
#include "complex.hpp"
#include "fourier.hpp"
using namespace std;

namespace py = pybind11;
using namespace calc;

static py::dict py_differentiate(const string& expr, const string& v, int order,
                                 const vector<string>& vars) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        ExprPtr result;

        if (!vars.empty()) {

            result = mixedPartial(e, vars);
        } else {

            result = differentiateOrder(e, v, order);
        }

        result = simplify(result);
        out["ok"] = true;
        out["result"] = toString(result);
        out["latex"] = toLatex(result);
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_integrate(const string& expr, const string& v) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        CalcResult r = integrate(e, v);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_definite_integral(const string& expr, const string& v,
                                     double lower, double upper) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        DefiniteResult r = definiteIntegral(e, v, lower, upper);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
        out["value"] = r.value;
        out["is_numeric"] = r.is_numeric;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_multi_integral(const string& expr, py::list specs_list) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        vector<IntegralSpec> specs;
        for (auto item : specs_list) {
            py::dict d = item.cast<py::dict>();
            IntegralSpec s;
            s.var = d["var"].cast<string>();
            s.lower = d.contains("lower") ? d["lower"].cast<double>() : NAN;
            s.upper = d.contains("upper") ? d["upper"].cast<double>() : NAN;
            specs.push_back(s);
        }
        MultiIntegralResult r = multiIntegral(e, specs);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
        out["is_numeric"] = r.is_numeric;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_vector_calc(const string& op,
                               const string& scalar_expr,
                               const vector<string>& vector_exprs,
                               const vector<string>& vars) {
    py::dict out;
    try {
        if (op == "gradient") {
            VectorCalcResult r = gradient(scalar_expr, vars);
            out["ok"] = r.ok;
            out["components"] = r.components;
            out["latex_components"] = r.latex_components;
            out["message"] = r.message;
        } else if (op == "divergence") {
            VectorCalcResult r = divergence(vector_exprs, vars);
            out["ok"] = r.ok;
            out["scalar"] = r.scalar;
            out["latex_scalar"] = r.latex_scalar;
            out["message"] = r.message;
        } else if (op == "curl") {
            VectorCalcResult r = curl(vector_exprs, vars);
            out["ok"] = r.ok;
            out["components"] = r.components;
            out["latex_components"] = r.latex_components;
            out["message"] = r.message;
        } else if (op == "laplacian") {
            VectorCalcResult r = laplacian(scalar_expr, vars);
            out["ok"] = r.ok;
            out["scalar"] = r.scalar;
            out["latex_scalar"] = r.latex_scalar;
            out["message"] = r.message;
        } else {
            out["ok"] = false;
            out["message"] = "Unknown operation: " + op;
        }
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_solve_ode(const string& dydx_expr,
                             double x0, double y0,
                             double xmin, double xmax, int steps) {
    py::dict out;
    try {

        ODEResult sym = solveODE(dydx_expr);
        out["symbolic_ok"] = sym.ok;
        out["symbolic_solution"] = sym.symbolic_solution;
        out["latex_solution"] = sym.latex_solution;
        out["method"] = sym.method;
        out["symbolic_message"] = sym.message;

        ODETrajectory traj = computeTrajectory(dydx_expr, x0, y0, xmin, xmax, steps);
        out["trajectory_ok"] = traj.ok;
        py::list points;
        for (auto& pt : traj.points) {
            py::list p;
            p.append(pt.first);
            p.append(pt.second);
            points.append(p);
        }
        out["trajectory"] = points;
        out["trajectory_message"] = traj.message;

        SlopeFieldData sf = computeSlopeField(dydx_expr, xmin, xmax,
                                              y0 - (xmax - xmin) / 2,
                                              y0 + (xmax - xmin) / 2, 20, 20);
        out["slope_field_ok"] = sf.ok;
        py::list sfPoints;
        for (auto& pt : sf.points) {
            if (pt.valid) {
                py::dict d;
                d["x"] = pt.x;
                d["y"] = pt.y;
                d["slope"] = pt.slope;
                sfPoints.append(d);
            }
        }
        out["slope_field"] = sfPoints;

        out["ok"] = sym.ok || traj.ok;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_solve_ode2(double a, double b, double c,
                              const string& d2ydx2_expr,
                              double x0, double y0, double v0,
                              double xmin, double xmax, int steps) {
    py::dict out;
    try {

        ODEResult sym = solveODE2(a, b, c);
        out["symbolic_ok"] = sym.ok;
        out["symbolic_solution"] = sym.symbolic_solution;
        out["latex_solution"] = sym.latex_solution;
        out["method"] = sym.method;
        out["symbolic_message"] = sym.message;

        ODETrajectory traj = computeTrajectory2(a, b, c, d2ydx2_expr, x0, y0, v0, xmin, xmax, steps);
        out["trajectory_ok"] = traj.ok;
        py::list points;
        for (auto& pt : traj.points) {
            py::list p;
            p.append(pt.first);
            p.append(pt.second);
            points.append(p);
        }
        out["trajectory"] = points;
        out["trajectory_message"] = traj.message;

        out["ok"] = sym.ok || traj.ok;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_limit(const string& expr, const string& v,
                         const string& pointStr) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        bool posInf = (pointStr == "inf" || pointStr == "+inf" || pointStr == "infinity");
        bool negInf = (pointStr == "-inf" || pointStr == "-infinity");
        double point = 0.0;
        if (!posInf && !negInf) point = stod(pointStr);
        LimitResult r = limit(e, v, point, posInf, negInf);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_solve(const string& expr, const string& v) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        SolveResult r = solve(e, v);
        out["ok"] = r.ok;
        out["roots"] = r.roots;
        out["latex_roots"] = r.latex_roots;
        out["message"] = r.message;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["roots"] = vector<string>{};
        out["latex_roots"] = vector<string>{};
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_solve_system(const vector<string>& equations,
                                const vector<string>& variables) {
    py::dict out;
    try {
        SystemSolveResult r = solveSystem(equations, variables);
        out["ok"] = r.ok;
        if (r.ok) {
            py::dict solution;
            for (auto const& [varName, value] : r.solution) {
                solution[py::str(varName)] = value;
            }
            out["solution"] = solution;
        } else {
            out["message"] = r.message;
        }
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_evaluate(const string& expr, py::dict values) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        map<string, double> vals;
        for (auto item : values) vals[py::str(item.first)] = item.second.cast<double>();
        double r = evaluate(e, vals);
        out["ok"] = true;
        out["result"] = r;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = 0.0;
        out["error"] = string(ex.what());
    }
    return out;
}

static py::dict py_batch_evaluate(const vector<string>& exprs,
                                  const string& var,
                                  double xmin, double xmax, int steps) {
    py::dict out;
    try {
        double dx = (xmax - xmin) / max(steps - 1, 1);
        py::list series;

        for (const auto& exprStr : exprs) {
            ExprPtr e = parse(exprStr);
            py::dict s;
            s["name"] = exprStr;
            py::list points;

            for (int i = 0; i < steps; i++) {
                double x = xmin + i * dx;
                map<string, double> vals;
                vals[var] = x;
                try {
                    double y = evaluate(e, vals);
                    if (isfinite(y)) {
                        py::list pt;
                        pt.append(x);
                        pt.append(y);
                        points.append(pt);
                    }
                } catch (...) {

                }
            }
            s["points"] = points;
            series.append(s);
        }

        out["ok"] = true;
        out["series"] = series;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_vector_field_grid(const string& P_expr, const string& Q_expr,
                                     double xmin, double xmax,
                                     double ymin, double ymax, int grid_n) {
    py::dict out;
    try {
        ExprPtr Pe = parse(P_expr);
        ExprPtr Qe = parse(Q_expr);

        double dx = (xmax - xmin) / max(grid_n - 1, 1);
        double dy = (ymax - ymin) / max(grid_n - 1, 1);

        py::list arrows;
        for (int i = 0; i < grid_n; i++) {
            for (int j = 0; j < grid_n; j++) {
                double x = xmin + i * dx;
                double y = ymin + j * dy;
                map<string, double> vals;
                vals["x"] = x;
                vals["y"] = y;
                try {
                    double px = evaluate(Pe, vals);
                    double qx = evaluate(Qe, vals);
                    if (isfinite(px) && isfinite(qx)) {
                        double mag = sqrt(px * px + qx * qx);
                        py::dict arrow;
                        arrow["x"] = x;
                        arrow["y"] = y;
                        arrow["dx"] = px;
                        arrow["dy"] = qx;
                        arrow["mag"] = mag;
                        arrows.append(arrow);
                    }
                } catch (...) {}
            }
        }

        out["ok"] = true;
        out["arrows"] = arrows;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_taylor_series(const string& expr, const string& var,
                                 double a, int terms) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        SeriesResult r = taylorSeries(e, var, a, terms);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
        out["message"] = r.message;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_laplace_transform(const string& expr, const string& t_var,
                                     const string& s_var) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        LaplaceResult r = laplaceTransform(e, t_var, s_var);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
        out["method"] = r.method;
        out["message"] = r.message;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_inverse_laplace_transform(const string& expr, const string& s_var,
                                             const string& t_var) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        LaplaceResult r = inverseLaplaceTransform(e, s_var, t_var);
        out["ok"] = r.ok;
        out["result"] = r.text;
        out["latex"] = r.latex;
        out["method"] = r.method;
        out["message"] = r.message;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_compute_residue(const string& expr, const string& var, double real_part, double imag_part) {
    py::dict out;
    try {
        complex<double> z0(real_part, imag_part);
        complex<double> res = compute_residue(expr, var, z0);
        out["ok"] = true;
        out["real"] = res.real();
        out["imag"] = res.imag();
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_contour_integral_circle(const string& expr, const string& var, double cx, double cy, double radius, int steps) {
    py::dict out;
    try {
        complex<double> center(cx, cy);
        complex<double> res = contour_integral_circle(expr, var, center, radius, steps);
        out["ok"] = true;
        out["real"] = res.real();
        out["imag"] = res.imag();
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_fourier_transform(const string& expr, const string& t_var, const string& w_var) {
    py::dict out;
    try {
        string res = fourier_transform(expr, t_var, w_var);
        out["ok"] = true;
        out["result"] = res;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_line_integral(const vector<string>& F_exprs, const vector<string>& path_exprs, const string& t_var, double t_start, double t_end, int steps) {
    py::dict out;
    try {
        double res = line_integral(F_exprs, path_exprs, t_var, t_start, t_end, steps);
        out["ok"] = true;
        out["value"] = res;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_surface_integral(const vector<string>& F_exprs, const vector<string>& surface_exprs, const vector<string>& uv_vars, double u_start, double u_end, double v_start, double v_end, int steps) {
    py::dict out;
    try {
        double res = surface_integral(F_exprs, surface_exprs, uv_vars, u_start, u_end, v_start, v_end, steps);
        out["ok"] = true;
        out["value"] = res;
    } catch (const exception& ex) {
        out["ok"] = false;
        out["message"] = string("Error: ") + ex.what();
    }
    return out;
}

static py::dict py_simplify(const string& expr) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        ExprPtr s = simplify(e);
        out["ok"] = true;
        out["result"] = toString(s);
        out["latex"] = toLatex(s);
    } catch (const exception& ex) {
        out["ok"] = false;
        out["result"] = string("Error: ") + ex.what();
        out["latex"] = "";
    }
    return out;
}

static py::dict py_to_latex(const string& expr) {
    py::dict out;
    try {
        ExprPtr e = parse(expr);
        out["ok"] = true;
        out["latex"] = toLatex(e);
    } catch (const exception& ex) {
        out["ok"] = false;
        out["latex"] = "";
    }
    return out;
}

PYBIND11_MODULE(calcengine, m) {
    m.doc() = "C++ symbolic calculator engine (pybind11 bindings) — modular architecture";

    m.def("differentiate", &py_differentiate,
          "Differentiate an expression (supports nth-order and mixed partials)",
          py::arg("expr"), py::arg("var") = "x",
          py::arg("order") = 1, py::arg("vars") = vector<string>{});

    m.def("integrate", &py_integrate,
          "Compute indefinite integral",
          py::arg("expr"), py::arg("var") = "x");

    m.def("definite_integral", &py_definite_integral,
          "Compute definite integral with bounds",
          py::arg("expr"), py::arg("var") = "x",
          py::arg("lower") = 0.0, py::arg("upper") = 1.0);

    m.def("multi_integral", &py_multi_integral,
          "Compute double/triple iterated integral",
          py::arg("expr"), py::arg("specs"));

    m.def("vector_calc", &py_vector_calc,
          "Vector calculus operations (gradient, divergence, curl, laplacian)",
          py::arg("op"), py::arg("scalar_expr") = "",
          py::arg("vector_exprs") = vector<string>{},
          py::arg("vars") = vector<string>{"x", "y", "z"});

    m.def("solve_ode", &py_solve_ode,
          "Solve ODE dy/dx = f(x,y) symbolically + numerically (RK4 + slope field)",
          py::arg("dydx_expr"), py::arg("x0") = 0.0, py::arg("y0") = 1.0,
          py::arg("xmin") = -5.0, py::arg("xmax") = 5.0, py::arg("steps") = 500);

    m.def("solve_ode2", &py_solve_ode2,
          "Solve 2nd-order ODE ay''+by'+cy=0 symbolically + numerically (RK4)",
          py::arg("a"), py::arg("b"), py::arg("c"),
          py::arg("d2ydx2_expr"), py::arg("x0") = 0.0, py::arg("y0") = 1.0, py::arg("v0") = 0.0,
          py::arg("xmin") = -5.0, py::arg("xmax") = 5.0, py::arg("steps") = 500);

    m.def("limit", &py_limit,
          "Compute a limit",
          py::arg("expr"), py::arg("var") = "x", py::arg("point") = "0");

    m.def("solve", &py_solve,
          "Solve expr = 0",
          py::arg("expr"), py::arg("var") = "x");

    m.def("solve_system", &py_solve_system,
          "Solve a system of linear equations",
          py::arg("equations"), py::arg("variables"));

    m.def("evaluate", &py_evaluate,
          "Numerically evaluate an expression",
          py::arg("expr"), py::arg("values"));

    m.def("batch_evaluate", &py_batch_evaluate,
          "Batch evaluate expressions over a range (for graphing)",
          py::arg("exprs"), py::arg("var") = "x",
          py::arg("xmin") = -10.0, py::arg("xmax") = 10.0, py::arg("steps") = 400);

    m.def("vector_field_grid", &py_vector_field_grid,
          "Compute 2D vector field grid for visualization",
          py::arg("P_expr"), py::arg("Q_expr"),
          py::arg("xmin") = -5.0, py::arg("xmax") = 5.0,
          py::arg("ymin") = -5.0, py::arg("ymax") = 5.0,
          py::arg("grid_n") = 20);

    m.def("simplify", &py_simplify, "Simplify an expression", py::arg("expr"));

    m.def("to_latex", &py_to_latex, "Convert expression to LaTeX", py::arg("expr"));

    m.def("taylor_series", &py_taylor_series,
          "Compute Taylor/Maclaurin series expansion",
          py::arg("expr"), py::arg("var") = "x", py::arg("a") = 0.0, py::arg("terms") = 5);

    m.def("laplace_transform", &py_laplace_transform,
          "Compute Forward Laplace Transform",
          py::arg("expr"), py::arg("t_var") = "t", py::arg("s_var") = "s");

    m.def("inverse_laplace_transform", &py_inverse_laplace_transform,
          "Compute Inverse Laplace Transform",
          py::arg("expr"), py::arg("s_var") = "s", py::arg("t_var") = "t");

    m.def("compute_residue", &py_compute_residue,
          "Compute residue of a complex function",
          py::arg("expr"), py::arg("var") = "z", py::arg("real_part"), py::arg("imag_part"));

    m.def("contour_integral_circle", &py_contour_integral_circle,
          "Compute contour integral over a circle in the complex plane",
          py::arg("expr"), py::arg("var") = "z", py::arg("cx"), py::arg("cy"), py::arg("radius"), py::arg("steps") = 1000);

    m.def("fourier_transform", &py_fourier_transform,
          "Compute Continuous Fourier Transform",
          py::arg("expr"), py::arg("t_var") = "t", py::arg("w_var") = "w");

    m.def("line_integral", &py_line_integral,
          "Compute vector line integral",
          py::arg("F_exprs"), py::arg("path_exprs"), py::arg("t_var") = "t", py::arg("t_start") = 0.0, py::arg("t_end") = 1.0, py::arg("steps") = 1000);

    m.def("surface_integral", &py_surface_integral,
          "Compute vector surface integral",
          py::arg("F_exprs"), py::arg("surface_exprs"), py::arg("uv_vars"), py::arg("u_start") = 0.0, py::arg("u_end") = 1.0, py::arg("v_start") = 0.0, py::arg("v_end") = 1.0, py::arg("steps") = 100);
}
