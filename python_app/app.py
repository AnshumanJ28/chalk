"""
Flask web app for the C++-powered symbolic calculator.

All the actual math (parsing, differentiation, integration, limits,
equation solving, evaluation, vector calculus, ODE solving) happens
in C++ (see cpp/), exposed to Python via the pybind11 module `calcengine`.
This file just wraps that engine and serves a small web UI.
"""
from flask import Flask, request, jsonify, render_template
import calcengine as ce

app = Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.config['TEMPLATES_AUTO_RELOAD'] = True

def clean_var(v):
    v = (v or "x").strip()
    return v if v else "x"

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/api/differentiate", methods=["POST"])
def api_differentiate():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    order = data.get("order", 1)
    vars_list = data.get("vars", [])
    try:
        r = ce.differentiate(expr, var, order, vars_list)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/integrate", methods=["POST"])
def api_integrate():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    order = int(data.get("order", 1))
    try:
        r = ce.integrate(expr, var)
        for i in range(1, order):
            if r.get("ok"):
                r = ce.integrate(r["result"], var)
            else:
                break
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/definite_integral", methods=["POST"])
def api_definite_integral():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    lower = float(data.get("lower", 0))
    upper = float(data.get("upper", 1))
    try:
        r = ce.definite_integral(expr, var, lower, upper)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/multi_integral", methods=["POST"])
def api_multi_integral():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    specs = data.get("specs", [])
    try:
        r = ce.multi_integral(expr, specs)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/vector_calc", methods=["POST"])
def api_vector_calc():
    data = request.get_json(force=True)
    op = data.get("op", "")
    scalar_expr = data.get("scalar_expr", "")
    vector_exprs = data.get("vector_exprs", [])
    vars_list = data.get("vars", ["x", "y", "z"])
    try:
        r = ce.vector_calc(op, scalar_expr, vector_exprs, vars_list)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/solve_ode", methods=["POST"])
def api_solve_ode():
    data = request.get_json(force=True)
    dydx_expr = data.get("expr", "")
    x0 = float(data.get("x0", 0))
    y0 = float(data.get("y0", 1))
    xmin = float(data.get("xmin", -5))
    xmax = float(data.get("xmax", 5))
    steps = int(data.get("steps", 500))
    try:
        r = ce.solve_ode(dydx_expr, x0, y0, xmin, xmax, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/limit", methods=["POST"])
def api_limit():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    point = (data.get("point", "0") or "0").strip()
    try:
        r = ce.limit(expr, var, point)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/solve", methods=["POST"])
def api_solve():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    try:
        r = ce.solve(expr, var)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "roots": [], "latex_roots": [], "message": f"Error: {ex}"})

@app.route("/api/solve_system", methods=["POST"])
def api_solve_system():
    data = request.get_json(force=True)
    equations = data.get("equations", [])
    variables = data.get("variables", [])
    try:
        variables = [v.strip() for v in variables if v.strip()]
        r = ce.solve_system(equations, variables)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/evaluate", methods=["POST"])
def api_evaluate():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    values = data.get("values", {}) or {}
    try:
        values = {k: float(v) for k, v in values.items()}
        r = ce.evaluate(expr, values)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": 0, "error": str(ex)})

@app.route("/api/batch_evaluate", methods=["POST"])
def api_batch_evaluate():
    data = request.get_json(force=True)
    exprs = data.get("exprs", [])
    var = data.get("var", "x")
    xmin = float(data.get("xmin", -10))
    xmax = float(data.get("xmax", 10))
    steps = int(data.get("steps", 400))
    try:
        r = ce.batch_evaluate(exprs, var, xmin, xmax, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/vector_field_grid", methods=["POST"])
def api_vector_field_grid():
    data = request.get_json(force=True)
    P_expr = data.get("P_expr", "")
    Q_expr = data.get("Q_expr", "")
    xmin = float(data.get("xmin", -5))
    xmax = float(data.get("xmax", 5))
    ymin = float(data.get("ymin", -5))
    ymax = float(data.get("ymax", 5))
    grid_n = int(data.get("grid_n", 20))
    try:
        r = ce.vector_field_grid(P_expr, Q_expr, xmin, xmax, ymin, ymax, grid_n)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/simplify", methods=["POST"])
def api_simplify():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    try:
        r = ce.simplify(expr)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "result": f"Error: {ex}", "latex": ""})

@app.route("/api/taylor_series", methods=["POST"])
def api_taylor_series():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "x"))
    a = float(data.get("a", 0))
    terms = int(data.get("terms", 5))
    try:
        r = ce.taylor_series(expr, var, a, terms)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/laplace", methods=["POST"])
def api_laplace():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    t_var = clean_var(data.get("t_var", "t"))
    s_var = clean_var(data.get("s_var", "s"))
    try:
        r = ce.laplace_transform(expr, t_var, s_var)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/inv_laplace", methods=["POST"])
def api_inv_laplace():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    s_var = clean_var(data.get("s_var", "s"))
    t_var = clean_var(data.get("t_var", "t"))
    try:
        r = ce.inverse_laplace_transform(expr, s_var, t_var)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/solve_ode2", methods=["POST"])
def api_solve_ode2():
    data = request.get_json(force=True)
    a = float(data.get("a", 1))
    b = float(data.get("b", 0))
    c = float(data.get("c", 0))
    expr = data.get("expr", "")
    x0 = float(data.get("x0", 0))
    y0 = float(data.get("y0", 1))
    v0 = float(data.get("v0", 0))
    xmin = float(data.get("xmin", -5))
    xmax = float(data.get("xmax", 5))
    steps = int(data.get("steps", 500))
    try:
        r = ce.solve_ode2(a, b, c, expr, x0, y0, v0, xmin, xmax, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": f"Error: {ex}"})

@app.route("/api/complex_residue", methods=["POST"])
def api_complex_residue():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "z"))
    real_part = float(data.get("real_part", 0))
    imag_part = float(data.get("imag_part", 0))
    try:
        r = ce.compute_residue(expr, var, real_part, imag_part)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": str(ex)})

@app.route("/api/complex_contour", methods=["POST"])
def api_complex_contour():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    var = clean_var(data.get("var", "z"))
    cx = float(data.get("cx", 0))
    cy = float(data.get("cy", 0))
    radius = float(data.get("radius", 1))
    steps = int(data.get("steps", 1000))
    try:
        r = ce.contour_integral_circle(expr, var, cx, cy, radius, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": str(ex)})

@app.route("/api/fourier", methods=["POST"])
def api_fourier():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    t_var = clean_var(data.get("t_var", "t"))
    w_var = clean_var(data.get("w_var", "w"))
    try:
        r = ce.fourier_transform(expr, t_var, w_var)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": str(ex)})

@app.route("/api/line_integral", methods=["POST"])
def api_line_integral():
    data = request.get_json(force=True)
    F_exprs = data.get("F_exprs", ["0", "0", "0"])
    path_exprs = data.get("path_exprs", ["0", "0", "0"])
    t_var = clean_var(data.get("t_var", "t"))
    t_start = float(data.get("t_start", 0))
    t_end = float(data.get("t_end", 1))
    steps = int(data.get("steps", 1000))
    try:
        r = ce.line_integral(F_exprs, path_exprs, t_var, t_start, t_end, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": str(ex)})

@app.route("/api/surface_integral", methods=["POST"])
def api_surface_integral():
    data = request.get_json(force=True)
    F_exprs = data.get("F_exprs", ["0", "0", "0"])
    surface_exprs = data.get("surface_exprs", ["0", "0", "0"])
    uv_vars = data.get("uv_vars", ["u", "v"])
    u_start = float(data.get("u_start", 0))
    u_end = float(data.get("u_end", 1))
    v_start = float(data.get("v_start", 0))
    v_end = float(data.get("v_end", 1))
    steps = int(data.get("steps", 100))
    try:
        r = ce.surface_integral(F_exprs, surface_exprs, uv_vars, u_start, u_end, v_start, v_end, steps)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "message": str(ex)})

@app.route("/api/to_latex", methods=["POST"])
def api_to_latex():
    data = request.get_json(force=True)
    expr = data.get("expr", "")
    try:
        r = ce.to_latex(expr)
        return jsonify(r)
    except Exception as ex:
        return jsonify({"ok": False, "latex": ""})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
