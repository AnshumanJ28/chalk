#include "vector_calc.hpp"
#include "parser.hpp"
#include "printing.hpp"
#include "differentiation.hpp"
using namespace std;

namespace calc {

VectorCalcResult gradient(const string& expr, const vector<string>& vars) {
    VectorCalcResult res;
    res.ok = false;

    try {
        ExprPtr e = parse(expr);

        for (const auto& v : vars) {
            ExprPtr d = simplify(differentiate(e, v));
            res.components.push_back(toString(d));
            res.latex_components.push_back(toLatex(d));
        }

        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing gradient: ") + ex.what();
    }
    return res;
}

VectorCalcResult divergence(const vector<string>& componentExprs,
                            const vector<string>& vars) {
    VectorCalcResult res;
    res.ok = false;

    if (componentExprs.size() != vars.size()) {
        res.message = "Number of component expressions must match number of variables.";
        return res;
    }

    try {
        ExprPtr sum = num(0);

        for (size_t i = 0; i < vars.size(); i++) {
            ExprPtr comp = parse(componentExprs[i]);
            ExprPtr d = simplify(differentiate(comp, vars[i]));
            sum = add(sum, d);
        }

        sum = simplify(sum);
        res.scalar = toString(sum);
        res.latex_scalar = toLatex(sum);
        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing divergence: ") + ex.what();
    }
    return res;
}

VectorCalcResult curl(const vector<string>& componentExprs,
                      const vector<string>& vars) {
    VectorCalcResult res;
    res.ok = false;

    if (componentExprs.size() != 3 || vars.size() != 3) {
        res.message = "Curl requires exactly 3 component expressions and 3 variables.";
        return res;
    }

    try {
        ExprPtr Fx = parse(componentExprs[0]);
        ExprPtr Fy = parse(componentExprs[1]);
        ExprPtr Fz = parse(componentExprs[2]);

        const string& x = vars[0];
        const string& y = vars[1];
        const string& z = vars[2];

        ExprPtr cx = simplify(sub(
            simplify(differentiate(Fz, y)),
            simplify(differentiate(Fy, z))
        ));

        ExprPtr cy = simplify(sub(
            simplify(differentiate(Fx, z)),
            simplify(differentiate(Fz, x))
        ));

        ExprPtr cz = simplify(sub(
            simplify(differentiate(Fy, x)),
            simplify(differentiate(Fx, y))
        ));

        res.components = {toString(cx), toString(cy), toString(cz)};
        res.latex_components = {toLatex(cx), toLatex(cy), toLatex(cz)};
        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing curl: ") + ex.what();
    }
    return res;
}

VectorCalcResult laplacian(const string& expr, const vector<string>& vars) {
    VectorCalcResult res;
    res.ok = false;

    try {
        ExprPtr e = parse(expr);
        ExprPtr sum = num(0);

        for (const auto& v : vars) {

            ExprPtr d1 = simplify(differentiate(e, v));
            ExprPtr d2 = simplify(differentiate(d1, v));
            sum = add(sum, d2);
        }

        sum = simplify(sum);
        res.scalar = toString(sum);
        res.latex_scalar = toLatex(sum);
        res.ok = true;
    } catch (const exception& ex) {
        res.message = string("Error computing laplacian: ") + ex.what();
    }
    return res;
}

double line_integral(const vector<string>& F_exprs, const vector<string>& path_exprs, const string& t_var, double t_start, double t_end, int steps) {
    if (F_exprs.size() != 3 || path_exprs.size() != 3) return 0.0;
    auto Fx = parse(F_exprs[0]), Fy = parse(F_exprs[1]), Fz = parse(F_exprs[2]);
    auto px = parse(path_exprs[0]), py = parse(path_exprs[1]), pz = parse(path_exprs[2]);
    
    double dt = (t_end - t_start) / steps;
    double sum = 0.0;
    
    for(int i=0; i<steps; ++i) {
        double t = t_start + i * dt;
        double x_val = evaluate(px, {{t_var, t}});
        double y_val = evaluate(py, {{t_var, t}});
        double z_val = evaluate(pz, {{t_var, t}});
        
        double x_next = evaluate(px, {{t_var, t + dt}});
        double y_next = evaluate(py, {{t_var, t + dt}});
        double z_next = evaluate(pz, {{t_var, t + dt}});
        
        double dx = x_next - x_val;
        double dy = y_next - y_val;
        double dz = z_next - z_val;
        
        map<string, double> env = {{"x", x_val}, {"y", y_val}, {"z", z_val}};
        double fx = evaluate(Fx, env);
        double fy = evaluate(Fy, env);
        double fz = evaluate(Fz, env);
        
        sum += fx*dx + fy*dy + fz*dz;
    }
    return sum;
}

double surface_integral(const vector<string>& F_exprs, const vector<string>& surface_exprs, const vector<string>& uv_vars, double u_start, double u_end, double v_start, double v_end, int steps) {
    if (F_exprs.size() != 3 || surface_exprs.size() != 3 || uv_vars.size() != 2) return 0.0;
    auto Fx = parse(F_exprs[0]), Fy = parse(F_exprs[1]), Fz = parse(F_exprs[2]);
    auto sx = parse(surface_exprs[0]), sy = parse(surface_exprs[1]), sz = parse(surface_exprs[2]);
    
    string u_var = uv_vars[0], v_var = uv_vars[1];
    double du = (u_end - u_start) / steps;
    double dv = (v_end - v_start) / steps;
    double sum = 0.0;
    
    for(int i=0; i<steps; ++i) {
        double u = u_start + i * du;
        for(int j=0; j<steps; ++j) {
            double v = v_start + j * dv;
            
            map<string, double> env = {{u_var, u}, {v_var, v}};
            double x_val = evaluate(sx, env);
            double y_val = evaluate(sy, env);
            double z_val = evaluate(sz, env);
            
            // Compute tangents
            map<string, double> env_u = {{u_var, u+du}, {v_var, v}};
            double dxdu = (evaluate(sx, env_u) - x_val)/du;
            double dydu = (evaluate(sy, env_u) - y_val)/du;
            double dzdu = (evaluate(sz, env_u) - z_val)/du;
            
            map<string, double> env_v = {{u_var, u}, {v_var, v+dv}};
            double dxdv = (evaluate(sx, env_v) - x_val)/dv;
            double dydv = (evaluate(sy, env_v) - y_val)/dv;
            double dzdv = (evaluate(sz, env_v) - z_val)/dv;
            
            // Cross product normal
            double nx = dydu*dzdv - dzdu*dydv;
            double ny = dzdu*dxdv - dxdu*dzdv;
            double nz = dxdu*dydv - dydu*dxdv;
            
            map<string, double> xyz_env = {{"x", x_val}, {"y", y_val}, {"z", z_val}};
            double fx = evaluate(Fx, xyz_env);
            double fy = evaluate(Fy, xyz_env);
            double fz = evaluate(Fz, xyz_env);
            
            sum += (fx*nx + fy*ny + fz*nz) * du * dv;
        }
    }
    return sum;
}

}
