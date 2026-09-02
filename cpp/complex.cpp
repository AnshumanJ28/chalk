#include "complex.hpp"
#include "parser.hpp"
#include <cmath>
#include <stdexcept>

namespace calc {

using namespace std;

complex<double> evaluate_complex(ExprPtr expr, const string& var, complex<double> val) {
    if (!expr) return complex<double>(0.0, 0.0);
    switch (expr->op) {
        case Op::NUM:
            return complex<double>(expr->num, 0.0);
        case Op::VAR:
            if (expr->var == var) return val;
            if (expr->var == "pi") return complex<double>(M_PI, 0.0);
            if (expr->var == "e") return complex<double>(M_E, 0.0);
            if (expr->var == "i" || expr->var == "j") return complex<double>(0.0, 1.0);
            return complex<double>(0.0, 0.0);
        case Op::ADD:
            return evaluate_complex(expr->args[0], var, val) + evaluate_complex(expr->args[1], var, val);
        case Op::SUB:
            return evaluate_complex(expr->args[0], var, val) - evaluate_complex(expr->args[1], var, val);
        case Op::MUL:
            return evaluate_complex(expr->args[0], var, val) * evaluate_complex(expr->args[1], var, val);
        case Op::DIV: {
            complex<double> r = evaluate_complex(expr->args[1], var, val);
            if (abs(r) < 1e-15) throw runtime_error("Division by zero in complex eval");
            return evaluate_complex(expr->args[0], var, val) / r;
        }
        case Op::POW:
            return pow(evaluate_complex(expr->args[0], var, val), evaluate_complex(expr->args[1], var, val));
        case Op::NEG:
            return -evaluate_complex(expr->args[0], var, val);
        case Op::SIN:
            return sin(evaluate_complex(expr->args[0], var, val));
        case Op::COS:
            return cos(evaluate_complex(expr->args[0], var, val));
        case Op::TAN:
            return tan(evaluate_complex(expr->args[0], var, val));
        case Op::EXP:
            return exp(evaluate_complex(expr->args[0], var, val));
        case Op::LN:
            return log(evaluate_complex(expr->args[0], var, val));
        case Op::LOG10:
            return log(evaluate_complex(expr->args[0], var, val)) / log(10.0);
        case Op::SQRT:
            return sqrt(evaluate_complex(expr->args[0], var, val));
        default:
            return complex<double>(0.0, 0.0);
    }
}

complex<double> compute_residue(const string& expr_str, const string& var, complex<double> z0) {
    auto f = parse(expr_str);
    
    // Res(f, z0) = lim_{z->z0} (z-z0) * f(z)
    // Approximate by taking z very close to z0.
    // E.g. z = z0 + 1e-8
    complex<double> dz(1e-8, 1e-8);
    complex<double> z_near = z0 + dz;
    
    complex<double> f_val = evaluate_complex(f, var, z_near);
    complex<double> res = dz * f_val;
    return res;
}

complex<double> contour_integral_circle(const string& expr_str, const string& var, complex<double> center, double radius, int steps) {
    auto f = parse(expr_str);
    complex<double> sum(0.0, 0.0);
    complex<double> i_comp(0.0, 1.0);
    
    double dt = 2.0 * 3.14159265358979323846 / steps;
    
    // Parameterization: z(t) = center + R*e^{it}
    // dz = i*R*e^{it} dt
    for (int k = 0; k < steps; ++k) {
        double t = k * dt;
        complex<double> e_it = exp(i_comp * complex<double>(t, 0.0));
        complex<double> z = center + complex<double>(radius, 0.0) * e_it;
        complex<double> dz = i_comp * complex<double>(radius, 0.0) * e_it * complex<double>(dt, 0.0);
        
        complex<double> f_val = evaluate_complex(f, var, z);
        sum += f_val * dz;
    }
    
    return sum;
}

}
