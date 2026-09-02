#include "fourier.hpp"
#include "parser.hpp"
#include "printing.hpp"

namespace calc {

using namespace std;

string fourier_transform(const string& expr_str, const string& t_var, const string& w_var) {
    auto f = parse(expr_str);
    if (!f) return "Error parsing expression";
    
    // Very basic pattern matching for Fourier transforms
    if (f->op == Op::NUM) {
        if (f->num == 1.0) return "2*pi*delta(" + w_var + ")";
        if (f->num == 0.0) return "0";
        return to_string(f->num) + " * 2*pi*delta(" + w_var + ")";
    }
    
    if (f->op == Op::VAR) {
        // F(t) = i*pi*delta'(w)
        return "Unsupported Fourier transform for plain " + f->var;
    }
    
    if (f->op == Op::SIN) {
        auto arg = f->args[0];
        if (arg->op == Op::VAR && arg->var == t_var) {
            return "i*pi*(delta(" + w_var + "+1) - delta(" + w_var + "-1))";
        }
        if (arg->op == Op::MUL) {
            if (arg->args[0]->op == Op::NUM) {
                string a = to_string(arg->args[0]->num);
                return "i*pi*(delta(" + w_var + "+" + a + ") - delta(" + w_var + "-" + a + "))";
            }
        }
    }
    
    if (f->op == Op::COS) {
        auto arg = f->args[0];
        if (arg->op == Op::VAR && arg->var == t_var) {
            return "pi*(delta(" + w_var + "-1) + delta(" + w_var + "+1))";
        }
        if (arg->op == Op::MUL) {
            if (arg->args[0]->op == Op::NUM) {
                string a = to_string(arg->args[0]->num);
                return "pi*(delta(" + w_var + "-" + a + ") + delta(" + w_var + "+" + a + "))";
            }
        }
    }
    
    if (f->op == Op::EXP) {
        auto arg = f->args[0];
        if (arg->op == Op::NEG) {
            auto inner = arg->args[0];
            if (inner->op == Op::VAR && inner->var == t_var) {
                return "1/(1 + i*" + w_var + ")";
            }
            if (inner->op == Op::MUL) {
                if (inner->args[0]->op == Op::NUM) {
                    string a = to_string(inner->args[0]->num);
                    return "1/(" + a + " + i*" + w_var + ")";
                }
            }
        }
    }
    
    return "Fourier Transform not implemented for this expression pattern.";
}

}
