#include "printing.hpp"
#include <sstream>
#include <cmath>
using namespace std;

namespace calc {

static int prec(Op op) {
    switch (op) {
        case Op::ADD: case Op::SUB: return 1;
        case Op::MUL: case Op::DIV: return 2;
        case Op::NEG: return 3;
        case Op::POW: return 4;
        default: return 5;
    }
}

static string fnName(Op op) {
    switch (op) {
        case Op::SIN: return "sin"; case Op::COS: return "cos"; case Op::TAN: return "tan";
        case Op::ASIN: return "asin"; case Op::ACOS: return "acos"; case Op::ATAN: return "atan";
        case Op::LN: return "ln"; case Op::LOG10: return "log"; case Op::EXP: return "exp";
        case Op::SQRT: return "sqrt"; case Op::ABS: return "abs";
        case Op::FRESNELS: return "FresnelS"; case Op::FRESNELC: return "FresnelC"; case Op::ERF: return "erf";
        default: return "?";
    }
}

static void toStringRec(const ExprPtr& e, ostringstream& out, int parentPrec) {
    switch (e->op) {
        case Op::NUM: out << numToStr(e->num); return;
        case Op::VAR: out << e->var; return;
        case Op::NEG: {
            out << "-";
            toStringRec(e->args[0], out, prec(Op::NEG));
            return;
        }
        case Op::ADD: case Op::SUB: case Op::MUL: case Op::DIV: case Op::POW: {
            int p = prec(e->op);
            bool wrap = p < parentPrec;
            string opstr = e->op == Op::ADD ? " + " : e->op == Op::SUB ? " - " :
                                e->op == Op::MUL ? "*" : e->op == Op::DIV ? "/" : "^";
            if (wrap) out << "(";
            toStringRec(e->args[0], out, p + (e->op == Op::SUB || e->op == Op::DIV ? 1 : 0));
            out << opstr;
            bool rightIsNeg = (e->args[1]->op == Op::NEG) ||
                              (e->args[1]->op == Op::NUM && e->args[1]->num < 0);
            if (rightIsNeg && e->op != Op::ADD) out << "(";
            toStringRec(e->args[1], out, p + (e->op == Op::POW ? 0 : 1));
            if (rightIsNeg && e->op != Op::ADD) out << ")";
            if (wrap) out << ")";
            return;
        }
        default: {
            out << fnName(e->op) << "(";
            toStringRec(e->args[0], out, 0);
            out << ")";
            return;
        }
    }
}

string toString(const ExprPtr& e) {
    ostringstream ss;
    toStringRec(e, ss, 0);
    return ss.str();
}

static string latexFnName(Op op) {
    switch (op) {
        case Op::SIN: return "\\sin";
        case Op::COS: return "\\cos";
        case Op::TAN: return "\\tan";
        case Op::ASIN: return "\\sin^{-1}";
        case Op::ACOS: return "\\cos^{-1}";
        case Op::ATAN: return "\\tan^{-1}";
        case Op::LN: return "\\ln";
        case Op::LOG10: return "\\log_{10}";
        case Op::EXP: return "e^";
        case Op::SQRT: return "\\sqrt";
        case Op::ABS: return "\\left|";
        case Op::FRESNELS: return "S";
        case Op::FRESNELC: return "C";
        case Op::ERF: return "\\text{erf}";
        default: return "?";
    }
}

static string latexNum(double v) {
    if (abs(v) < 1e-6) return "0";
    if (abs(v - M_PI) < 1e-9) return "\\pi";
    if (abs(v + M_PI) < 1e-9) return "-\\pi";
    if (abs(v - M_E) < 1e-9) return "e";
    if (abs(v + M_E) < 1e-9) return "-e";
    if (abs(v - round(v)) < 1e-9) {
        long long r = llround(v);
        return to_string(r);
    }
    ostringstream ss;
    ss.precision(10);
    ss << v;
    return ss.str();
}

static void toLatexRec(const ExprPtr& e, ostringstream& out, int parentPrec) {
    switch (e->op) {
        case Op::NUM:
            out << latexNum(e->num);
            return;

        case Op::VAR:
            out << e->var;
            return;

        case Op::NEG: {
            bool wrapChild = (e->args[0]->op == Op::ADD || e->args[0]->op == Op::SUB);
            out << "-";
            if (wrapChild) out << "\\left(";
            toLatexRec(e->args[0], out, prec(Op::NEG));
            if (wrapChild) out << "\\right)";
            return;
        }

        case Op::ADD: {
            int p = prec(Op::ADD);
            bool wrap = p < parentPrec;
            if (wrap) out << "\\left(";
            toLatexRec(e->args[0], out, p);
            out << " + ";
            toLatexRec(e->args[1], out, p);
            if (wrap) out << "\\right)";
            return;
        }

        case Op::SUB: {
            int p = prec(Op::SUB);
            bool wrap = p < parentPrec;
            if (wrap) out << "\\left(";
            toLatexRec(e->args[0], out, p);
            out << " - ";
            toLatexRec(e->args[1], out, p + 1);
            if (wrap) out << "\\right)";
            return;
        }

        case Op::MUL: {
            int p = prec(Op::MUL);
            bool wrap = p < parentPrec;
            if (wrap) out << "\\left(";

            bool leftIsNum = (e->args[0]->op == Op::NUM);
            bool rightIsNum = (e->args[1]->op == Op::NUM);

            toLatexRec(e->args[0], out, p);

            if (leftIsNum && rightIsNum) {
                out << " \\cdot ";
            } else if (leftIsNum) {

                out << " ";
            } else {
                out << " \\cdot ";
            }

            toLatexRec(e->args[1], out, p + 1);
            if (wrap) out << "\\right)";
            return;
        }

        case Op::DIV: {

            out << "\\frac{";
            toLatexRec(e->args[0], out, 0);
            out << "}{";
            toLatexRec(e->args[1], out, 0);
            out << "}";
            return;
        }

        case Op::POW: {

            if (e->args[1]->op == Op::DIV &&
                isNum(e->args[1]->args[0], 1) && isNum(e->args[1]->args[1], 2)) {
                out << "\\sqrt{";
                toLatexRec(e->args[0], out, 0);
                out << "}";
                return;
            }

            bool baseNeedsWrap = (e->args[0]->op == Op::ADD || e->args[0]->op == Op::SUB ||
                                  e->args[0]->op == Op::MUL || e->args[0]->op == Op::DIV ||
                                  e->args[0]->op == Op::NEG || e->args[0]->op == Op::POW);
            if (baseNeedsWrap) out << "\\left(";
            toLatexRec(e->args[0], out, prec(Op::POW) + 1);
            if (baseNeedsWrap) out << "\\right)";

            out << "^{";
            toLatexRec(e->args[1], out, 0);
            out << "}";
            return;
        }

        case Op::SQRT: {
            out << "\\sqrt{";
            toLatexRec(e->args[0], out, 0);
            out << "}";
            return;
        }

        case Op::ABS: {
            out << "\\left|";
            toLatexRec(e->args[0], out, 0);
            out << "\\right|";
            return;
        }

        case Op::EXP: {
            out << "e^{";
            toLatexRec(e->args[0], out, 0);
            out << "}";
            return;
        }

        case Op::SIN: case Op::COS: case Op::TAN:
        case Op::ASIN: case Op::ACOS: case Op::ATAN:
        case Op::LN: case Op::LOG10:
        case Op::FRESNELS: case Op::FRESNELC: case Op::ERF: {
            out << latexFnName(e->op) << "\\!\\left(";
            toLatexRec(e->args[0], out, 0);
            out << "\\right)";
            return;
        }

        default:
            out << "?";
            return;
    }
}

string toLatex(const ExprPtr& e) {
    ostringstream ss;
    toLatexRec(e, ss, 0);
    return ss.str();
}

}
