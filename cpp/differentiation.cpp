#include "differentiation.hpp"
using namespace std;

namespace calc {

ExprPtr differentiate(const ExprPtr& e, const string& v) {
    switch (e->op) {
        case Op::NUM: return num(0);
        case Op::VAR: return num(e->var == v ? 1 : 0);
        case Op::ADD: return add(differentiate(e->args[0], v), differentiate(e->args[1], v));
        case Op::SUB: return sub(differentiate(e->args[0], v), differentiate(e->args[1], v));
        case Op::NEG: return neg(differentiate(e->args[0], v));

        case Op::MUL: {

            auto& a = e->args[0]; auto& b = e->args[1];
            return add(mul(differentiate(a, v), b), mul(a, differentiate(b, v)));
        }

        case Op::DIV: {

            auto& a = e->args[0]; auto& b = e->args[1];
            return divi(
                sub(mul(differentiate(a, v), b), mul(a, differentiate(b, v))),
                powr(b, num(2))
            );
        }

        case Op::POW: {
            auto& a = e->args[0]; auto& b = e->args[1];
            bool aHas = contains(a, v), bHas = contains(b, v);
            if (!bHas) {

                return mul(mul(b, powr(a, sub(b, num(1)))), differentiate(a, v));
            } else if (!aHas) {

                return mul(mul(powr(a, b), fn(Op::LN, a)), differentiate(b, v));
            } else {

                return mul(powr(a, b),
                    add(mul(differentiate(b, v), fn(Op::LN, a)),
                        mul(b, divi(differentiate(a, v), a))));
            }
        }

        case Op::SIN:
            return mul(fn(Op::COS, e->args[0]), differentiate(e->args[0], v));
        case Op::COS:
            return neg(mul(fn(Op::SIN, e->args[0]), differentiate(e->args[0], v)));
        case Op::TAN:
            return divi(differentiate(e->args[0], v),
                       powr(fn(Op::COS, e->args[0]), num(2)));

        case Op::ASIN:
            return divi(differentiate(e->args[0], v),
                       fn(Op::SQRT, sub(num(1), powr(e->args[0], num(2)))));
        case Op::ACOS:
            return neg(divi(differentiate(e->args[0], v),
                           fn(Op::SQRT, sub(num(1), powr(e->args[0], num(2))))));
        case Op::ATAN:
            return divi(differentiate(e->args[0], v),
                       add(num(1), powr(e->args[0], num(2))));

        case Op::LN:
            return divi(differentiate(e->args[0], v), e->args[0]);
        case Op::LOG10:
            return divi(differentiate(e->args[0], v),
                       mul(e->args[0], num(log(10.0))));
        case Op::EXP:
            return mul(fn(Op::EXP, e->args[0]), differentiate(e->args[0], v));

        case Op::SQRT:
            return divi(differentiate(e->args[0], v),
                       mul(num(2), fn(Op::SQRT, e->args[0])));
        case Op::ABS:
            return mul(differentiate(e->args[0], v),
                       divi(e->args[0], fn(Op::ABS, e->args[0])));

        case Op::FRESNELS:
            return mul(fn(Op::SIN, mul(divi(num(M_PI), num(2)), powr(e->args[0], num(2)))), differentiate(e->args[0], v));
        case Op::FRESNELC:
            return mul(fn(Op::COS, mul(divi(num(M_PI), num(2)), powr(e->args[0], num(2)))), differentiate(e->args[0], v));
        case Op::ERF:
            return mul(mul(divi(num(2), num(sqrt(M_PI))), fn(Op::EXP, neg(powr(e->args[0], num(2))))), differentiate(e->args[0], v));

        default: return num(0);
    }
}

ExprPtr differentiateOrder(const ExprPtr& e, const string& v, int order) {
    if (order <= 0) return e;
    ExprPtr result = e;
    for (int i = 0; i < order; i++) {
        result = simplify(differentiate(result, v));
    }
    return result;
}

ExprPtr mixedPartial(const ExprPtr& e, const vector<string>& vars) {
    if (vars.empty()) return e;
    ExprPtr result = e;
    for (const auto& v : vars) {
        result = simplify(differentiate(result, v));
    }
    return result;
}

}
