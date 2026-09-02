#include "parser.hpp"
#include <cctype>
#include <stdexcept>
using namespace std;

namespace calc {

namespace {

enum class Tok { NUMBER, IDENT, PLUS, MINUS, STAR, SLASH, CARET, LP, RP, COMMA, END };

struct Token {
    Tok type;
    double num = 0;
    string text;
};

struct Lexer {
    string s;
    size_t i = 0;

    Lexer(const string& src) : s(src) {}

    void skipWs() {
        while (i < s.size() && isspace((unsigned char)s[i])) i++;
    }

    Token next() {
        skipWs();
        if (i >= s.size()) return {Tok::END};

        char c = s[i];

        if (isdigit((unsigned char)c) || c == '.') {
            size_t start = i;
            while (i < s.size() && (isdigit((unsigned char)s[i]) || s[i] == '.')) i++;
            return {Tok::NUMBER, stod(s.substr(start, i - start)), ""};
        }

        if (isalpha((unsigned char)c) || c == '_') {
            size_t start = i;
            while (i < s.size() && (isalnum((unsigned char)s[i]) || s[i] == '_')) i++;
            return {Tok::IDENT, 0, s.substr(start, i - start)};
        }

        i++;
        switch (c) {
            case '+': return {Tok::PLUS};
            case '-': return {Tok::MINUS};
            case '*': return {Tok::STAR};
            case '/': return {Tok::SLASH};
            case '^': return {Tok::CARET};
            case '(': return {Tok::LP};
            case ')': return {Tok::RP};
            case ',': return {Tok::COMMA};
            default: throw runtime_error(string("Unexpected character: ") + c);
        }
    }
};

struct Parser {
    Lexer lex;
    Token cur;

    Parser(const string& s) : lex(s) { cur = lex.next(); }
    void adv() { cur = lex.next(); }

    bool startsPrimary() {
        return cur.type == Tok::NUMBER || cur.type == Tok::IDENT || cur.type == Tok::LP;
    }

    ExprPtr parseExpr() {
        ExprPtr e = parseTerm();
        while (cur.type == Tok::PLUS || cur.type == Tok::MINUS) {
            bool isPlus = cur.type == Tok::PLUS;
            adv();
            ExprPtr rhs = parseTerm();
            e = isPlus ? add(e, rhs) : sub(e, rhs);
        }
        return e;
    }

    ExprPtr parseTerm() {
        ExprPtr e = parseUnary();
        for (;;) {
            if (cur.type == Tok::STAR) { adv(); e = mul(e, parseUnary()); }
            else if (cur.type == Tok::SLASH) { adv(); e = divi(e, parseUnary()); }
            else if (startsPrimary()) { e = mul(e, parseUnary()); }
            else break;
        }
        return e;
    }

    ExprPtr parseUnary() {
        if (cur.type == Tok::MINUS) { adv(); return neg(parseUnary()); }
        if (cur.type == Tok::PLUS) { adv(); return parseUnary(); }
        return parsePow();
    }

    ExprPtr parsePow() {
        ExprPtr base = parsePrimary();
        if (cur.type == Tok::CARET) {
            adv();
            ExprPtr exp = parseUnary();
            return powr(base, exp);
        }
        return base;
    }

    static Op fnOp(const string& name, bool& ok) {
        ok = true;
        if (name == "sin") return Op::SIN;
        if (name == "cos") return Op::COS;
        if (name == "tan") return Op::TAN;
        if (name == "asin" || name == "arcsin") return Op::ASIN;
        if (name == "acos" || name == "arccos") return Op::ACOS;
        if (name == "atan" || name == "arctan") return Op::ATAN;
        if (name == "ln") return Op::LN;
        if (name == "log") return Op::LOG10;
        if (name == "exp") return Op::EXP;
        if (name == "sqrt") return Op::SQRT;
        if (name == "abs") return Op::ABS;
        if (name == "FresnelS" || name == "S") return Op::FRESNELS;
        if (name == "FresnelC" || name == "C") return Op::FRESNELC;
        if (name == "erf") return Op::ERF;
        ok = false;
        return Op::NUM;
    }

    bool tryInverseTrigPattern(const string& name, Op& resultOp) {
        if ((name == "sin" || name == "cos" || name == "tan") &&
            cur.type == Tok::CARET) {

            size_t savedPos = lex.i;
            Token savedCur = cur;

            adv();

            bool isNegOne = false;

            if (cur.type == Tok::MINUS) {

                Token savedMinus = cur;
                size_t savedPos2 = lex.i;
                adv();
                if (cur.type == Tok::NUMBER && abs(cur.num - 1.0) < 1e-12) {
                    adv();
                    isNegOne = true;
                } else {

                    lex.i = savedPos;
                    cur = savedCur;
                    return false;
                }
            } else if (cur.type == Tok::LP) {

                size_t savedPos2 = lex.i;
                Token savedLP = cur;
                adv();
                if (cur.type == Tok::MINUS) {
                    adv();
                    if (cur.type == Tok::NUMBER && abs(cur.num - 1.0) < 1e-12) {
                        adv();
                        if (cur.type == Tok::RP) {
                            adv();
                            isNegOne = true;
                        } else {
                            lex.i = savedPos;
                            cur = savedCur;
                            return false;
                        }
                    } else {
                        lex.i = savedPos;
                        cur = savedCur;
                        return false;
                    }
                } else {
                    lex.i = savedPos;
                    cur = savedCur;
                    return false;
                }
            } else {

                lex.i = savedPos;
                cur = savedCur;
                return false;
            }

            if (isNegOne) {
                if (name == "sin") resultOp = Op::ASIN;
                else if (name == "cos") resultOp = Op::ACOS;
                else resultOp = Op::ATAN;
                return true;
            }
            return false;
        }
        return false;
    }

    ExprPtr parsePrimary() {
        if (cur.type == Tok::NUMBER) {
            double v = cur.num;
            adv();
            return num(v);
        }

        if (cur.type == Tok::LP) {
            adv();
            ExprPtr e = parseExpr();
            expect(Tok::RP);
            return e;
        }

        if (cur.type == Tok::IDENT) {
            string name = cur.text;
            adv();

            if (name == "pi") return num(M_PI);
            if (name == "e") return num(M_E);

            Op invOp;
            if (tryInverseTrigPattern(name, invOp)) {
                expect(Tok::LP);
                ExprPtr arg = parseExpr();
                expect(Tok::RP);
                return fn(invOp, arg);
            }

            bool ok;
            Op fo = fnOp(name, ok);
            if (ok) {
                expect(Tok::LP);
                ExprPtr arg = parseExpr();
                expect(Tok::RP);
                return fn(fo, arg);
            }

            return var(name);
        }

        throw runtime_error("Unexpected token in expression");
    }

    void expect(Tok t) {
        if (cur.type != t)
            throw runtime_error("Malformed expression (missing parenthesis or operand)");
        adv();
    }
};

}

ExprPtr parse(const string& text) {
    Parser p(text);
    ExprPtr e = p.parseExpr();
    if (p.cur.type != Tok::END)
        throw runtime_error("Unexpected trailing tokens in expression");
    return e;
}

}
