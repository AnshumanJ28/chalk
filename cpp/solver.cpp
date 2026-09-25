#include "solver.hpp"
#include "parser.hpp"
#include "printing.hpp"
#include "differentiation.hpp"
#include <cmath>
#include <set>
#include <algorithm>
#include <map>
using namespace std;

namespace calc {

SolveResult solve(const ExprPtr& e, const string& v) {
    ExprPtr se = simplify(e);
    SolveResult res;
    res.ok = false;

    ExprPtr d1 = simplify(differentiate(se, v));
    ExprPtr d2 = simplify(differentiate(d1, v));
    ExprPtr d3 = simplify(differentiate(d2, v));

    bool isPolyUpTo2 = (d3->op == Op::NUM && abs(d3->num) < 1e-9);

    map<string, double> zero{{v, 0.0}};

    if (isPolyUpTo2) {
        double c0, c1, c2;
        try {
            c0 = evaluate(se, zero);
            c1 = evaluate(d1, zero);
            c2 = evaluate(d2, zero) / 2.0;
        } catch (...) { c0 = c1 = c2 = NAN; }

        if (!isnan(c0) && !isnan(c1) && !isnan(c2)) {
            if (abs(c2) < 1e-9) {

                if (abs(c1) < 1e-9) {
                    res.ok = true;
                    if (abs(c0) < 1e-9) res.message = "Infinitely many solutions (identity).";
                    else res.message = "No solution.";
                    return res;
                }
                res.ok = true;
                string root = numToStr(-c0 / c1);
                res.roots.push_back(root);
                res.latex_roots.push_back(root);
                return res;
            } else {

                double disc = c1 * c1 - 4 * c2 * c0;
                res.ok = true;
                if (disc >= 0) {
                    double sq = sqrt(disc);
                    double r1 = (-c1 + sq) / (2 * c2);
                    double r2 = (-c1 - sq) / (2 * c2);
                    string sr1 = numToStr(r1);
                    string sr2 = numToStr(r2);
                    res.roots.push_back(sr1);
                    res.latex_roots.push_back(sr1);
                    if (abs(r1 - r2) > 1e-9) {
                        res.roots.push_back(sr2);
                        res.latex_roots.push_back(sr2);
                    }
                } else {
                    double re = -c1 / (2 * c2);
                    double im = sqrt(-disc) / (2 * c2);
                    string r1 = numToStr(re) + " + " + numToStr(im) + "i";
                    string r2 = numToStr(re) + " - " + numToStr(im) + "i";
                    string lr1 = numToStr(re) + " + " + numToStr(im) + "i";
                    string lr2 = numToStr(re) + " - " + numToStr(im) + "i";
                    res.roots.push_back(r1);
                    res.roots.push_back(r2);
                    res.latex_roots.push_back(lr1);
                    res.latex_roots.push_back(lr2);
                }
                return res;
            }
        }
    }

    set<long long> seen;
    vector<double> roots;
    for (double seed = -20.0; seed <= 20.0; seed += 0.5) {
        double x = seed;
        bool converged = false;
        for (int iter = 0; iter < 60; iter++) {
            map<string, double> m{{v, x}};
            double fx, dfx;
            try {
                fx = evaluate(se, m);
                dfx = evaluate(d1, m);
            } catch (...) { break; }
            if (!isfinite(fx) || !isfinite(dfx) || abs(dfx) < 1e-12) break;
            double nx = x - fx / dfx;
            if (abs(nx - x) < 1e-9) { x = nx; converged = true; break; }
            x = nx;
            if (abs(x) > 1e8) break;
        }
        if (converged) {
            map<string, double> m{{v, x}};
            double check;
            try { check = evaluate(se, m); } catch (...) { continue; }
            if (abs(check) < 1e-6) {
                long long key = llround(x * 1e4);
                if (!seen.count(key)) { seen.insert(key); roots.push_back(x); }
            }
        }
    }

    if (!roots.empty()) {
        sort(roots.begin(), roots.end());
        res.ok = true;
        for (double r : roots) {
            string s = numToStr(r);
            res.roots.push_back(s);
            res.latex_roots.push_back(s);
        }
        return res;
    }

    res.ok = false;
    res.message = "No real roots found in range [-20, 20]. The equation may have no real solutions, or roots outside this range.";
    return res;
}

static pair<string, string> splitEquation(const string& eqStr) {
    size_t eqPos = eqStr.find('=');
    if (eqPos == string::npos) return {eqStr, "0"};
    if (eqStr.find('=', eqPos + 1) != string::npos)
        throw runtime_error("Equation contains multiple '=' signs.");
    return {eqStr.substr(0, eqPos), eqStr.substr(eqPos + 1)};
}

static bool solveMatrix(vector<vector<double>>& A,
                        vector<double>& B,
                        vector<double>& X) {
    int N = A.size();
    X.assign(N, 0.0);

    for (int i = 0; i < N; i++) {
        int pivot = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[pivot][i])) pivot = j;
        }
        if (pivot != i) {
            swap(A[i], A[pivot]);
            swap(B[i], B[pivot]);
        }
        if (abs(A[i][i]) < 1e-12) return false;

        for (int j = i + 1; j < N; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) A[j][k] -= factor * A[i][k];
            B[j] -= factor * B[i];
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) sum += A[i][j] * X[j];
        X[i] = (B[i] - sum) / A[i][i];
    }
    return true;
}

SystemSolveResult solveSystem(const vector<string>& equations,
                              const vector<string>& variables) {
    SystemSolveResult res;
    res.ok = false;
    int N = equations.size();

    if (N == 0 || N != (int)variables.size()) {
        res.message = "Number of equations must match number of variables and be non-zero.";
        return res;
    }
    if (N > 4) {
        res.message = "System solver only supports up to 4 equations.";
        return res;
    }

    vector<vector<double>> A(N, vector<double>(N, 0.0));
    vector<double> B(N, 0.0);

    try {
        for (int i = 0; i < N; i++) {
            auto parts = splitEquation(equations[i]);
            ExprPtr lhs = parse(parts.first);
            ExprPtr rhs = parse(parts.second);
            ExprPtr eqExpr = simplify(sub(lhs, rhs));

            for (int j = 0; j < N; j++) {
                const auto& varName = variables[j];
                ExprPtr diffExpr = simplify(differentiate(eqExpr, varName));

                for (const auto& v : variables) {
                    if (contains(diffExpr, v)) {
                        res.message = "The equations must be linear in the variables.";
                        return res;
                    }
                }

                map<string, double> zeroVals;
                for (const auto& v : variables) zeroVals[v] = 0.0;
                double coeff = evaluate(diffExpr, zeroVals);
                A[i][j] = coeff;
            }

            map<string, double> zeroVals;
            for (const auto& v : variables) zeroVals[v] = 0.0;
            double c0 = evaluate(eqExpr, zeroVals);
            B[i] = -c0;
        }
    } catch (const exception& ex) {
        res.message = string("Error parsing/evaluating equations: ") + ex.what();
        return res;
    }

    vector<double> X(N, 0.0);
    if (!solveMatrix(A, B, X)) {
        res.message = "The system has no unique solution (matrix is singular or inconsistent).";
        return res;
    }

    res.ok = true;
    for (int j = 0; j < N; j++) {
        res.solution[variables[j]] = X[j];
    }
    return res;
}

}
