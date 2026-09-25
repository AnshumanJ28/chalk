#pragma once
#include "ast.hpp"
#include <string>
#include <vector>
#include <utility>
using namespace std;

namespace calc {

struct ODEResult {
    bool ok;
    string symbolic_solution;
    string latex_solution;
    string method;
    string message;
};

ODEResult solveODE(const string& dydx_expr);

ODEResult solveODE2(double a, double b, double c);

struct ODETrajectory {
    bool ok;
    vector<pair<double, double>> points;
    string message;
};

ODETrajectory computeTrajectory(const string& dydx_expr,
                                double x0, double y0,
                                double xmin, double xmax,
                                int steps = 500);

ODETrajectory computeTrajectory2(double a, double b, double c,
                                 const std::string& d2ydx2_expr,
                                 double x0, double y0, double v0,
                                 double xmin, double xmax,
                                 int steps = 500);

struct SlopeFieldPoint {
    double x, y;
    double slope;
    bool valid;
};

struct SlopeFieldData {
    bool ok;
    vector<SlopeFieldPoint> points;
    string message;
};

SlopeFieldData computeSlopeField(const string& dydx_expr,
                                 double xmin, double xmax,
                                 double ymin, double ymax,
                                 int grid_nx = 25, int grid_ny = 25);

}
