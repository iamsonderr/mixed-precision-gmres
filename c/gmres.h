#pragma once
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;     // using Eigen::MatrixXd; 
using namespace std;
int GmresUnpreconditionedDouble(const MatrixXd *A, const VectorXd *b, const VectorXd *x0,
								const int restart_m, const double tol,
								VectorXd *xm);