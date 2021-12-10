#pragma once
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;     // using Eigen::MatrixXd; 
using namespace std;

// function declaration
int GmresUnpreconditionedInFloat(const MatrixXf *A, const VectorXf *b, const VectorXf *x0, 
								const int restart_m, const float tol,
								VectorXf *xm, int &inner_iteration_counts);
						
int GmresUnpreconditionedInMixedPrecision(const MatrixXd *A, const MatrixXf *A_float, const VectorXd *b, const VectorXd *x0, 
								const int restart_m, const double tol,
								VectorXd *xm, int &inner_iteration_counts);

int GmresUnpreconditionedInDouble(const MatrixXd *A, const VectorXd *b, const VectorXd *x0,
								const int restart_m, const double tol,
								VectorXd *xm, int &inner_iteration_counts);



