#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "gmres.h"
#include "matrix_loader.h"

using namespace Eigen;     // or using Eigen::MatrixXd; 
using namespace std;


int main()
{
	MatrixXd A;
	MatrixLoader("./matrix_collection/494_bus.mtx",&A);
	MatrixXf A_float = A.cast<float>();
	//cout << A << endl;
	VectorXd b = VectorXd::Zero(A.rows());
	VectorXf b_float = b.cast<float>();
	VectorXd x0 = VectorXd::Zero(A.rows());x0(0) = 1.0;
	VectorXf x0_float = x0.cast<float>();
	VectorXd xm = VectorXd::Zero(A.rows());
	VectorXf xm_float = xm.cast<float>();

	int restart_m = 100;
	double tol_double = 1e-10;float tol_float = 1e-10; 
	int inner_iteration_counts = 0;

	GmresUnpreconditionedInFloat(&A_float, &b_float, &x0_float, restart_m, tol_float, &xm_float, inner_iteration_counts);
	//GmresUnpreconditionedInMixedPrecision(&A, &A_float, &b, &x0, restart_m, tol_double, &xm, inner_iteration_counts);
	//GmresUnpreconditionedInDouble(&A, &b, &x0, restart_m, tol_double, &xm, inner_iteration_counts);
	cout << "xm = " << xm_float << endl;
	cout << "inner_iteration_counts = " << inner_iteration_counts << endl;
	system("pause");
	return 0;
}
