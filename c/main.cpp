#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "gmres.h"
using namespace Eigen;     // or using Eigen::MatrixXd; 
using namespace std;
int main()
{
	// matrix load from .mtx
	std::ifstream file("./matrix_collection/494_bus.mtx");
	int num_row, num_col, num_lines;
	
	// Ignore comments headers
	while (file.peek() == '%') file.ignore(2048, '\n');

	// Read number of rows and columns
	file >> num_row >> num_col >> num_lines;
	
	// Create 2D array and fill with zeros
	MatrixXd A = MatrixXd::Zero(num_row, num_col);

	// fill the matrix with data
	for (int l = 0; l < num_lines; l++)
	{
		double data;
		int row, col;
		file >> row >> col >> data;
		A(row - 1,col - 1) = data;
	}

 	file.close();
	//cout << A << endl;

	//MatrixXd A = MatrixXd::Identity(3,3);
	VectorXd b = VectorXd::Zero(A.rows());
	VectorXd x0 = VectorXd::Zero(A.rows());
	VectorXd xm = VectorXd::Zero(A.rows());

	x0(0) = 1;
	int restart_m = 100;
	double tol = 1e-10; 

	int inner_iteration_counts = GmresUnpreconditionedDouble(&A, &b, &x0, restart_m, tol, &xm);
	cout << "xm = " << xm << endl;
	cout << "inner_iteration_counts = " << inner_iteration_counts << endl;
	system("pause");
	return 0;
}
