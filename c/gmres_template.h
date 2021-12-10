#pragma once
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;     // using Eigen::MatrixXd; 
using namespace std;
// function declaration
template <typename T_matrix, typename T_verctor>
int GmresUnpreconditioned(const T_matrix *A, const T_verctor *b, const T_verctor *x0,
								const int restart_m, const double tol,
								T_verctor *xm);

// function definition
template <typename T_matrix, typename T_verctor>
int GmresUnpreconditioned(const T_matrix *A, const T_verctor *b, const T_verctor *x0,
								const int restart_m, const double tol,
								T_verctor *xm)
{
	int inner_iteration_counts = 0;
	int row_A = (*A).rows();
	
	T_matrix H = T_matrix::Zero(restart_m + 1, restart_m);
	T_matrix V = T_matrix::Zero(row_A, restart_m + 1);
	T_matrix Hm_bar;
	T_matrix Vm;
	T_matrix G;
	T_matrix Rm;
	T_matrix gm;

	T_verctor R(restart_m);
	T_verctor r0;
	T_verctor beta_e1;
	T_verctor ym;
	T_verctor rm;

	float beta = 0;
	float down = 0.0, c = 0.0, s = 0.0;
	int real_m = 0;
	int j = 0, i = 0, k = 0;
	int Rm_size = 0;
	while (1)
	{
		// Arnoldi Process
		r0 = (*b) - (*A) * (*x0);
		beta = r0.norm();
		V.col(0) = r0 / beta;
		for (j = 0; j < restart_m; j++)
		{
			inner_iteration_counts++;
			R = (*A) * V.col(j);
			for (i = 0; i < j; i++)
			{
				H(i, j) = R.transpose() * V.col(i);
				R = R - H(i, j) * V.col(i);
			}
			H(j + 1, j) = R.norm();
			if (abs(H(j + 1, j)) < 1e-10)
			{
				j++;// here j don't increment automatically;
				break;
			}
			else
				V.col(j + 1) = R / H(j + 1, j);
		}
		//cout << "j=" << j << endl;
		Hm_bar = H.topLeftCorner(j + 1, j);
		Vm = V.leftCols(j);
		//cout << Hm_bar << endl;
		//cout << "restart_m=" << restart_m << endl;




		//Givens Rotation
		real_m = Hm_bar.cols();
		beta_e1 = T_verctor::Zero(real_m + 1); beta_e1(0) = beta;
		//cout << beta_e1 << endl; cout << endl;
		for (k = 0; k < real_m; k++)
		{
			G = T_matrix::Identity(real_m + 1, real_m + 1);
			// double sqrt(double), otherwise float sqrtf(float)
			down = sqrt((Hm_bar(k, k) * Hm_bar(k, k) + Hm_bar(k + 1, k) * Hm_bar(k + 1, k)));
			s = Hm_bar(k + 1, k) / down;
			c = Hm_bar(k, k) / down;
			G.block(k, k, 2, 2) << c, s, -s, c;
			Hm_bar = G * Hm_bar;
			beta_e1 = G * beta_e1;
		}



		Rm = Hm_bar.topLeftCorner(real_m, real_m);
		gm = beta_e1.head(real_m);
		// BackwardUpperTriangular Process, Rm*ym = gm;
		Rm_size = Rm.cols();
		ym = T_verctor::Zero(Rm_size);
		ym(ym.size() - 1) = gm(gm.size() - 1) / Rm(Rm_size - 1, Rm_size - 1);
		for (k = Rm_size - 2; k > -1; k--)
		{
			ym(k) = (gm(k) - Rm.row(k).tail(Rm_size - k - 1).dot(ym.tail(Rm_size - k - 1))) / Rm(k, k);
		}




		// x = x0 + V * y
		// real_solution = inv(A) * b;
		// xm = x0 + Vm * ym;
		(*xm) = (*x0) + Vm * ym;
		//cout << Vm << endl; cout << endl;
		//cout << Rm << endl; cout << endl;
		//cout << gm << endl; cout << endl;
		//cout << ym << endl; cout << endl;

		// judge whether to restart
		rm = (*b) - (*A) * (*xm);
		//cout << rm.norm() << endl;
		if (rm.norm() <= tol)
		{
			//cout << "xm =\n" << xm << endl;
			break;
		}
			
		x0 = xm;
		
	}
	
	return inner_iteration_counts;

}


