#include "gmres.h"

int GmresUnpreconditionedInMixedPrecision(const MatrixXd *A, const MatrixXf *A_float, const VectorXd *b, const VectorXd *x0, 
								const int restart_m, const double tol,
								VectorXd *xm, int &inner_iteration_counts)
{
	inner_iteration_counts = 0;
	int row_A = (*A).rows();
	
	MatrixXf H = MatrixXf::Zero(restart_m + 1, restart_m);
	MatrixXf V = MatrixXf::Zero(row_A, restart_m + 1);
	MatrixXf Hm_bar;
	MatrixXf Vm;
	MatrixXf G;
	MatrixXf Rm;
	MatrixXf gm;

	VectorXf R(restart_m);
	VectorXd r0;// r0 is stored in double precison.
	VectorXf beta_e1;
	VectorXf um;
	VectorXf ym;
	VectorXd rm; // rm is stored in double precison.

	float beta = 0;// beta is stored in float precision.
	double down = 0.0, c = 0.0, s = 0.0;
	int real_m = 0;
	int j = 0, i = 0, k = 0;
	int Rm_size = 0;
	while (1)
	{
		// Arnoldi Process
		r0 = (*b) - (*A) * (*x0);// calculate in double precision.
		beta = r0.norm();
		V.col(0) = r0.cast<float>() / beta;
		for (j = 0; j < restart_m; j++)
		{
			inner_iteration_counts++;
			R = (*A_float) * V.col(j);
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
		beta_e1 = VectorXf::Zero(real_m + 1); beta_e1(0) = beta; // beta is in float precision.
		for (k = 0; k < real_m; k++)
		{
			G = MatrixXf::Identity(real_m + 1, real_m + 1);
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
		ym = VectorXf::Zero(Rm_size);
		ym(ym.size() - 1) = gm(gm.size() - 1) / Rm(Rm_size - 1, Rm_size - 1);
		for (k = Rm_size - 2; k > -1; k--)
		{
			ym(k) = (gm(k) - Rm.row(k).tail(Rm_size - k - 1).dot(ym.tail(Rm_size - k - 1))) / Rm(k, k);
		}


		// xm = x0 + Vm * ym;
		um = Vm * ym;
		(*xm) = (*x0) + um.cast<double>();// update x in double precison.
		

		// judge whether to restart
		rm = (*b) - (*A) * (*xm);
		if (rm.norm() <= tol)
		{
			break;
		}
			
		x0 = xm;
		
	}
	
	return 0;

}