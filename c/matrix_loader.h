#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <type_traits>
using namespace Eigen;     // using Eigen::MatrixXd; 
using namespace std;
// function declaration
template <typename T>
int MatrixLoader(string matrix_path,T *matrix);

// function definition
template <typename T>
int MatrixLoader(string matrix_path,T *matrix)
{
    // matrix load from .mtx
	std::ifstream file(matrix_path);
	int num_row, num_col, num_lines;
	
	// Ignore comments headers
	while (file.peek() == '%') file.ignore(2048, '\n');

	// Read number of rows and columns
	file >> num_row >> num_col >> num_lines;

	// Create 2D array and fill with zeros
    
	
    *matrix = T::Zero(num_row, num_col);//T is a class name, and we can use functions in it using T::.
    

	// fill the matrix with data
	for (int l = 0; l < num_lines; l++)
	{
        T data = T::Zero(1,1); //as per T, we construct a matrix which just have 1 elements to reserve the value read from file.
		int row, col;
		file >> row >> col >> data(0,0);
        (*matrix)(row - 1,col - 1) = data(0,0);
        
		
	}

 	file.close();
    return 0;
}