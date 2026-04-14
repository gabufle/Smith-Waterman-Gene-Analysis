// Smith Waterman Alignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>

int main() {
	std::string seq1 = "CGTGAATTCG";
	std::string seq2 = "ACTGAATTCC";

	int cols = seq1.length() + 1;
	int rows = seq2.length() + 1;

	int match_score = 3;
	int mismatch_pen = -3;
	int gap_pen = -2;

	//memory alocation 
	int** matrix = new int* [rows];

	for (int i = 0; i < rows; i++) {
		matrix[i] = new int[cols];
	}

	//initialize 
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			matrix[i][j] = 0;
		}
	}

	std::cout << "matrix successfully made in memory: " << rows << " rows by: " << cols << " columns" << std::endl;

	//memory cleanup
	for (int i = 0; i < rows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;


	return 0;
}