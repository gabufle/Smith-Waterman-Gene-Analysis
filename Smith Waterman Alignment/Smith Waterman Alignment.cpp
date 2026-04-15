// Smith Waterman Alignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <algorithm>

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

	

	// make the scoring system 
	//keep track of the max score
	int max_score = 0; 

	for (int i = 1; i < rows; i++) {
		for (int j = 1; j < cols; j++) {

			//Option 1: Diagonal match or mismatch
			int diagonal_score = matrix[i - 1][j - 1];
			if (seq2[i - 1] == seq1[j - 1]) {
				diagonal_score += match_score; 
			}
			else {
				diagonal_score += mismatch_pen; 
			}

			// options 2 and 3: Gap conditions
			int up_score = matrix[i - 1][j] + gap_pen; 
			int left_score = matrix[i][j - 1] + gap_pen; 

			// highest score
			matrix[i][j] = std::max({ 0, diagonal_score, up_score, left_score }); 

			//Highest total score overall grid
			if (matrix[i][j] > max_score) {
				max_score = matrix[i][j];
			}

		}
	}
	std::cout << "The highest local alignment score is: " << max_score << std::endl;

	//memory cleanup
	for (int i = 0; i < rows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	return 0;
}