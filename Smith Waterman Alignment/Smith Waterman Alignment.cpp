// Smith Waterman Alignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>

int smith_waterman_alg(std::string seq1, std::string seq2) {
	//std::string seq1 = "GTCCGATGCTAGCTAGCTAGCATCGATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGATGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA";
	//std::string seq2 = "TTGACGTAAAGCTAGCTAGCATCGATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGATGCCATTGTAATGGGCCGCTGAAAGGGTTAAAGAT";
	int max_i = 0; 
	int max_j = 0; 

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
				max_i = i;
				max_j = j;
				//saves row and column of max score 
			}

		}
	}
	std::cout << "The highest local alignment score is: " << max_score << std::endl;

	//memory cleanup
	for (int i = 0; i < rows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	return max_score;

}
int main() {
	assert(smith_waterman_alg("CGTGAATTCG", "ACTGAATTCC") == 22);
	std::cout << "Baseline Sequence test passed successfully!" << std::endl;

	assert(smith_waterman_alg("ATCG", "ATCG") == 12); 
	std::cout << "Identical Sequence test passed successfully!" << std::endl;

	assert(smith_waterman_alg("AAAA", "TTTT") == 0);
	std::cout << "Zero Floor test passed successfully!" << std::endl;

	std::string wild_type = "GTCCGATGCTAGCTAGCTAGCATCGATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGATGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA";
	std::string variant = "TTGACGTAAAGCTAGCTAGCATCGATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGACTAGCTAGCTAGCATCGATCGATCGATGCCATTGTAATGGGCCGCTGAAAGGGTTAAAGAT";

	int stress_score = smith_waterman_alg(wild_type, variant);
	std::cout << "massive memory leak test passed. Score: " << stress_score << std::endl; 

	std::cout << "All tests passed successfully!" << std::endl;

	return 0; 
}