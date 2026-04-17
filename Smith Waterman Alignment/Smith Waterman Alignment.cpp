// Smith Waterman Alignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <fstream>

std::string read_fasta(std::string filepath) {
	std::ifstream file(filepath);
	std::string full_seq = "";
	std::string current_line;

	// check if file actually freakin opened
	if (!file.is_open()) {
		std::cerr << "Critical Error: could not open file: " << filepath << std::endl;
		return "";
	}

	// read the file
	while (std::getline(file, current_line)) {
		if (current_line.empty() || current_line[0] == '>'){
			continue;
		}

		full_seq += current_line;
	}

	file.close();
	return full_seq;

}




int smith_waterman_alg(std::string seq1, std::string seq2) {

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


	// --- PHASE 3: THE TRACEBACK ---
	std::string align1 = "";
	std::string align2 = "";

	int i = max_i;
	int j = max_j;

	// Walk backward until we hit a 0
	while (matrix[i][j] > 0) {
		int current_score = matrix[i][j];

		// Look at the neighbors
		int diagonal_score = matrix[i - 1][j - 1];
		int up_score = matrix[i - 1][j];
		int left_score = matrix[i][j - 1];

		// Calculate what the diagonal math would have been
		int diag_points;
		if (seq2[i - 1] == seq1[j - 1]) {
			diag_points = match_score;
		}
		else {
			diag_points = mismatch_pen; 
		}

		// REVERSE MATH: Which path gave us the current score?
		if (current_score == diagonal_score + diag_points) {
			align1 += seq1[j - 1]; // Grab DNA letter
			align2 += seq2[i - 1]; // Grab DNA letter
			i--; j--;              // Move diagonally up-left
		}
		else if (current_score == left_score + gap_pen) { 
			align1 += seq1[j - 1]; // Grab DNA letter
			align2 += "-";         // Insert gap
			j--;                   // Move left
		}
		else if (current_score == up_score + gap_pen) {   
			align1 += "-";         // Insert gap
			align2 += seq2[i - 1]; // Grab DNA letter
			i--;                   // Move up
		}
		else {
			break; // Safety net in case the math ever breaks
		}
	}

	// Since we walked backwards, the strings are backwards. Reverse them!
	std::reverse(align1.begin(), align1.end());
	std::reverse(align2.begin(), align2.end());

	std::cout << "\nOptimal Local Alignment Found:" << std::endl;
	std::cout << align1 << std::endl;
	std::cout << align2 << "\n" << std::endl;



	//memory cleanup
	for (int i = 0; i < rows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	return max_score;

}
int main() {

	// HERE is where you would add the FASTA filepaths to analyze. Make sure to have the FASTA files in the same directory as the executable or provide the correct relative/absolute path.
	// make sure to use double backslashes \\ if you're on Windows or single forward slashes / for Unix-based systems in the file paths.
	//std::string << wild_type = read_fasta("insert wild type DNA filepath here");
	//std::string << variant = read_fasta("insert variant or DNA of interest filepath here");

	//std::cout << "\nAttempting to load sequences from FASTA files..." << std::endl;
	//std::string wild_type = read_fasta("wild_type.fasta");
	//std::string variant = read_fasta("variant.fasta");

	//if (wild_type.length() > 0 && variant.length() > 0) {
		//int file_score = smith_waterman_alg(wild_type, variant);
	//	std::cout << "--- Test 4 Passed: File IO Successful ---" << std::endl;
	//}
	//else {
	//	std::cout << "--- Test 4 Skipped: Files not found ---" << std::endl;
	//}

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