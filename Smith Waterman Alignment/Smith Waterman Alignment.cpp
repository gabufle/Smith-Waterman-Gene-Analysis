#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <chrono>

struct SWParams {
	int match = 3;
	int mismatch = -3;
	int gap_open = -4;
	int gap_extend = -1;
};



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
		if (current_line.empty() || current_line[0] == '>') {
			continue;
		}

		full_seq += current_line;
	}

	file.close();
	return full_seq;

}




int smith_waterman_alg(std::string seq1, std::string seq2, SWParams params = SWParams{}) {

	int max_i = 0;
	int max_j = 0;

	int cols = seq1.length() + 1;
	int rows = seq2.length() + 1;

	//memory alocation 
	int** M = new int* [rows];
	int** X = new int* [rows];
	int** Y = new int* [rows];

	for (int i = 0; i < rows; i++) {
		M[i] = new int[cols];
		X[i] = new int[cols];
		Y[i] = new int[cols];
	}

	//initialize 
	const int NEG_INF = -1000000;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			M[i][j] = 0;
			X[i][j] = (j == 0) ? NEG_INF : 0;
			Y[i][j] = (i == 0) ? NEG_INF : 0;
		}
	}

	std::cout << "matrix successfully made in memory: " << rows << " rows by: " << cols << " columns" << std::endl;



	// make the scoring system 
	//keep track of the max score
	int max_score = 0;

	for (int i = 1; i < rows; i++) {
		for (int j = 1; j < cols; j++) {

			//Diagonal match or mismatch
			int sub = (seq2[i - 1] == seq1[j - 1]) ? params.match : params.mismatch;
			M[i][j] = std::max(0, std::max({ M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1] }) + sub);

			//Gap condition: between X and M matrices
			X[i][j] = std::max(0, std::max({ M[i][j - 1] + params.gap_open, X[i][j - 1] + params.gap_extend }));

			//Gap condition: between Y and M matrices
			Y[i][j] = std::max(0, std::max({ M[i - 1][j] + params.gap_open, Y[i - 1][j] + params.gap_extend }));

			int best = std::max({ M[i][j], X[i][j], Y[i][j] });
			if (best > max_score) {
				max_score = best;
				max_i = i;
				max_j = j;
			}
		}
	}
	std::cout << "The highest local alignment score is: " << max_score << std::endl;


	// --- PHASE 3: THE TRACEBACK ---
	std::string align1 = "";
	std::string align2 = "";
	std::string match_line = "";

	int i = max_i;
	int j = max_j;

	// Figure out which matrix the max score came from
	enum State { IN_M, IN_X, IN_Y };

	State current_state;
	if (M[max_i][max_j] >= X[max_i][max_j] && M[max_i][max_j] >= Y[max_i][max_j])
		current_state = IN_M;
	else if (X[max_i][max_j] >= Y[max_i][max_j])
		current_state = IN_X;
	else
		current_state = IN_Y;


	// Walk backward until we hit a 0
	while (i > 0 && j > 0 && std::max({ M[i][j], X[i][j], Y[i][j] }) > 0) {

		if (current_state == IN_M) {
			// Came diagonally — always a match/mismatch
			align1 += seq1[j - 1];
			align2 += seq2[i - 1];

			// Bring back the visual match logic!
			if (seq1[j - 1] == seq2[i - 1]) {
				match_line += "|";
			}
			else {
				match_line += "*";
			}

			// Figure out where M[i][j] came from
			int sub = (seq2[i - 1] == seq1[j - 1]) ? params.match : params.mismatch;
			int prev_best = M[i][j] - sub; // what the score was before adding sub

			if (prev_best == M[i - 1][j - 1])      current_state = IN_M;
			else if (prev_best == X[i - 1][j - 1]) current_state = IN_X;
			else                               current_state = IN_Y;

			i--; j--;

		}
		else if (current_state == IN_X) {
			// In X means gap in seq2 — consuming seq1, moving left
			align1 += seq1[j - 1];
			align2 += "-";
			match_line += " ";

			// Did we open a fresh gap from M, or extend an existing X gap?
			if (M[i][j - 1] + params.gap_open >= X[i][j - 1] + params.gap_extend)
				current_state = IN_M;
			else
				current_state = IN_X;

			j--;

		}
		else { // IN_Y
			// In Y means gap in seq1 — consuming seq2, moving up
			align1 += "-";
			align2 += seq2[i - 1];
			match_line += " "; // Add space for gaps

			// Did we open a fresh gap from M, or extend an existing Y gap?
			if (M[i - 1][j] + params.gap_open >= Y[i - 1][j] + params.gap_extend)
				current_state = IN_M;
			else
				current_state = IN_Y;

			i--;
		}
	}

	std::reverse(align1.begin(), align1.end());
	std::reverse(match_line.begin(), match_line.end()); // Reverse the new line
	std::reverse(align2.begin(), align2.end());

	std::cout << "\nOptimal Local Alignment Found:" << std::endl;
	std::cout << "WT:  " << align1 << std::endl;
	std::cout << "     " << match_line << std::endl;
	std::cout << "VAR: " << align2 << "\n" << std::endl;



	//memory cleanup
	for (int i = 0; i < rows; i++) {
		delete[] M[i];
		delete[] X[i];
		delete[] Y[i];
	}
	delete[] M;
	delete[] X;
	delete[] Y;;

	return max_score;

}

void run_unit_tests() {
	std::cout << "Running internal unit tests..." << std::endl;

	// Test 1: Baseline
	assert(smith_waterman_alg("CGTGAATTCG", "ACTGAATTCC") == 20);

	// Test 2: Identical
	assert(smith_waterman_alg("ATCG", "ATCG") == 12);

	// Test 3: Zero Floor
	assert(smith_waterman_alg("AAAA", "TTTT") == 0);

	std::cout << "All unit tests passed. Engine is mathematically stable.\n" << std::endl;
}

void run_integration_tests() {
	std::cout << "\n--- RUNNING INTEGRATION TEST ---" << std::endl;
	std::cout << "Attempting to load real sequences from local FASTA files..." << std::endl;

	std::string wild_type = read_fasta("C:\\Users\\Gabriel\\Desktop\\code\\my-shit\\Smith Waterman Alignment\\wt.fasta.txt");
	std::string variant = read_fasta("C:\\Users\\Gabriel\\Desktop\\code\\my-shit\\Smith Waterman Alignment\\variant.fasta.txt");

	if (wild_type.length() > 0 && variant.length() > 0) {
		std::cout << "Files loaded! Running Affine Smith-Waterman...\n" << std::endl;
		smith_waterman_alg(wild_type, variant);
	}
	else {
		std::cerr << "Integration Test Skipped: Could not find the files." << std::endl;
	}
}

int main(int argc, char* argv[]) {
	// developer testing bloc
	if (argc == 1) {
		std::cout << "Starting Developer Test Suite...\n" << std::endl;
		run_unit_tests();
		run_integration_tests();
		return 0;
	}

	// MODE 2: Web API Mode (Exactly 2 files passed in)
	if (argc == 3) {
		std::string wild_type = read_fasta(argv[1]);
		std::string variant = read_fasta(argv[2]);

		if (wild_type.length() == 0 || variant.length() == 0) {
			std::cerr << "ERROR: One or both files were empty or unreadable." << std::endl;
			return 1;
		}

		// 1. Start the clock!
		auto start_time = std::chrono::high_resolution_clock::now();

		// 2. Run the engine silently for the web app
		smith_waterman_alg(wild_type, variant);

		// 3. Stop the clock!
		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> elapsed = end_time - start_time;

		// 4. Calculate the matrix payload
		long long matrix_cells = (long long)wild_type.length() * (long long)variant.length();

		// 5. Print the Flex Stats
		std::cout << "\n================ ENGINE STATS ================\n";
		std::cout << "Sequence 1 Length : " << wild_type.length() << " bp\n";
		std::cout << "Sequence 2 Length : " << variant.length() << " bp\n";
		std::cout << "Matrix Cells Calc : " << matrix_cells << "\n";
		std::cout << "Execution Time    : " << elapsed.count() << " ms\n";
		std::cout << "==============================================\n";
		return 0;
	}
	// ERROR CATCHER: Wrong number of arguments
	std::cerr << "ERROR: Invalid usage." << std::endl;
	std::cerr << "Usage for Web API: aligner.exe <reference.fasta> <variant.fasta>" << std::endl;
	std::cerr << "Usage for Dev Test: aligner.exe (with no arguments)" << std::endl;
	return 1;



	return 0;
}


