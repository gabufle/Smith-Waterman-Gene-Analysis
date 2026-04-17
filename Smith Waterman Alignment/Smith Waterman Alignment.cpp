#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <fstream>

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
		if (current_line.empty() || current_line[0] == '>'){
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
			X[i][j] = std::max(0, std::max({M[i][j - 1] + params.gap_open, X[i][j - 1] + params.gap_extend}));

			//Gap condition: between Y and M matrices
			Y[i][j] = std::max(0, std::max({M[i - 1][j] + params.gap_open, Y[i - 1][j] + params.gap_extend}));

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

	assert(smith_waterman_alg("CGTGAATTCG", "ACTGAATTCC") == 21);
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