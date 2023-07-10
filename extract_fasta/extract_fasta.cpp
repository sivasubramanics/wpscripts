#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: extract_fasta genome.fasta chr2:13242-143412" << std::endl;
        return 1;
    }

    std::string fasta_file = argv[1];
    std::string region = argv[2];
    std::string identifier = region.substr(0, region.find(':'));
    std::string range = region.substr(region.find(':') + 1);
    int start = std::stoi(range.substr(0, range.find('-')));
    int end = std::stoi(range.substr(range.find('-') + 1));

    std::ifstream file(fasta_file);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << fasta_file << std::endl;
        return 1;
    }

    bool found = false;
    bool in_range = false;
    int pos = 0;
    std::string line;
    std::ostringstream sequence;

    while (std::getline(file, line)) {
        if (line[0] == '>') { // Header line
            if (found) break; // Already found the target sequence
            if (line.substr(1) == identifier) found = true;
        }
        else if (found) { // Sequence line
            for (char c : line) {
                pos++;
                if (pos == start) in_range = true;
                if (in_range) sequence << c;
                if (pos == end) in_range = false;
            }
        }
        if (pos > end) break;
    }

    if (!found) {
        std::cerr << "Failed to find sequence with identifier: " << identifier << std::endl;
        return 1;
    }
    if (pos < start || pos < end) {
        std::cerr << "Warning: specified range extends beyond the end of the sequence" << std::endl;
    }

    std::cout << ">" << region << std::endl;
    std::cout << sequence.str() << std::endl;

    return 0;
}

