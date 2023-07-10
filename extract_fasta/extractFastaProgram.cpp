#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct Region {
    std::string id;
    int start;
    int end;
};

std::vector<Region> parseRegionsFile(const std::string& path) {
    std::ifstream file(path);
    std::vector<Region> regions;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        Region region;
        lineStream >> region.id >> region.start >> region.end;
        regions.push_back(region);
    }

    return regions;
}

Region parseRegionString(const std::string& str) {
    std::istringstream lineStream(str);
    Region region;
    lineStream >> region.id >> region.start >> region.end;
    return region;
}

void parseFastaFile(const std::string& fastaFile, const std::string& outputFile, 
                    const std::vector<Region>& regions) {
    std::ifstream file(fastaFile);
    std::ofstream outFile(outputFile);
    std::vector<Region> regionsToSearch = regions; // Regions we're searching for
    std::string sequenceIdentifier;
    bool in_range = false;
    int posInSequence = 0; // Position in sequence
    int outputLineLength = 0; 

    std::string line; 
    while(std::getline(file, line)) {
        if (line[0] == '>') { // This is a sequence identifier
            sequenceIdentifier = line.substr(1); // Get the identifier from the line
        }
        else { // This is a sequence line
            for(char c : line) {
                posInSequence++;
                for(auto it = regionsToSearch.begin(); it != regionsToSearch.end();) {
                    Region &region = *it;
                    if(sequenceIdentifier == region.id) {
                        if(posInSequence == region.start) in_range = true;
                        if(in_range) {
                            outFile << c;
                            outputLineLength++;
                            if(outputLineLength % 60 == 0) outFile << '\n';
                        }
                        if(posInSequence == region.end) {
                            in_range = false;
                            it = regionsToSearch.erase(it); // remove it- we're done with this region
                            continue;
                        }
                    }
                    ++it;
                }
            }
        }
        if(regionsToSearch.empty()) break;
    }

    // If there are still items in regionsToSearch, we didn't find everything we needed to.
    if(!regionsToSearch.empty()) {
        std::cerr << "Failed to find sequence for the following regions:\n";
        for(const Region& region : regionsToSearch) {
            std::cerr << region.id << ": " << region.start << "-" << region.end << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: extract_fasta -i genome.fasta -r chr2:13242-143412 -o out.fasta" << std::endl;
        std::cerr << "or: extract_fasta -i genome.fasta -b in.bed -o out.fasta" << std::endl;
        return 1;
    }

    std::string fastaFile, outputFile;
    std::vector<Region> regions;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-i") {
            fastaFile = argv[++i];
        } else if (arg == "-r") {
            regions.push_back(parseRegionString(argv[++i]));
        } else if (arg == "-b") {
            regions = parseRegionsFile(argv[++i]);
        } else if (arg == "-o") {
            outputFile = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }
    
    parseFastaFile(fastaFile, outputFile, regions);

    return 0;
}

