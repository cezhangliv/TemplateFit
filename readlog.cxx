#include <iostream>
#include <fstream>
#include <string>

double SearchForPattern(const char* filename, const std::string& pattern) {
    std::ifstream file(filename);
    std::string line;
    int lineNumber = 0;

    double number = -1;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            lineNumber++;

            // Search for the pattern in the line
            size_t pos = line.find(pattern);
            if (pos != std::string::npos) {
                // Extract the number following the pattern
                pos += pattern.length(); // Move past the pattern
                size_t numStart = std::string::npos;
                size_t numEnd = std::string::npos;
                for (size_t i = pos; i < line.length(); ++i) {
                    if (std::isdigit(line[i]) || line[i] == '.') {
                        numStart = i;
                        break;
                    }
                }
                if (numStart != std::string::npos) {
                    for (size_t i = numStart; i < line.length(); ++i) {
                        if (!std::isdigit(line[i]) && line[i] != '.') {
                            numEnd = i;
                            break;
                        }
                    }
                }
                if (numStart != std::string::npos && numEnd != std::string::npos) {
                    std::string numberStr = line.substr(numStart, numEnd - numStart);
                    number = std::stod(numberStr);
                    std::cout << "Pattern found in line " << lineNumber << ": " << line << std::endl;
                    std::cout << "Extracted number: " << number << std::endl;
                }
            }
        }
        file.close();
        
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
    return number;
}

//void MyMacro() {
//    const char* filename = "data.txt"; // Change this to your file name
//    std::string pattern = "Weights normalization"; // Your pattern
//    SearchForPattern(filename, pattern);
//}


void readlog() {
    
    const char* filename = "/scratch/cezhang/Simulation_3cm_HitSharing1/files/exampleProductionJob_3cm_1708536731_328.root_/stat"; 
    
    std::string pattern = "Weights normalization"; 
    double weights = SearchForPattern(filename, pattern);
    
    std::string pattern = "CPU time to generate  ";
    double Nevents = SearchForPattern(filename, pattern);
    
}


