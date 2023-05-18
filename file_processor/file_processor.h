#ifndef FILEPROCESSOR_H
#define FILEPROCESSOR_H

#include <string>
#include <vector>

bool fileExists(const std::string& filename);
std::vector<std::string> splitString(const std::string& input, const std::string& delimiter);
std::vector<std::string> getFilesInDirectory(const std::string& directory);
//std::string findBraindir(const std::string& brain, const std::string& modality);
std::vector<std::string> findResdir(const std::string& brain_dir);
std::string convertIntToString(int number, int width);

#endif
