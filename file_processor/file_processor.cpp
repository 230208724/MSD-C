#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

#include "file_processor.h"

std::vector<std::string> splitString(const std::string& input, const std::string& delimiter=' ') {
    std::cout << "split " << input << "into ";

    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = input.find(delimiter);

    while (end != std::string::npos) {
        std::string token = input.substr(start, end - start);
        tokens.push_back(token);
        std::cout << token <<" ";
        start = end + delimiter.length();
        end = input.find(delimiter, start);
    }

    //output
    std::string lastToken = input.substr(start);
    tokens.push_back(lastToken);
    std::cout << lasttToken << " ";
    std::cout << std::endl;

    return tokens;
}

std::vector<std::string> getFilesInDirectory(const std::string& directory) {
    std::vector<std::string> files;
    std::ifstream dirStream(directory.c_str());
    if (dirStream) {
        std::string file;
        while (dirStream >> file) {
            if (file.back() == '0') {
                files.push_back(file.substr(0, file.length() - 1));
            }
        }
        dirStream.close();
    }
    return files;
}

std::string convertIntToString(int number, int width=-1) {
    if width == -1:
        return number.to_string()
    else:
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(width) << number;
        return oss.str();
}

std::string findBraindir(const std::string& brain, const std::string& modality="U19_Zeng") {
    std::cout << "brain " << brain << "has the directory ";
    std::string brain_dir = '';

    if (modality == "U19_Zeng") {
        std::vector<std::string> brain_dirs = {};
        std::string modality_dir = "/PBshare/TeraconvertedBrain";
        
        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find("mouse" + brain + "_teraconvert") != std::string::npos) {
                brain_dirs.push_back(entry.path().string());
            }
        }

        if (brain_dirs.empty()) {
            modality_dir = "/PBshare/BrainRaw/Unzipped_Brains/";
            for (const auto& entry : fs::directory_iterator(modality_dir)) {
                std::string filename = entry.path().filename().string();

                if (filename.find("mouse" + brain + "_teraconvert") != std::string::npos) {
                    brain_dirs.push_back(entry.path().string());
                }
            }
        }

        if (!brain_dirs.empty()) {
            brain_dir = brain_dirs[0];
        }
    }
    else if (modality == "U19_Huang") {
        std::string modality_dir = "/PBshare/Huang_Brains";

        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find("mouse" + brain) != std::string::npos) {
                brain_dir = entry.path().string();
                break;
            }
        }
    }
    else if (modality == "Huang") {
        std::string modality_dir = "/PBshare/Huang_Brains";

        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find(brain + "_JH_HK") != std::string::npos && filename.find("processed") != std::string::npos) {
                brain_dir = entry.path().string();
                break;
            }
        }
    }
    else if (modality == "Wu") {
        std::string modality_dir = "/PBshare/Zhuhao_Wu";

        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find(brain) != std::string::npos) {
                brain_dir = entry.path().string();
                break;
            }
        }
    }
    else if (modality == "Osten") {
        std::string modality_dir = "/PBshare/SEU-ALLEN/Projects/A2_A3_MouseBrain";
        std::string brain_type = "A2";
        if (brain.find("A3") != std::string::npos) {
            brain_type = "A3";
        }

        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find(brain_type) != std::string::npos && filename.find("terafly") != std::string::npos) {
                brain_dir = entry.path().string();
                break;
            }
        }
    }
    else { // Modality: Dong
        std::string modality_dir = "/PBshare/Dong"
        std::string modality_dir = "/PBshare/DongHW_Brains/20220315_SW220203_03_LS_6x_1000z";

        for (const auto& entry : fs::directory_iterator(modality_dir)) {
            std::string filename = entry.path().filename().string();

            if (filename.find(brain + "_TeraFly") != std::string::npos) {
                brain_dir = entry.path().string();
                break;
            }
        }
    }
    
    std::cout << brain_dir <<std::endl;
    return brain_dir;
}

std::vector<std::string> findResdir(const std::string& brain_dir) {
    std::cout << "brain dir " << brain_dir << "has subdir list ";

    std::vector<std::string> resdir_list;

    for (const auto& entry : fs::directory_iterator(brain_dir)) {
        std::string filename = entry.path().filename().string();
        std::cout << filename << " ";

        if (filename.substr(0, 3) == "RES") {
            resdir_list.push_back(entry.path().string());
        }
    }
    std::cout<<std::endl;

    std::sort(resdir_list.begin(), resdir_list.end(), [&](const std::string& a, const std::string& b) {
        int xsize_a = std::stoi(a.substr(a.find('x') + 1));
        int xsize_b = std::stoi(b.substr(b.find('x') + 1));
        return xsize_a < xsize_b;
    });

    std::cout<<"has final res list ";
    for (std::string res: resdir_list){
        std::cout << res << " ";
    }
    std::cout<<std::endl;

    return resdir_list;
}



int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "用法：" << argv[0] << " -i <输入字符串>" << std::endl;
        return 1;
    }

    std::string option = argv[1];
    std::string input = argv[2];

    if (option != "-i") {
        std::cout << "无效选项。请使用 -i 指定输入字符串。" << std::endl;
        return 1;
    }

    const std::string delimeter = '_';
    std::vector<std::string> splits = splitString(input,delimeter) 

    // 打印数值
    for (std::string var : splits) {
        std::cout << var << std::endl;
    }

    return 0;
}