#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
//#include <filesystem>
#include <unordered_map>

#include "../config.h"
#include "../file_processor/file_processor.h"

#include "swc_basic_process.h"

void parseSwcfile(const std::string& swc_file, std::vector<Node>& tree) {
    std::ifstream fp(swc_file);
    std::string line;
    while (std::getline(fp, line)) {
        line = line.substr(0, line.find('#')); // Remove comments
        if (line.empty()) continue;

        std::stringstream ss(line);
        Node node;
        ss >> node.idx >> node.type >> node.x >> node.y >> node.z >> node.r >> node.parent;
        tree.push_back(node);
    }
}

bool checkNodestree(const std::string& swc_file, std::vector<Node>& tree){
    // Check if file exists
    if (!fileExists(swc_file)){
        std::cout << "File does not exist: " << swc_file << std::endl;
        return false;
    }

    // Parse SWC file
    parseSwcfile(swc_file, tree);
    if (tree.size() == 0){
        std::cout << "Failed to parse SWC file: " << swc_file << std::endl;
        return false;
    }

    std::cout << swc_file << " Nodes size is: " << tree.size() << std::endl;
    return true;
}

void writeNodestree(const std::vector<Node>& tree, const std::string& swc_file, const std::vector<std::string>& header) {
    std::ofstream fp(swc_file);
    for (const auto& s : header) {
        fp << s << "\n";
    }
    fp << "##n type x y z r parent\n";
    for (const Node& leaf : tree) {
        fp << leaf.idx << " " << leaf.type << " " << leaf.x << " " << leaf.y << " " << leaf.z << " " << leaf.r << " " << leaf.parent << "\n";
    }
}

std::vector<Marker> parseMarkerfile(const std::string& marker_file) {
    std::vector<Marker> tree;
    std::ifstream fp(marker_file);
    std::string line;
    while (std::getline(fp, line)) {
        line = line.substr(0, line.find('#')); // Remove comments
        
        if (line.empty()) continue;

        Marker marker;

        // Split the string into individual components
        std::vector<std::string> components;
        size_t startPos = 0;
        size_t endPos = line.find(",");
        while (endPos != std::string::npos) {
            components.push_back(line.substr(startPos, endPos - startPos));
            startPos = endPos + 1;
            endPos = line.find(",", startPos);
        }
        components.push_back(line.substr(startPos, line.length() - startPos));

        // Assign the components to the marker's variables
        size_t componentCount = components.size();
        if (componentCount >= 3) {
            marker.x = std::stod(components[0]);
            marker.y = std::stod(components[1]);
            marker.z = std::stod(components[2]);
        }
        if (componentCount >= 4) {
            marker.r = std::stod(components[3]);
        }
        if (componentCount >= 5) {
            marker.shape = std::stoi(components[4]);
        }
        if (componentCount >= 6) {
            marker.name = components[5];
        }
        if (componentCount >= 7) {
            marker.comment = components[6];
        }
        if (componentCount >= 8) {
            marker.red = std::stoi(components[7]);
        }
        if (componentCount >= 9) {
            marker.green = std::stoi(components[8]);
        }
        if (componentCount >= 10) {
            marker.blue = std::stoi(components[9]);
        }
        tree.push_back(marker);
    }
    return tree;
}

bool writeMarker(const Marker& marker, const std::string& markerfile, const std::vector<std::string>& header) {
    std::ofstream fp(markerfile);
    for (const auto& s : header) {
        fp << s << "\n";
    }
    fp << "##x,y,z,radius,shape,name,comment,color_r,color_g,color_b\n";
    
    fp << marker.x << "," << marker.y << "," << marker.z << "," << marker.r << "," << marker.shape << "," << marker.name << "," << marker.comment << "," << marker.red << "," << marker.green << "," << marker.blue<< std::endl;
    return true;
}

bool writeMarkerstree(const std::vector<Marker>& markerstree, const std::string& markerfile, const std::vector<std::string>& header) {
    std::ofstream fp(markerfile);
    for (const auto& s : header) {
        fp << s << "\n";
    }
    fp << "##x,y,z,radius,shape,name,comment,color_r,color_g,color_b\n";
    
    for (const Marker& marker : markerstree) {
         fp << marker.x << "," << marker.y << "," << marker.z << "," << marker.r << "," << marker.shape << "," << marker.name << "," << marker.comment << "," << marker.red << "," << marker.green << "," << marker.blue<< std::endl;
    }
    return true;
}

Node find_first_soma(const std::vector<Node>& tree, int p_soma) {
    for (const auto& leaf : tree) {
        if (leaf.parent == p_soma) {
            return leaf;
        }
    }
    throw std::runtime_error("Could not find the soma node!");
}

bool is_in_box(double x, double y, double z, const std::vector<int>& xyzshape) {
    if (x < 0 || y < 0 || z < 0 ||
        x > xyzshape[0] - 1 ||
        y > xyzshape[1] - 1 ||
        z > xyzshape[2] - 1) {
        return false;
    }
    return true;
}

bool is_in_bbox(double x, double y, double z, const std::vector<std::vector<int>>& xyzxyz) {
    int xmin = xyzxyz[0][0];
    int ymin = xyzxyz[0][1];
    int zmin = xyzxyz[0][2];
    int xmax = xyzxyz[1][0];
    int ymax = xyzxyz[1][1];
    int zmax = xyzxyz[1][2];
    if (x < xmin || y < ymin || z < zmin ||
        x > xmax ||
        y > ymax ||
        z > zmax) {
        return false;
    }
    return true;
}

std::vector<Node> crop_bbox_tree(const std::vector<Node>& tree, const std::vector<std::vector<int>>& xyzxyz) {
    std::vector<Node> new_tree;
    for (const auto& leaf : tree) {
        double x = leaf.x;
        double y = leaf.y;
        double z = leaf.z;
        if (is_in_bbox(x, y, z, xyzxyz)) {
            new_tree.push_back(leaf);
        }
    }
    return new_tree;
}

std::unordered_map<int, Node> get_pos_dict(const std::vector<Node>& tree) {
    std::unordered_map<int, Node> pos_dict;
    for (const auto& leaf : tree) {
        pos_dict[leaf.idx] = leaf;
    }
    return pos_dict;
}

std::unordered_map<int, std::vector<int>> get_child_dict(const std::vector<Node>& tree){
    std::unordered_map<int, std::vector<int>> child_dict;
    for (const auto& leaf : tree) {
        child_dict[leaf.parent].push_back(leaf.idx);
    }
    return child_dict;
}

std::vector<Marker> get_furcation(const std::vector<Node>& tree, std::string& furfile) {
    std::vector<Marker> furcation;
    std::unordered_map<int, Node> pos_dict = get_pos_dict(tree);
    std::unordered_map<int, std::vector<int>> child_dict = get_child_dict(tree);
    for (const auto& entry : child_dict) {
        int parent = entry.first;
        const std::vector<int>& children = entry.second;
        if (children.size() > 1) {
            Marker marker;
            Node node = pos_dict[parent];
            marker.x = node.x;
            marker.y = node.y;
            marker.z = node.z;
            marker.r = node.r;
            furcation.push_back(marker);
        }
    }
    std::cout << "furcation number is: " << furcation.size() << std::endl;
    bool write = writeMarkerstree(furcation, furfile);
    if (write){std::cout << "done written furcation file." << std::endl;}
    return furcation;
}
    

