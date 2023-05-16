#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <unordered_map>

#include "../config.h"

std::vector<Node> parseSwcfile(const std::string& swc_file)
void writeNodestree(const std::vector<Node>& tree, const std::string& swc_file, const std::vector<std::string>& header = {}) 
std::vector<Marker> parseMarkerfile(const std::string& marker_file)
void writeMarkerstree(const std::vector<Marker>& markertree, const std::string& markerfile, const std::vector<std::string>& header = {}) 
Node find_first_soma(const std::vector<Node>& tree, int p_soma = -1) 
bool is_in_box(double x, double y, double z, const std::vector<int>& xyzshape) 
bool is_in_bbox(double x, double y, double z, const std::vector<std::vector<int>>& xyzxyz)
std::vector<Node> crop_bbox_tree(const std::vector<Node>& tree, const std::vector<std::vector<int>>& xyzxyz)
std::unordered_map<int, Node> get_pos_dict(const std::vector<Node>& tree) 
std::unordered_map<int, std::vector<int>> get_child_dict(const std::vector<Node>& tree)
std::vector<std::vector<double>> get_furcation(const std::vector<Node>& tree) 