#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>

bool equalizationImagefileToU8(std::string& imgfile, float& p, std::string& newfile, cv::Mat& equaliout);
bool equalizationImageToU8(cv::Mat& image, float p = 0);
bool getMipfileFromImagefile(const std::string& imgfile, int& axis, const std::string& newfile);
bool getMipFromImage(cv::Mat& image, int axis, cv::Mat& mip);

#endif  // PREPROCESS_H
