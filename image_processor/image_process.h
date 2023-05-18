#ifndef IMAGE_PROCESS_H
#define IMAGE_PROCESS_H

#include <iostream>
#include <string>

#include "/usr/include/x86_64-linux-gnu/tiffio.h"

//#include "opencv2/opencv.hpp"
#include "/usr/include/opencv2/opencv.hpp"

void openU8(const char* file,std::vector<cv::Mat> &buffer,int *size);
void openU16(const char* file, uint16 **buffer,int *size);
void saveTiff(const char *path,uint16 *buffer,int *size);
cv::Mat extractSubImage(const cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end);
bool setSubImage(cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end);
double calculatePercentile(cv::Mat image, double percentile);
void findMinMaxLocValueIterative(const cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, int& minValue, int& maxValue);
void findMinMaxLocValueSliced(const cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, int& minValue, int& maxValue);
bool printImage3DInfo(const cv::Mat& image3D);
bool equalizationImagefileToU8(std::string& imgfile, float& p, std::string& newfile, cv::Mat& equaliout);
bool equalizationImageToU8(cv::Mat& image, float p = 0);
//bool getMipfileFromImagefile(const std::string& imgfile, int& axis, const std::string& newfile);
//bool getMipFromImage(cv::Mat& image, int axis, cv::Mat& mip);

#endif  // IMAGE_PROCESS_H
