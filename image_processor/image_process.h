#ifndef IMAGE_PROCESS_H
#define IMAGE_PROCESS_H

#include <iostream>
#include <string>

#include "/usr/include/x86_64-linux-gnu/tiffio.h"

//#include "opencv2/opencv.hpp"
#include "/usr/include/opencv2/opencv.hpp"

uint32_t countValueNumberU8(cv::Mat& image3D, uint8_t value);
int SumMatU8(cv::Mat& image3D);
void copyMat16(cv::Mat& newimage, cv::Mat& image3D);
void copyMat8(cv::Mat& newimage, cv::Mat& image3D);
void saveU8(std::string& imagefile, cv::Mat& image3D, bool& save);
void infoPrint(cv::Mat& image);
void findMinMaxLocValueIterative(cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, uint16_t& minValue, uint16_t& maxValue);
void findMinMaxLocValueIterativeU8(cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, uint8_t& minValue, uint8_t& maxValue);
void findMinMaxLocValueIterative2D(cv::Mat& image2D, cv::Point2i& minLoc, cv::Point2i& maxLoc, uint16_t& minValue, uint16_t& maxValue);
void findMinMaxLocValueByZU8(cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues);
void findMinMaxLocValueByYU8(cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues);
void findMinMaxLocValueByXU8(cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues);
void flatten(cv::Mat& image3D, cv::Mat& flattenedImage);
void calculatePercentile(cv::Mat& image, double& percentile, uint16_t& interpolatedValue, bool& per);
void equalizationImageToU8(cv::Mat& imagein, cv::Mat& imageout, double& pin, bool& equal);
cv::Mat extractSubImage(cv::Mat& image, cv::Mat& subImage, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end);
void setPixelValue(cv::Mat& image, int x, int y, int z, uint8_t value=0);
void setSubImage(cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end, uint8_t value=0);

#endif  // IMAGE_PROCESS_H
