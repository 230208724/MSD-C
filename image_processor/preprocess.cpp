#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <opencv2/opencv.hpp> // Make sure to install OpenCV and link it correctly

#include "preprocess.h"

bool equalizationImagefileToU8(std::string& imgfile, float& p, std::string& newfile, cv::Mat& equalout) {
    cv::Mat image = cv::imread(imgfile);
    double imax, imin;
    cv::minMaxLoc(image, &imin, &imax);
    imin = cv::percentile(image, p);
    image = (image - imin) * 255 / (imax - imin);
    cv::threshold(image, image, 0, 0, cv::THRESH_TOZERO);
    image.convertTo(image, CV_8U);
    cv::imwrite(newfile, image);
    equalout = image;
    return true;
}

//Note that the C++ conversion uses the OpenCV library to load and save images, 
//perform min-max calculation and percentile calculation, as well as thresholding and type conversion.
//Also note that this is just one possible conversion, 
bool equalizationImageToU8(cv::Mat& image, float p=0) {
    double imax, imin;
    cv::minMaxLoc(image, &imin, &imax);
    imin = cv::percentile(image, p);
    image = (image - imin) * 255 / (imax - imin);
    cv::threshold(image, image, 0, 0, cv::THRESH_TOZERO);
    // for (auto& val : image) {
    //     val = (val - imin) * 255 / (imax - imin);
    //     val = max(0.0, val);
    // }
    image.convertTo(image, CV_8U);
    return true;
}

bool getMipfileFromImagefile(const std::string& imgfile, int& axis, const std::string& mipfile) {
    cv::Mat image;
    image = cv::imread(imgfile, cv::IMREAD_GRAYSCALE);
    std::cout << "image shape: " << image.size() << std::endl;
    cv::Mat mip;
    cv::reduce(image, mip, axis, cv::REDUCE_MAX);
    std::cout << "MIP shape: " << mip.size() << std::endl;
    cv::imwrite(mipfile, mip);
    return true;
}

bool getMipFromImage(cv::Mat& image, int axis, cv::Mat& mip) {
    std::cout << "image shape: " << image.size() << std::endl;

    cv::reduce(image, mip, axis, cv::REDUCE_MAX);
    std::cout << "MIP shape: " << mip.size() << std::endl;
    
    return true;
}

int main() {
    std::string imgfile = "path/to/image.jpg";
    float p = 10;
    std::string outdir = "path/to/outdir";
    std::string newfile = equalization_u8(imgfile, p, outdir);
    std::cout << "Generated file: " << newfile << std::endl;
    return 0;
}



