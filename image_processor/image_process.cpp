#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "/usr/include/x86_64-linux-gnu/tiffio.h"
//#include "opencv2/opencv.hpp" // Make sure to install OpenCV and link it correctly
#include "/usr/include/opencv2/opencv.hpp"

#include "image_process.h"

// Function to print 3D image information
bool printImage3DInfo(const cv::Mat& image3D) {
    // Check if image is empty
    if (image3D.empty()) {
        std::cout << "unable to read image" << std::endl;
        return false;
    }

    // Get image dimensions
    int numSlices = image3D.size[0];
    int height = image3D.size[1];
    int width = image3D.size[2];
    int numChannels = image3D.channels();
    int dataType = image3D.type();
    int pixelSize = image3D.elemSize();

    // Print image information
    std::cout << "slice: " << numSlices << std::endl;
    std::cout << "height: " << height << " voxe;" << std::endl;
    std::cout << "width: " << width << " voxel" << std::endl;
    std::cout << "channel number: " << numChannels << std::endl;
    std::cout << "voxeldatatype: " << dataType << std::endl;
    std::cout << "voxelpixelsize: " << pixelSize << std::endl;
    std::cout << "imagepixelsize: " << image3D.total() * pixelSize << std::endl;
    return true;
}

cv::Mat openU16(const std::string& file)
{
    //-------------使用TIFFOpen函数以只读形式打开图像。
	TIFF *tif = TIFFOpen(file.c_str(), "r");      
	if (tif == nullptr)
	{
		std::cout << "unable to read image";
        cv::Mat Nimage;
		return Nimage;
	}
    //-------------获取单帧图像的长高
	int nTotalFrame = TIFFNumberOfDirectories(tif);
	int width, height;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	uint16_t datatype;
	TIFFGetField(tif,TIFFTAG_SAMPLEFORMAT,&datatype);
    size_t stripSize = TIFFStripSize(tif);
    //TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &channel);
    
    //---------------创建cv2 mat
	int size[3] = {nTotalFrame,height,width};
    cv::Mat image3D(3,size,CV_16UC1,cv::Scalar(0));
    uint16 maxValue = std::numeric_limits<uint16>::lowest();
    cv::Point3d maxLoc(0,0,0);    
    for (int z = 0; z < nTotalFrame; z++) {
        // 申请单帧内存，此处注意16位图像在此对应的是uint16，申请内存空间大小按TIFFStripSize的返回值进行申请
        uint16* sliceBuffer = new uint16[stripSize];
       
        for (int y = 0; y < height; y++)
        {
            // 每次都是从当前slice第一个字节开始写入
            TIFFReadScanline(tif, (&(sliceBuffer)[0] + y * int(width)), y); //---按行读取
        }
       
        int ss = 0;
        for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				uint16 value =  sliceBuffer[ss];                          
				//std::cout << value << std::endl;
                image3D.at<int>(z,y,x) = value;
                if (maxValue<value) 
                    {
                        //std::cout << maxValue << std::endl;
                        maxValue = value;
                        maxLoc.z = z;
                        maxLoc.y = y;
                        maxLoc.x = x;
                    }
                ss += 1;
            }
        }
        //std::cout << "for z" << z << " : ss = " << ss << std::endl;
        TIFFReadDirectory(tif);
    }

	TIFFClose(tif);

	// info print
    // Function to print 3D image information
    // Get image dimensions
    int numSlices = image3D.size[0];
    int imageheight = image3D.size[1];
    int imagewidth = image3D.size[2];
    int numChannels = image3D.channels();
    int dataType = image3D.type();
    int pixelSize = image3D.elemSize();

    // Print image information
    std::cout << "元数据类型：" << datatype << std::endl;
    std::cout << "切片数: " << numSlices << std::endl;
    std::cout << "图像高度: " << imageheight << " 像素" << std::endl;
    std::cout << "图像宽度: " << imagewidth << " 像素" << std::endl;
    std::cout << "图像通道数: " << numChannels << std::endl;
    std::cout << "图像数据类型: " << dataType << std::endl;
    std::cout << "图像最大像素值: " << maxValue << std::endl;
    std::cout << "图像最大像素位置: " << maxLoc.z << " " << maxLoc.y << " " << maxLoc.x << std::endl;
    std::cout << "图像每个像素字节数: " << pixelSize << std::endl;
    std::cout << "图像总字节数: " << image3D.total() * pixelSize << std::endl;
    return image3D;
}

void saveTiff(const char *path,uint16 *buffer,int *size)
{
	int width = size[0];
	int height = size[1];
	int slice = size[2];
 
	TIFF* out = TIFFOpen(path, "w");
	if (out)
	{
		int N_size = 0;
		size_t nCur = 0;
		//UChar den = (sizeof(T) == 1) ? 1 : 4;
		do{
			TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
			TIFFSetField(out, TIFFTAG_PAGENUMBER, slice);
			TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32)width);
			TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32)height);
			//TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, 2);
			/*TIFFSetField(out, TIFFTAG_YRESOLUTION, 196.0f);
			TIFFSetField(out, TIFFTAG_XRESOLUTION, 204.0f);*/
			TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
			// 
			TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 16);    //根据图像位深填不同的值
			TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
			TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
			TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
			TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
			TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
 
 
			for (int m = 0; m < height; m++)
			{
				TIFFWriteScanline(out, &buffer[N_size + width*m], m, 0);
			}
			//TIFFWriteEncodedStrip(out, 0, &buffer[N_size], width * height);      //另一种写入方法
 
			++nCur;
			N_size = N_size + width*height;
		} while (TIFFWriteDirectory(out) && nCur < slice);
		TIFFClose(out);
 
		std::cout << "save over";
	}
}

cv::Mat extractSubImage(const cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end)
{
    // 确保起点和终点在图像范围内
    int x_min = std::max(0, std::min(x_start, x_end));
    int y_min = std::max(0, std::min(y_start, y_end));
    int z_min = std::max(0, std::min(z_start, z_end));
    int x_max = std::min(image.cols - 1, std::max(x_start, x_end));
    int y_max = std::min(image.rows - 1, std::max(y_start, y_end));
    int z_max = std::min(image.size[2] - 1, std::max(z_start, z_end));

    // 计算子图像的宽度、高度和深度
    int width = x_max - x_min + 1;
    int height = y_max - y_min + 1;
    int depth = z_max - z_min + 1;

    // 创建子图像
    cv::Mat subImage(height, width, image.type());

    // 复制子图像的像素值
    for (int z = z_min; z <= z_max; ++z)
    {
        for (int y = y_min; y <= y_max; ++y)
        {
            for (int x = x_min; x <= x_max; ++x)
            {
                subImage.at<uchar>(y - y_min, x - x_min) = image.at<uchar>(y, x, z);
            }
        }
    }

    return subImage;
}

void setPixelValue(cv::Mat& image, int x, int y, int z, uchar value)
{
    if (image.empty())
    {
        // 图像为空，无法设置像素值
        return;
    }

    if (x >= 0 && x < image.cols && y >= 0 && y < image.rows && z >= 0 && z < image.size[2])
    {
        // 修改像素值
        image.at<uchar>(y, x, z) = value;
    }
    else
    {
        // 像素坐标超出图像范围
        // 可以选择抛出异常或输出错误信息
        // 在这个示例中，我们选择输出错误信息
        std::cout << "candidate out of image range!" << std::endl;
    }
}

bool setSubImage(cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end)
{
    // 确保起点和终点在图像范围内
    int x_min = std::max(0, std::min(x_start, x_end));
    int y_min = std::max(0, std::min(y_start, y_end));
    int z_min = std::max(0, std::min(z_start, z_end));
    int x_max = std::min(image.cols - 1, std::max(x_start, x_end));
    int y_max = std::min(image.rows - 1, std::max(y_start, y_end));
    int z_max = std::min(image.size[2] - 1, std::max(z_start, z_end));

    // 计算子图像的宽度、高度和深度
    int width = x_max - x_min + 1;
    int height = y_max - y_min + 1;
    int depth = z_max - z_min + 1;

    // 更改子图像的像素值为0
    for (int z = z_min; z <= z_max; ++z)
    {
        for (int y = y_min; y <= y_max; ++y)
        {
            for (int x = x_min; x <= x_max; ++x)
            {
                setPixelValue(image, x, y, z, 0);
            }
        }
    }
    return true;
}

double calculatePercentile(cv::Mat image, double percentile) {
    // Flatten the image to a 1D array
    cv::Mat flattenedImage = image.reshape(1, 1);

    // Sort the flattened array
    cv::Mat sortedImage;
    cv::sort(flattenedImage, sortedImage, cv::SORT_ASCENDING);

    // Compute the index corresponding to the percentile
    double index = percentile / 100.0 * (sortedImage.cols - 1);
    std::cout << "image percentile is " << std::to_string(index) << "out of image shape" << std::to_string(sortedImage.cols) << std::endl;

    // Interpolate the pixel value at the computed index
    int lowerIndex = static_cast<int>(std::floor(index));
    int upperIndex = static_cast<int>(std::ceil(index));
    double interpolatedValue = sortedImage.at<uchar>(0, lowerIndex) +
                               (index - lowerIndex) * (sortedImage.at<uchar>(0, upperIndex) - sortedImage.at<uchar>(0, lowerIndex));
    std::cout << "image percentile number is " << std::to_string(interpolatedValue) <<std::endl;

    return interpolatedValue;
}

// Function to find the minimum and maximum value positions and values in a 3D image (Z/Y/X)
void findMinMaxLocValueIterative(const cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, int& minValue, int& maxValue) {
    minValue = std::numeric_limits<int>::max();
    maxValue = std::numeric_limits<int>::lowest();

    // Iterate over the 3D image
    for (int z = 0; z < image3D.size[0]; ++z) {
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                double pixelValue = image3D.at<double>(z, y, x);

                // Update minLoc and minValue if current value is smaller
                if (pixelValue < minValue) {
                    minValue = pixelValue;
                    minLoc = cv::Point3i(x, y, z);
                }

                // Update maxLoc and maxValue if current value is larger
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    maxLoc = cv::Point3i(x, y, z);
                }
            }
        }
    }

    std::cout << "image min/max is " << std::to_string(minValue) << "/" << std::to_string(maxValue) << std::endl;
}

// Function to find the minimum and maximum value positions and values in a 3D image (Z/Y/X)
void findMinMaxLocValueSliced(const cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, int& minValue, int& maxValue) {
    minValue = std::numeric_limits<int>::max();
    maxValue = std::numeric_limits<int>::lowest();

    // Iterate over each 2D slice of the 3D image
    for (int z = 0; z < image3D.size[0]; ++z) {
        cv::Mat slice = image3D.row(z).reshape(1, image3D.size[1]);  // Extract a 2D slice from the 3D image

        double sliceMinValue, sliceMaxValue;
        cv::Point minSliceLoc, maxSliceLoc;

        // Find the minimum and maximum value positions in the 2D slice
        cv::minMaxLoc(slice, &sliceMinValue, &sliceMaxValue, &minSliceLoc, &maxSliceLoc);

        // Update minLoc and minValue if the minimum value of the slice is smaller
        if (sliceMinValue < minValue) {
            minValue = sliceMinValue;
            minLoc = cv::Point3i(minSliceLoc.x, minSliceLoc.y, z);
        }

        // Update maxLoc and maxValue if the maximum value of the slice is larger
        if (sliceMaxValue > maxValue) {
            maxValue = sliceMaxValue;
            maxLoc = cv::Point3i(maxSliceLoc.x, maxSliceLoc.y, z);
        }
    }

    std::cout << "image min/max is " << std::to_string(minValue) << "/" << std::to_string(maxValue) << std::endl;
}

bool equalizationImagefileToU8(std::string& imgfile, float& p, std::string& newfile, cv::Mat& equalout) {
    //cv::Mat image = cv::imread(imgfile,cv::IMREAD_ANYDEPTH | cv::IMREAD_ANYCOLOR | cv::IMREAD_UNCHANGED);
    cv::Mat image = openU16(imgfile);
    //bool info = printImage3DInfo(image);
    bool info = image.empty();
    if (info==false){return false;}

    return true;
    int imax, imin;
    cv::Point3i minLoc, maxLoc;
    findMinMaxLocValueIterative(image, minLoc, maxLoc, imin, imax);
    imin = calculatePercentile(image, p);
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
bool equalizationImageToU8(cv::Mat& image, float p) {
    double imax, imin;
    cv::minMaxLoc(image, &imin, &imax);
    imin = calculatePercentile(image, p);
    image = (image - imin) * 255 / (imax - imin);
    cv::threshold(image, image, 0, 0, cv::THRESH_TOZERO);
    // for (auto& val : image) {
    //     val = (val - imin) * 255 / (imax - imin);
    //     val = max(0.0, val);
    // }
    image.convertTo(image, CV_8U);
    return true;
}

/* unsettled
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

    cv::reduce(image, mip, axis, [](const cv::Mat& a, const cv::Mat& b){return cv::max(a, b);});
    std::cout << "MIP shape: " << mip.size() << std::endl;
    
    return true;
}
*/




