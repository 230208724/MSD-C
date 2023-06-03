#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "/usr/include/x86_64-linux-gnu/tiffio.h"
#include "/usr/include/opencv2/opencv.hpp"

#include "image_process.h"

uint32_t countValueNumberU8(cv::Mat& image3D, uint8_t value){
    uint32_t countNumber = 0;
    int slice = image3D.size[0];
	int height = image3D.size[1];
	int width = image3D.size[2];
    for (int z = 0; z < slice; z++){
        for (int y = 0; y < height; y++){
            for (int x = 0; x < width; x++){
                if (image3D.ptr<uint8_t>(z,y)[x] > value){
                    countNumber += 1;
                }
            }
        }
    }
    return countNumber;
}

int SumMatU8(cv::Mat& image3D){
    uint8_t sum = 0;
    int slice = image3D.size[0];
	int height = image3D.size[1];
	int width = image3D.size[2];
    for (int z = 0; z < slice; z++){
        for (int y = 0; y < height; y++){
            for (int x = 0; x < width; x++){
                sum += image3D.ptr<uint8_t>(z,y)[x];
            }
        }
    }
    return static_cast<int>(sum);
}

/*------------------uint16 mat copy----------------*/
void copyMat16(cv::Mat& newimage, cv::Mat& image3D){
    int slice = image3D.size[0];
	int height = image3D.size[1];
	int width = image3D.size[2];
    for (int z = 0; z < slice; z++){
        for (int y = 0; y < height; y++){
            for (int x = 0; x < width; x++){
                newimage.ptr<uint16_t>(z,y)[x] = image3D.ptr<uint16_t>(z,y)[x];
            }
        }
    }
}
/*---------*/


/*---------------uint8 mat copy-------------------*/
void copyMat8(cv::Mat& newimage, cv::Mat& image3D){
    int slice = image3D.size[0];
	int height = image3D.size[1];
	int width = image3D.size[2];
    for (int z = 0; z < slice; z++){
        for (int y = 0; y < height; y++){
            for (int x = 0; x < width; x++){
                newimage.ptr<uint8_t>(z,y)[x] = static_cast<uint8_t>(image3D.ptr<uint16_t>(z,y)[x]);
            }
        }
    }
}
/*-----------*/


/*------------save cv mat 3d with tiff-----------------*/
void saveU8(std::string& imagefile, cv::Mat& image3D, bool& save)
{
	int slice = image3D.size[0];
	int height = image3D.size[1];
	int width = image3D.size[2];
    uint8_t *buffer = new uint8_t[slice*height*width];
    uint32_t index;
    uint8_t value;
    for (int z = 0; z < slice; z++){
        for (int y = 0; y < height; y++){
            for (int x = 0; x < width; x++){
                index = z * (height*width) + y * width + x;
                value = image3D.ptr<uint8_t>(z,y)[x];
                buffer[index] = value;
                }
        }
    }

	TIFF* out = TIFFOpen(imagefile.c_str(), "w");
    if (out)
	{
        int n;
		int N_size = 0;
		int nCur = 0;
		do{
			TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
			TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
            TIFFSetField(out,TIFFTAG_RESOLUTIONUNIT,2);
            TIFFSetField(out,TIFFTAG_PLANARCONFIG,1);
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(out,TIFFTAG_BITSPERSAMPLE,8);
            TIFFSetField(out,TIFFTAG_SAMPLESPERPIXEL,1);
            TIFFSetField(out,TIFFTAG_IMAGEWIDTH,width);
            TIFFSetField(out,TIFFTAG_IMAGELENGTH,height);

            n = TIFFWriteEncodedStrip(out, 0, &buffer[N_size], width * height);  
            if (n < 0)
                std::cout << "error writing TIFF file" << std::endl;

            TIFFWriteDirectory(out);
            
            //TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
			//TIFFSetField(out, TIFFTAG_PAGENUMBER, slice);
			//TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
			// TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
			// for (int m = 0; m < height; m++)
			// {
            //     TIFFWriteScanline(out, &buffer[N_size + width*m], m , (uint16)0U);
			// }
			
            //TIFFSetDirectory(out, nCur);

            N_size = N_size + width*height;
            ++nCur;
            
		} while (nCur < slice);
		
        TIFFClose(out);
		//std::cout << "save over" << std::endl;
	}
    save = true;
}
/*-----*/



/*------------------Function to find the minimum and maximum value positions and values in a 3D image (Z/Y/X)--------------------*/
void findMinMaxLocValueIterative( cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, uint16_t& minValue, uint16_t& maxValue) {
    minValue = std::numeric_limits<uint16_t>::max();
    maxValue = std::numeric_limits<uint16_t>::lowest();

    // Iterate over the 3D image
    for (int z = 0; z < image3D.size[0]; ++z) {
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                uint16_t pixelValue = image3D.ptr<uint16_t>(z, y)[x];

                // Update minLoc and minValue if current value is smaller
                if (pixelValue < minValue) {
                    minValue = pixelValue;
                    minLoc.z = z;
                    minLoc.y = y;
                    minLoc.x = x;
                }

                // Update maxLoc and maxValue if current value is larger
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    maxLoc.z = z;
                    maxLoc.y = y;
                    maxLoc.x = x;
                }
            }
        }
    }
    std::cout << "U16 3D min/max and location: " << minValue << "/" << maxValue << ", ";
    std::cout << minLoc.z << " " << minLoc.y << " " << minLoc.x;
    std::cout << "/" << maxLoc.z << " " << maxLoc.y << " " << maxLoc.x << std::endl;
}

void findMinMaxLocValueIterativeU8( cv::Mat& image3D, cv::Point3i& minLoc, cv::Point3i& maxLoc, uint8_t& minValue, uint8_t& maxValue) {
    minValue = std::numeric_limits<uint8_t>::max();
    maxValue = std::numeric_limits<uint8_t>::lowest();

    // Iterate over the 3D image
    for (int z = 0; z < image3D.size[0]; ++z) {
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                uint8_t pixelValue = image3D.ptr<uint8_t>(z, y)[x];

                // Update minLoc and minValue if current value is smaller
                if (pixelValue < minValue) {
                    minValue = pixelValue;
                    minLoc.z = z;
                    minLoc.y = y;
                    minLoc.x = x;
                }

                // Update maxLoc and maxValue if current value is larger
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    maxLoc.z = z;
                    maxLoc.y = y;
                    maxLoc.x = x;
                }
            }
        }
    }
    std::cout << "U8 3D min/max and location: " << std::to_string(minValue) << "/" << std::to_string(maxValue) << ", ";
    std::cout << minLoc.z << " " << minLoc.y << " " << minLoc.x;
    std::cout << "/" << maxLoc.z << " " << maxLoc.y << " " << maxLoc.x << std::endl;
}

void findMinMaxLocValueIterative2D( cv::Mat& image2D, cv::Point2i& minLoc, cv::Point2i& maxLoc, uint16_t& minValue, uint16_t& maxValue) {
    minValue = std::numeric_limits<uint16_t>::max();
    maxValue = std::numeric_limits<uint16_t>::lowest();

    // Iterate over the 2D image
    for (int y = 0; y < image2D.size[0]; ++y) {
        for (int x = 0; x < image2D.size[1]; ++x) {
            uint16_t pixelValue = image2D.ptr<uint16_t>(y)[x];

            // Update minLoc and minValue if current value is smaller
            if (pixelValue < minValue) {
                minValue = pixelValue;
                minLoc.y = y;
                minLoc.x = x;
            }

            // Update maxLoc and maxValue if current value is larger
            if (pixelValue > maxValue) {
                maxValue = pixelValue;
                maxLoc.y = y;
                maxLoc.x = x;
            }
        }
    }
    std::cout << "U16 2D min/max and location: " << minValue << "/" << maxValue << ", ";
    std::cout << minLoc.y << " " << minLoc.x;
    std::cout << "/" << maxLoc.y << " " << maxLoc.x << std::endl;
}


void findMinMaxLocValueByZU8( cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues) {
    // Iterate over the 3D image
    for (int z = 0; z < image3D.size[0]; ++z) {
        uint8_t minValue = std::numeric_limits<uint8_t>::max();
        uint8_t maxValue = std::numeric_limits<uint8_t>::lowest();
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                uint8_t pixelValue = image3D.ptr<uint8_t>(z, y)[x];

                if (pixelValue < minValue) {
                    minValue = pixelValue;
                }
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    }
            }
        }
        minValues[z] = minValue;
        maxValues[z] = maxValue;
    }
    for (int zz = 0; zz < image3D.size[0] ; zz++)
    {
        std::cout << std::to_string(minValues[zz]) << " " << std::to_string(maxValues[zz]) << ", ";
    }
}
void findMinMaxLocValueByYU8( cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues) {
    // Iterate over the 3D image
    for (int y = 0; y < image3D.size[1]; ++y) {
        uint8_t minValue = std::numeric_limits<uint8_t>::max();
        uint8_t maxValue = std::numeric_limits<uint8_t>::lowest();
        for (int z = 0; z < image3D.size[0]; ++z) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                uint8_t pixelValue = image3D.ptr<uint8_t>(z, y)[x];

                if (pixelValue < minValue) {
                    minValue = pixelValue;
                }
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    }
            }
        }
        minValues[y] = minValue;
        maxValues[y] = maxValue;
    }
    for (int yy = 0; yy < image3D.size[1] ; yy++)
    {
        std::cout << std::to_string(minValues[yy]) << " " << std::to_string(maxValues[yy]) << ", ";
    }
}
void findMinMaxLocValueByXU8( cv::Mat& image3D, uint8_t* minValues, uint8_t* maxValues) {
    // Iterate over the 3D image
    for (int x = 0; x < image3D.size[2]; ++x) {
        uint8_t minValue = std::numeric_limits<uint8_t>::max();
        uint8_t maxValue = std::numeric_limits<uint8_t>::lowest();
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int z = 0; z < image3D.size[0]; ++z) {    
                uint8_t pixelValue = image3D.ptr<uint8_t>(z, y)[x];

                if (pixelValue < minValue) {
                    minValue = pixelValue;
                }
                if (pixelValue > maxValue) {
                    maxValue = pixelValue;
                    }
            }
        }
        minValues[x] = minValue;
        maxValues[x] = maxValue;
    }
    for (int xx = 0; xx < image3D.size[2] ; xx++)
    {
        std::cout << std::to_string(minValues[xx]) << " " << std::to_string(maxValues[xx]) << ", ";
    }
}
/*---*/


/*-----------info print-------*/
void infoPrint(cv::Mat& image){
    int numSlices = image.size[0];
    int imageheight = image.size[1];
    int imagewidth = image.size[2];
    int numChannels = image.channels();
    int dataType = image.type();
    int pixelSize = image.elemSize();
    int imagePixelSize = image.total() * pixelSize;
    std::cout << "mat size(slice height width)/channel number/voxel datatype/voxelpixelsize/imagevoxelsize: ";
    std::cout << numSlices << "voxel " << imageheight << "voxel " << imagewidth << "voxel";
    std::cout << "/" << numChannels << "/" << dataType << "/" << pixelSize << "/" << imagePixelSize << std::endl;
    cv::Point3i minLoc, maxLoc;  
    if (dataType==0){
        uint8_t minValue = std::numeric_limits<uint8_t>::max();
        uint8_t maxValue = std::numeric_limits<uint8_t>::lowest();
        findMinMaxLocValueIterativeU8(image,minLoc,maxLoc,minValue,maxValue);    
    }
    else{
        uint16_t minValue = std::numeric_limits<uint16_t>::max();
        uint16_t maxValue = std::numeric_limits<uint16_t>::lowest();
        findMinMaxLocValueIterative(image,minLoc,maxLoc,minValue,maxValue);
    }
}
/*------*/



/*----------flatten image size to ------------------
---1,image.size[0]*image.size[1]*image.size[2]-----*/
void flatten(cv::Mat& image3D, cv::Mat& flattenedImage){
    std::cout << "start flatten" << std::endl;
    int slice = image3D.size[1] * image3D.size[2];
    int index;
    uint16 value = 0;
    for (int z = 0; z < image3D.size[0]; ++z) {
        for (int y = 0; y < image3D.size[1]; ++y) {
            for (int x = 0; x < image3D.size[2]; ++x) {
                index = z * (slice) + y * image3D.size[2] + x;
                flattenedImage.ptr<uint16_t>(0)[index] = image3D.ptr<uint16_t>(z, y)[x];
            }
        }
    }
    std::cout << "done flatten to image of size ";
    std::cout << flattenedImage.size[0] << "/" << flattenedImage.size[1] << std::endl;
}
/*-----------------------------------------------------------------------------------*/



/*---------------------------------------calculate percentile ----------------------------------*/
void calculatePercentile(cv::Mat& image, double& percentile, uint16_t& interpolatedValue, bool& per) {
    std::cout << "start percentile" << std::endl;

    // Flatten the image to a 1D array int cn, note reshape declaration: Mat reshape(int cn, const int newndims, const int* newsz);
    int imagesize[2] = {1,image.size[0]*image.size[1]*image.size[2]};
    cv::Mat flattenedImage(2,imagesize,CV_16UC1,cv::Scalar(0));
    flatten(image,flattenedImage);
    cv::Point2i _,__;
    uint16_t ___,____;
    findMinMaxLocValueIterative2D(flattenedImage,_,__,___,____);

    // Sort the flattened array and Compute the index corresponding to the percentile
    std::cout << "start sort" << std::endl;
    cv::Mat sortedImage(2,imagesize,CV_16UC1,cv::Scalar(0));
    cv::sort(flattenedImage, sortedImage, cv::SORT_ASCENDING);
    findMinMaxLocValueIterative2D(sortedImage,_,__,___,____);
    double index = percentile / 100.0 * (sortedImage.size[1]) - 1;
    int lowerIndex = static_cast<int>(std::floor(index));
    std::cout << "image percentile is id " << lowerIndex << " out of " << sortedImage.size[1] - 1 << std::endl;
    interpolatedValue = sortedImage.ptr<uint16_t>(0)[lowerIndex];
    std::cout << "image percentile number is " << interpolatedValue <<std::endl;
    std::cout << "done percentile" << std::endl;
    per = true;
}
/*--------------*/


/*----------------------output equlization image and save image to outfile-------------------*/
void equalizationImageToU8(cv::Mat& imagein, cv::Mat& imageout, double& pin, bool& equal) {
    std::cout << "start equalization" << std::endl;

    //find the percentile
    uint16_t imin;
    bool per;
    calculatePercentile(imagein, pin, imin, per);

    //find the minimum and maximum value positions and values in a 3D image (Z/Y/X)
    uint16_t imax, imin_;
    cv::Point3i minLoc, maxLoc;
    findMinMaxLocValueIterative(imagein, minLoc, maxLoc, imin_, imax);   

    //equealization
    imagein = (imagein - imin) * 255 / (imax - imin);
    //cv::threshold(imagein, imagein, 0, 0, cv::THRESH_TOZERO);
    copyMat8(imageout,imagein);
    equal = true;    
    std::cout << "done equlization" << std::endl;
}
/*------------------*/



/*---------------------------------extract and return sub image within certain block-------------------------------------*/
cv::Mat extractSubImage(cv::Mat& image, cv::Mat& subImage, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end)
{   
    for (int z = z_start; z <= z_end; ++z)
    {
        for (int y = y_start; y <= y_end; ++y)
        {
            for (int x = x_start; x <= x_end; ++x)
            {
                subImage.ptr<uint8_t>(z - z_start, y - y_start)[x - x_start] = image.ptr<uint8_t>(z,y)[x];
            }
        }
    }
    return subImage;
}
/*-----------------*/


/*-----------------set certain block to certain value-----------------*/
void setPixelValue(cv::Mat& image, int x, int y, int z, uint8_t value)
{
    if (x >= 0 && x < image.size[2] && y >= 0 && y < image.size[1] && z >= 0 && z < image.size[0])
    {
        // change value
        image.ptr<uint8_t>(z, y)[x] = value;
    }
    else
    {
        std::cout << "candidate z "<< std::to_string(z) << " y " << std::to_string(y) << " x " << std::to_string(x);
        std::cout << "out of image range z/y/x: " << image.size[0] << "/" << image.size[1] << "/" << image.size[2];
    }
}

void setSubImage(cv::Mat& image, int x_start, int y_start, int z_start, int x_end, int y_end, int z_end, uint8_t value)
{
    // set start and end point
    int x_min = std::max(0, std::min(x_start, x_end));
    int y_min = std::max(0, std::min(y_start, y_end));
    int z_min = std::max(0, std::min(z_start, z_end));
    int x_max = std::min(image.size[2] - 1, std::max(x_start, x_end));
    int y_max = std::min(image.size[1] - 1, std::max(y_start, y_end));
    int z_max = std::min(image.size[0] - 1, std::max(z_start, z_end));

    // reassign value to 0
    for (int z = z_min; z <= z_max; ++z)
    {
        for (int y = y_min; y <= y_max; ++y)
        {
            for (int x = x_min; x <= x_max; ++x)
            {
                setPixelValue(image, x, y, z, value);
            }
        }
    }
}
/*---------*/