#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <array>
#include <fstream>
#include "/usr/include/x86_64-linux-gnu/tiffio.h"
#include "/usr/include/opencv2/opencv.hpp" // Make sure to install OpenCV and link it correctly

#include "config.h"

#include "image_processor/image_process.h"

#include "image_processor/reconstruction.h"

#include "swc_processor/swc_basic_process.h"


/*NOTE that marker.y is filipped 
    because vaa3d pluginin read uint8 with y flipped
    so that true node.y will not hava value*/

template<typename NodeType>
void flipY(NodeType& node, int ySize){
    node.y = ySize-node.y; 
}    

template<typename NodeType>
void flipY(std::vector<NodeType>& nodeTree, int ySize){
    for (int in = 0; in < nodeTree.size(); in++){
        nodeTree[in].y = ySize-nodeTree[in].y;
    } 
} 
/*--*/   



template<typename NodeType>
void Mask(cv::Mat& image, NodeType& node, int& s, bool& mask) {
    // Related values set to 0
    int x = static_cast<int>(std::ceil(node.x));
    int y = static_cast<int>(std::ceil(node.y));
    int z = static_cast<int>(std::ceil(node.z));
    //std::cout << "mask z/y/x " << z << "/" << y << "/" << x << " with s " << s << std::endl;
    int x_start = x - s;
    int y_start = y - s;
    int z_start = z - s;
    int x_end = x + s;
    int y_end = y + s;
    int z_end = z + s;
    uint8_t value = 0;
    setSubImage(image, x_start, y_start, z_start, x_end, y_end, z_end, value);
    mask = true;
}

void Mask(cv::Mat& image, std::vector<Node>& tree, int& s, bool& mask){   
    bool nodecheck;
    for (Node node: tree){
        Mask(image,node,s,nodecheck);
    } 
    mask = true;   
} 



int calc_mass(int x, int y, int z, int r, cv::Mat& image) {
    //note this shape is z/y/x
    int zz = image.size[0];
    int yy = image.size[1];
    int xx = image.size[2];
    
    int x_start = std::max(0, x - r);
    int y_start = std::max(0, y - r);
    int z_start = std::max(0, z - r);
    int x_end = std::min(x + r, xx);
    int y_end = std::min(y + r, yy);
    int z_end = std::min(z + r, zz);
    
    //std:: cout << "subimage z/y/x start: " << z_end << "/" << y_end << "/" << x_end << std::endl;
    //std:: cout << "subimage z/y/x end: " << z_end << "/" << y_end << "/" << x_end << std::endl;
    int width = x_end - x_start + 1;
    int height = y_end - y_start + 1;
    int depth = z_end - z_start + 1;
    int newsize[3] = {depth,height,width};
    
    cv::Mat subImage(3,newsize,CV_8UC1,cv::Scalar(0));
    extractSubImage(image, subImage, x_start, y_start, z_start, x_end, y_end, z_end);
    
    int mass = SumMatU8(subImage);
    return mass;
}



template<typename NodeType>
double calc_distance(NodeType& marker1,NodeType& marker2){
    double dx = std::abs(marker1.x - marker2.x);
    double dy = std::abs(marker1.y - marker2.y);
    double dz = std::abs(marker1.z - marker2.z);
    double ds = std::sqrt(std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2));
    return ds;
}   

template<typename NodeType>
bool find_nearest(std::vector<NodeType>& markers, NodeType& marker, NodeType& nearest){
    double mds = __INTMAX__;
    for (NodeType m : markers){
        double ds = calc_distance(m,marker);
        if (mds>ds){
            nearest.x = m.x;
            nearest.y = m.y;
            nearest.z = m.z;
            mds = ds;
            }
    }   
    return true;
}


bool get_zone_mass_center(cv::Mat &image, 
                            std::vector<Marker> &furcations,
                            std::string& centerfile, 
                            Marker& center, 
                            Marker& nearest) {
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double m_ = 0.0;

    for (Marker marker : furcations) {
        int x = static_cast<int>(marker.x);
        int y = static_cast<int>(marker.y);
        int z = static_cast<int>(marker.z);
        int r = static_cast<int>(marker.r);

        //std::cout << "x/y/z/r: " << x << "/" << y << "/" << z << "/" << r << std::endl;
        int mass = calc_mass(x, y, z, r, image);
        //std::cout << "mass: " << mass << std::endl;
        x_ += x * mass;
        y_ += y * mass;
        z_ += z * mass;
        m_ += mass;
        //std::cout << "sum mass: " << m_ << std::endl;
    }
    
    if (m_==0){
        std::cout << "Done write empty centerfile: " << centerfile << std::endl;
        std::vector<Marker> markers;
        writeMarkerstree(markers, centerfile);
        return false;
        }  
    else{
        Marker center;
        center.x = x_/ m_; 
        center.y = y_ / m_; 
        center.z = z_ / m_; 
        std::cout << "center z/y/x is : " << center.z << "/" << center.y << "/" << center.x << std::endl;
    
        bool near = find_nearest(furcations, center, nearest);
        std::cout << "nearest z/y/x is : " << nearest.z << "/" << nearest.y << "/" << nearest.x << std::endl;
     
        std::cout << "Done write centerfile: " << centerfile << std::endl;
        writeMarker(nearest, centerfile);
        return true;
    }   
}

//int main(int argc, char* argv[]) {
int main() {
    int argc = 3;
    std::string* argv = new std::string[argc];

    argv[0] = "msd";
    argv[1] = "-i";
    argv[2] = "test/112200_046080_048450.tif";

    /*----------------parameter for input decoding----------------------*/    
    if (argc < 3) {
        std::cout << "usage: " << argv[0] << " -i <imgfile>" << std::endl;
        return 1;
    }
    std::string option = argv[1];
    std::string input = argv[2];
    if (option != "-i") {
        std::cout << "invalid option error. please use " << "<-i> to input" << std::endl;
        return 1;
    }
    /*-----------*/
    

    /*----------imagefile loading------------*/
    std::string imgfile = input;
	TIFF *tif = TIFFOpen(imgfile.c_str(), "r");      
	if (tif == nullptr)
	{
		std::cout << "unable to read image";
        TIFFClose(tif);
        return 1;
	}
    
    //raw image info
    std::cout<< "input imgfile is " << imgfile << std::endl;
    int nchannel;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nchannel);
    int nTotalFrame = TIFFNumberOfDirectories(tif);//TIFFNumberOfDirectories（）means slice number of image
    int nheight, nwidth;
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &nheight);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &nwidth);
    uint16_t ndatatype;
    TIFFGetField(tif,TIFFTAG_SAMPLEFORMAT,&ndatatype);
    size_t nstripSize = TIFFStripSize(tif);//TIFFStripSize（）means voxel number of every slice
    std::cout << "raw tif channel/frame/height/width/datatype/stripsize: " << nchannel << "/" << nTotalFrame << "/";
    std::cout << nheight << "/" << nwidth << "/" << ndatatype << "/" << nstripSize << std::endl;
    
    //create cv mat3d
    int size[3] = {nTotalFrame,nheight,nwidth};
    cv::Mat image(3,size,CV_16UC1,cv::Scalar(0));
    uint16_t maxValue = std::numeric_limits<uint16_t>::lowest();
    uint16_t minValue = std::numeric_limits<uint16_t>::max();
    cv::Point3i maxLoc(0,0,0);    
    cv::Point3i minLoc(0,0,0);    

    //assign
    for (int z = 0; z < nTotalFrame; z++) {
        //request memory for single slice with uint16 with size of TIFFStripSize
        uint16_t* nsliceBuffer = new uint16_t[nstripSize];
        uint16_t nvalue;
        
        //read by row and write from first byte of current slice
        for (int y = 0; y < nheight; y++)
        {
            TIFFReadScanline(tif, (&(nsliceBuffer)[0] + y * int(nwidth)), y);
        }
    
        //assign and find max loc with slice buffer
        int ss = 0;
        for (int y = 0; y < nheight; y++)
        {
            uint16_t* frame_ptr = image.ptr<uint16_t>(z,y);
            for (int x = 0; x < nwidth; x++)
            {
                nvalue = nsliceBuffer[ss];                          
                frame_ptr[x] = nvalue;
                if (maxValue<nvalue) 
                    {
                        maxValue = nvalue;
                        maxLoc.z = z;
                        maxLoc.y = y;
                        maxLoc.x = x;
                    }
                if (minValue>nvalue) 
                    {
                        minValue = nvalue;
                        minLoc.z = z;
                        minLoc.y = y;
                        minLoc.x = x;
                    }
                ss += 1;
            }
        }
        TIFFReadDirectory(tif);
    }
    TIFFClose(tif);

    std::cout << "raw tif max voxel value: " << maxValue;
    std::cout << " and location: " << maxLoc.z << " " << maxLoc.y << " " << maxLoc.x << std::endl;
    std::cout << "raw tif min voxel value: " << minValue;
    std::cout << " and location: " << minLoc.z << " " << minLoc.y << " " << minLoc.x << std::endl;
    std::cout << std::endl;
    
    //info print
    bool print = true;
    if (print){
        std::cout << "new cv" << std::endl;
        infoPrint(image);
        std::cout << std::endl;
    }
    /*------------------------*/



    /*----------equlization------*/
    double p = 80;
    bool equal;
    cv::Mat equalin(3,image.size,CV_16UC1,cv::Scalar(0));
    copyMat16(equalin,image);
    std::cout<< "make a copy of image and do equalization" << std::endl;;
    std::cout<< "equalization parameter is " << p << std::endl;
    findMinMaxLocValueIterative(image,minLoc,maxLoc,minValue,maxValue);
    std::cout << std::endl;
    
    cv::Mat equalout(3,image.size,CV_8UC1,cv::Scalar(0));
    equalizationImageToU8(equalin, equalout, p, equal); 
    if (equal==false){
        std::cout<< "equalization output is " << std::boolalpha << equal << std::endl;
        return 1;
    }
    std::cout << std::endl;
    
    std::cout << "new equalized cv" << std::endl;
    infoPrint(equalout);     
    uint8_t *minValues = new uint8_t[equalout.size[0]];
    uint8_t *maxValues = new uint8_t[equalout.size[0]];
    
    std::string equalfile = imgfile + "_eq" + std::to_string(static_cast<int>(p)) + ".tif";
    bool save;
    saveU8(equalfile, equalout, save);
    std::cout<< "done saving equalization imgfile: " << equalfile << std::endl;  
    std::cout << std::endl;
    /*--------------------*/
    
    
    /*-----------------------gsdt----------------------*/
    

    /*---*/


    /*------------------nms-------------------*/  
    std::cout << "Start processing" << std::endl; 
    std::vector<Marker> markers;
    std::vector<Marker> centers;
    uint32_t cn = 0;
    bool nms = true;
    int iNms = 0;
    uint8_t imax, imin;
    cv::Point3i Lmin, Lmax;
    uint32_t vn = 0;

    //iterative
    std::cout << "Iteratively find centers..." << std::endl;  
    std::cout << std::endl;   
    while (nms){
        iNms += 1;  
        std::cout << "round " << iNms << ": " << std::endl; 
        findMinMaxLocValueIterativeU8(equalout, Lmin, Lmax, imin, imax); 
        cn = countValueNumberU8(equalout, 0);
        std::cout << "image>0 number is " << cn << std::endl;
        nms = (cn >= 60*60*60*7 & imax >= 100); 
        if (!nms){
            std::cout << "Drop this image!" << std::endl;
            std::cout << std::endl;
            break;
            } 
        std::cout << std::endl;

        //reconstruction
        Marker marker;
        marker.x = Lmax.x;
        marker.y = Lmax.y;
        marker.z = Lmax.z;
        std::cout << std::endl;
        markers.push_back(marker);
        flipY(marker, equalout.size[1]);
        std::string markerfile = equalfile + "_i_" + std::to_string(iNms) +
                                "_x_" + std::to_string(static_cast<int>(marker.x)) + 
                                "_y_" +  std::to_string(static_cast<int>(marker.y)) + 
                                "_z_" + std::to_string(static_cast<int>(marker.z)) + "_app2.marker";
        bool writemarker = writeMarker(marker, markerfile);

        std::string swcfile = equalfile + "_i_" + std::to_string(iNms) +
                                "_x_" + std::to_string(static_cast<int>(marker.x)) + 
                                "_y_" +  std::to_string(static_cast<int>(marker.y)) + 
                                "_z_" + std::to_string(static_cast<int>(marker.z)) + "_app2.swc";
        bool recon =  App2OnImagefile(equalfile, swcfile, markerfile);  
        std::cout << std::endl;
        std::cout << "Done reconstruction" << std::endl;
        std::cout << std::endl;
                
        //get center
        bool mask;
        std::string iniTreeFile = equalfile + "_ini.swc";
        std::vector<Node> iniTree;
        bool trueIniTree = checkNodestree(iniTreeFile, iniTree);
        vn += static_cast<uint32_t>(iniTree.size());
        flipY(iniTree,equalout.size[1]);

        if (trueIniTree){
            //parse swc tree
            std::vector<Node> tree;
            bool trueTree = checkNodestree(swcfile, tree);
            flipY(tree,equalout.size[1]);
            
            //check and return furcation
            if (trueTree){
                std::string furfile = swcfile + "_furcation.marker";
                std::vector<Marker> furcations = get_furcation(tree,furfile);
                std::cout << std::endl;
                
                //get center furcation
                if (furcations.size()>0){
                    std::string centerfile = furfile + "_center.marker";
                    Marker center;
                    Marker nearest;
                    bool cen = get_zone_mass_center(equalout, furcations, centerfile, center, nearest);
                    
                    if (cen){
                        centers.push_back(nearest);
                    }
                    std::cout << std::endl;

                }

            } 

        }
        
        //mask
        if (trueIniTree){
            int s = 5;
            std::cout << "start masking ini tree with s=5" << std::endl;
            Mask(equalout, iniTree, s, mask); 
        }else{
            int s = 30;
            std::cout << "start masking marker with s=30" << std::endl;
            Mask(equalout, marker, s, mask); 
        }

        //iterative
        equalfile = equalfile + "_i_" + std::to_string(iNms) + "_mask.tif";
        std::cout << "new image file: " << equalfile << std::endl;
        saveU8(equalfile, equalout, save);
        std::cout<< "done saving masked image file: " << equalfile << std::endl;  
        std::cout << std::endl;
        
    }

    //write all center to files
    if (centers.size()){

        /*----NOTE remember to remove that!----*/
        
        flipY(centers,image.size[1]);
        std::string centerfile = imgfile + "_allcenter.marker";
        bool write = writeMarkerstree(centers, centerfile); 
    }

    delete[] argv;
    return 0;  
}
