#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <array>
#include <fstream>
//#include <filesystem>
#include "/usr/include/opencv2/opencv.hpp" // Make sure to install OpenCV and link it correctly

// try-except grammer:
// try {
// } catch (const std::exception& ex) {
//     std::cerr << "Error: " << ex.what() << std::endl;
// }

#include "config.h"

#include "image_processor/image_process.h"

#include "image_processor/reconstruction.h"

#include "swc_processor/swc_basic_process.h"

// bool MaxOnImage(cv::Mat& image, Marker& marker){
//     double minValue, maxValue;
//     cv::Point minLoc, maxLoc;
//     cv::minMaxLoc(image, &minValue, &maxValue, &minLoc, &maxLoc);//find max value on 3d image 'equalout' with cv2
//     if (maxValue<50){
//         return false;
//     }
//     else{
//         marker.x = maxLoc[0];
//         marker.y = maxLoc[1];
//         marker.z = maxLoc[2];
//         return true;
//     }
// }
bool MaxOnImage(cv::Mat& image, Marker& marker){
    float maxVal = std::numeric_limits<float>::min();
    cv::Point3i maxLoc;

    for (int z = 0; z < image.size[0]; ++z) {
        for (int y = 0; y < image.size[1]; ++y) {
            for (int x = 0; x < image.size[2]; ++x) {
                float value = image.at<float>(z, y, x);
                if (value > maxVal) {
                    maxVal = value;
                    maxLoc = cv::Point3i(x, y, z);
                }
            }
        }
    }

    if (maxVal < 50){
        return false;
    }
    else{
        marker.x = maxLoc.x;
        marker.y = maxLoc.y;
        marker.z = maxLoc.z;
        return true;
    }
}

template<typename NodeType>
bool MaskByMarker(cv::Mat& image, NodeType& node, int s) {
    // Related values set to 0
    int x = static_cast<int>(std::ceil(node.x));
    int y = static_cast<int>(std::ceil(node.y));
    int z = static_cast<int>(std::ceil(node.z));
    int x_start = x - s;
    int y_start = y - s;
    int z_start = z - s;
    int x_end = x + s;
    int y_end = y + s;
    int z_end = z + s;
    setSubImage(image, x_start, y_start, z_start, x_end, y_end, z_end);
    return true;
}

bool Mask(cv::Mat& image, std::string& inifile, Marker& marker){   
    std::vector<Node> initree;
    bool check = checkNodestree(inifile, initree);
    if (check){
        int s = 1;
        for (auto node: initree){
            bool mask = MaskByMarker(image,node,s);
        }
    }
    else{
        int s = 20;
        MaskByMarker(image,marker,s);
    }    
} 

int calc_mass(int x, int y, int z, int r, const cv::Mat& image) {
    //note this shape is z/y/x
    int zz = image.size[0];
    int yy = image.size[1];
    int xx = image.size[2];
    
    int x_start = std::max(0, x - r);
    int y_start = std::max(0, y - r);
    int z_start = std::max(0, z - r);
    int x_end = std::min(x + r + 1, xx);
    int y_end = std::min(y + r + 1, yy);
    int z_end = std::min(z + r + 1, zz);
    
    cv::Mat subimage = extractSubImage(image, x_start, y_start, z_start, x_end, y_end, z_end);
    
    int mass = cv::sum(subimage)[0];
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
                            const std::string& centerfile, 
                            Marker& center) {
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double m_ = 0.0;

    for (const Marker marker : furcations) {
        double x = marker.x;
        double y = marker.y;
        double z = marker.z;
        double r = marker.r;
        double mass = calc_mass(x, y, z, r, image);
        x_ += x * mass;
        y_ += y * mass;
        z_ += z * mass;
        m_ += mass;
    }

    if (m_==0){
        std::vector<Marker> markers;
        writeMarkerstree(markers, centerfile);
        return false;
        }  
    else{
        Marker center;
        center.x = x_/ m_; 
        center.y = y_ / m_; 
        center.z = z_ / m_; 

        Marker nearest;
        bool near = find_nearest(furcations, center, nearest);

        center.x = nearest.x; 
        center.y = nearest.y; 
        center.z = nearest.z; 
        //std::cout << cx << " " << cy << " " << cz << std::endl;
        writeMarkerstree(center, centerfile);
        return true;
    }   
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "usage：" << argv[0] << " -i <imgfile>" << std::endl;
        return 1;
    }

    std::string option = argv[1];
    std::string input = argv[2];

    if (option != "-i") {
        std::cout << "invalid option error. please use " << "<-i> to input" << std::endl;
        return 1;
    }

    //equalization
    std::string imgfile = input;
    float p = 80;
    std::string equalfile = imgfile + "_eq" + std::to_string(static_cast<int>(p)) + ".tif";
    cv::Mat equalout;
    bool equal = equalizationImagefileToU8(input, p, equalfile, equalout);

    //gsdt

    //nms  
    std::vector<Marker> markers;
    std::vector<Marker> centers;
    int i = 0;
    while (i<3){
        Marker marker;
        bool nms;
        nms = MaxOnImage(equalout, marker);
        i += 1;   
        if (nms){            
            markers.push_back(marker);
            //reconstruction
            std::string swcfile = imgfile + "_i_" + std::to_string(i) +"_x_" + std::to_string(marker.x) + "_y_" +  std::to_string(marker.y) + "_z_" + std::to_string(marker.z) + "_app2.swc";
            std::string markerfile = imgfile + "_i_" + std::to_string(i) + "_x_" + std::to_string(marker.x) + "_y_" + std::to_string(marker.y) + "_z_" + std::to_string(marker.z) + "_app2.marker";
            bool recon =  App2OnImagefile(equalfile, swcfile, markerfile);   
            //mask
            std::string inifile = imgfile + "_ini.swc";
            bool mask = Mask(equalout, inifile, marker);    
            //parse
            std::vector<Node> tree = parseSwcfile(swcfile);
            //check and return furcation
            std::string furfile = swcfile + "_furcation.marker";
            std::vector<Marker> furcations = get_furcation(tree,furfile);
            if (furcations.size()>0){
                //get center furcation
                std::string centerfile = furfile + "_center.marker";
                Marker center;
                bool cen = get_zone_mass_center(equalout, furcations, centerfile, center);
                if (cen){
                    centers.push_back(center);
                }
            }
        }            
    }
    //write all center to files
    std::string centerfile = imgfile + "_allcenter.marker";
    bool write = writeMarkerstree(centers, centerfile); 

    return 0;  
}



/* unsettled parts

// Function to calculate contract values
vector<double> calc_contract_values(const string& imgfile, const string& swcfile) {
    ImageType image;
    SWCType tree;

    if (!imgfile.empty()) {
        // Assuming you have a function to load the image, replace the following line with the appropriate code
        image = loadImage(imgfile, false);
    }

    if (!swcfile.empty()) {
        // Assuming you have a function to parse SWC files, replace the following line with the appropriate code
        tree = parseSWC(swcfile);
    }

    vector<vector<int>> array(tree.begin(), tree.end());
    vector<vector<int>> xyz;
    for (const auto& row : array) {
        xyz.push_back({ row[2], row[3], row[4] });
    }

    // Remove rows with negative values
    xyz.erase(remove_if(xyz.begin(), xyz.end(), [](const vector<int>& row) {
        return any_of(row.begin(), row.end(), [](int val) {
            return val < 0;
        });
    }), xyz.end());

    // Remove rows with values greater than the image shape
    const vector<int> shape = { image.zsize(), image.xsize(), image.ysize() };
    xyz.erase(remove_if(xyz.begin(), xyz.end(), [&shape](const vector<int>& row) {
        return any_of(row.begin(), row.end(), [&shape](int val, int i) {
            return val > shape[i];
        });
    }), xyz.end());

    // Get pixel values at the coordinates
    vector<double> values;
    for (const auto& coord : xyz) {
        values.push_back(image(coord[0], coord[1], coord[2]));
    }

    return values;
}

// Function to get image statistics
map<string, double> get_statistics(const string& imgfile) {
    ImageType image;

    if (!imgfile.empty()) {
        // Assuming you have a function to load the image, replace the following line with the appropriate code
        image = loadImage(imgfile, false);
    }

    map<string, double> stat_dict;

    // Compute aggregate statistics
    vector<string> agglist = { "max", "min", "mean", "std", "median" };
    for (const auto& agg : agglist) {
        if (agg == "median") {
            double median = 0.0;
            vector<double> sorted_values(image.begin(), image.end());
            sort(sorted_values.begin(), sorted_values.end());
            if (sorted_values.size() % 2 == 0) {
                median = (sorted_values[sorted_values.size() / 2 - 1] + sorted_values[sorted_values.size() / 2]) / 2.0;
            } else {
                median = sorted_values[sorted_values.size() / 2];
            }
            stat_dict[agg] = median;
        } else {
            double value = 0.0;
            if (agg == "max") {
                value = *max_element(image.begin(), image.end());
            } else if (agg == "min") {
                value = *min_element(image.begin(), image.end());
            } else if (agg == "mean") {
                value = accumulate(image.begin(), image.end(), 0.0) / image.size();
            } else if (agg == "std") {
                double mean = accumulate(image.begin(), image.end(), 0.0) / image.size();
                double sum = 0.0;
                for (const auto& val : image) {
                sum += pow(val - mean, 2);
                }
                value = sqrt(sum / image.size());
            }
            stat_dict[agg] = value;
            // Compute percentiles
            vector<string> percentiles = { "99", "95", "90", "80", "70", "60", "50", "20", "10", "5", "1" };
            for (const auto& percentile : percentiles) {
                double value = 0.0;
                if (!image.empty()) {
                    int rank = ceil((stod(percentile) / 100.0) * (image.size() - 1));
                    vector<double> sorted_values(image.begin(), image.end());
                    nth_element(sorted_values.begin(), sorted_values.begin() + rank, sorted_values.end());
                    value = sorted_values[rank];
                }
                stat_dict[percentile] = value;
            }
            return stat_dict;
        }
    }
}

*/

