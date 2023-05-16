#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <array>
#include <fstream>
#include <opencv2/opencv.hpp> // Make sure to install OpenCV and link it correctly

// try-except grammer:
// try {
// } catch (const std::exception& ex) {
//     std::cerr << "Error: " << ex.what() << std::endl;
// }

#include "../config.h"

#include "image_processor/preprocess.h"

#include "image_processor/reconstruction.h"

#include "swc_processor/swc_basic_process.h"

bool MaxOnImage(cv2::Mat& equalout, Marker& marker){
    int val = cv2::function(&&equalout);//find max value on 3d image 'equalout' with cv2
    if (val<50){
        return false;
    }
    else{
        vector<int> xyz cv2::function(&&equalout, &val);////find the first max value point on 3d image 'equalout' with cv2
        marker.x = xyz[0];
        marker.y = xyz[1];
        marker.z = xyz[2];
        return true;
    }
}

bool MaskByMarker(cv2:Mat& image, auto& node, int s){
    //related values set to 0
    image[node.z-s:node.z+s,node.z-s:node.z+s,node.z-s:node.z+s] = 0
    return true;
}

bool Mask(cv2:Mat& image, std::string& inifile, Marker& marker){   
    std::vector<Node> initree;
    bool check = checkNodestree(&inifile, &initree);
    if (check){
        for (auto &node: initree){
            MaskByMarker(&image,&node,int s=1);
        }
    else{
        MaskByMarker(&image,&marker,int s=20);
    }    
    } 
}


int calc_mass(int x, int y, int z, int r, cv::Mat image) {
    int zz = image.size[0];
    int yy = image.size[1];
    int xx = image.size[2];
    cv::Range z_range(std::max(0, z - r), std::min(z + r + 1, zz));
    cv::Range y_range(std::max(0, y - r), std::min(y + r + 1, yy));
    cv::Range x_range(std::max(0, x - r), std::min(x + r + 1, xx));
    cv::Mat sub_image = image(z_range, y_range, x_range);
    int mass = cv::sum(sub_image)[0];
    return mass;
}

double calc_distance(auto& marker1,auto& marker2){
    double dx = (marker1.x - marker2.x)>0 ? : (marker1.x - marker2.x) : (marker2.x - marker1.x);
    double dy = (marker1.y - marker2.y)>0 ? : (marker1.y - marker2.y) : (marker2.y - marker1.y);
    double dz = (marker1.z - marker2.z)>0 ? : (marker1.z - marker2.z) : (marker2.z - marker1.z);
    double ds = dx*dx + dy*dy + dz*dz;
    return ds;
}   

bool find_nearest(vector<Marker>& markers, Marker& marker, Marker& nearest){
    double mds = __INTMAX__;
    for (m : markers){
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

bool get_zone_mass_center(cv2::Mat &equalout, 
                            std::vector<Marker> &furcations,
                            const std::string& centerfile, 
                            Marker& center) {
double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double m_ = 0.0;

    for (const Marker& marker : furcations) {
        double x = marker.x;
        double y = marker.y;
        double z = marker.z;
        double r = marker.r;
        double mass = calc_mass(x, y, z, r, &equalout);
        x_ += x * mass;
        y_ += y * mass;
        z_ += z * mass;
        m_ += mass;
    }

    if (m_==0){
        vector<Marker> markers;
        writeMarkerstree(&markers, &centerfile);
        return false;
        }  
    else{
        Marker center;
        center.x = x_/ m_; 
        center.y = y_ / m_; 
        center.z = z_ / m_; 

        Marker nearest;
        bool near = find_nearest(&furcations, &center, &nearest);

        center.x = nearest.x; 
        center.y = nearest.y; 
        center.z = nearest.z; 
        //std::cout << cx << " " << cy << " " << cz << std::endl;
        writeMarkerstree(&center, &centerfile);
        return true;
    }   
}

bool main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "usage：" << argv[0] << " -i <imgfile>" << std::endl;
        return false;
    }

    std::string option = argv[1];
    std::string input = argv[2];

    if (option != "-i") {
        std::cout << "invalid option error. please use " << "<-i> to input" << std::endl;
        return false;
    }

    //equalization
    std::string imgfile = input;
    float p = 0.8;
    std::string equalfile = imgfile + "_eq" + std::to_string(p) + ".tif";
    cv::Mat equalout;
    bool equal = equalizationImagefileToU8(&input, &p, &equalfile, equalout);

    //gsdt

    //nms  
    std::vector<Marker> markers;
    std::vector<Marker> centers;
    int i = 0;
    while (i<3){
        Marker marker;
        bool nms;
        nms = MaxOnImage(cv2::Mat& equalout, Marker& marker);
        i += 1;   
        if nms:            
            markers.push_back(marker);
            //reconstruction
            std::string swcfile = imgfile + "_i_" + i.to_string() +"_x_" + marker.x + "_y_" +  marker.y + "_z_" +  marker.z + "_app2.swc";
            std::string markerfile = imgfile + "_i_" + i.to_string() + "_x_" + marker.x + "_y_" +  marker.y + "_z_" +  marker.z + "_app2.marker";
            bool recon =  App2OnImagefile(&equalfile, &swcfile, &markerfile);   
            //mask
            std::string inifile = imgfile + '_ini.swc';
            bool mask = Mask(cv2:Mat& equalout, std::string& inifile, Marker& marker);    
            //parse
            std::vector<Node> tree = parseSwcfile(&swcfile)
            //check and return furcation
            std::string furfile = swcfile + "_furcation.marker";
            std::vector<Marker> furcations = get_furcation(&tree,&furfile);
            if (furcations.size>0){
                //get center furcation
                std::string centerfile = furfile + "_center.marker";
                Marker center;
                bool cen = get_zone_mass_center(&equalout, &furcations, &centerfile, &center);
                if (cen){
                    centers.push_back(center);
                }
            }            
    }
    //write all center to files
    std::string centerfile = imgfile + "_allcenter.marker";
    writeMarkerstree(&centers, &centerfile);   
}



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




