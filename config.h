// config.h

#ifndef CONFIG_H
#define CONFIG_H

// 定义全局变量
#include <string>
std::string __PROJECTPATH__ = "/media/lyx/LaCie/MSD-C";
std::string __V3DSO__ = "/home/lyx/software/v3d_external/bin/vaa3d";
std::string __APP2SO__ = "/home/lyx/software/v3d_external/bin/plugins/neuron_tracing/Vaa3D_Neuron2/libvn2.so";
std::string __GSDTSO__ = "/home/lyx/software/v3d_external/bin/plugins/image_filters/Grayscale_Image_Distance_Transform/libgsdt.so";
double __INTMAX__ = 2147483647;

struct Node {
    int idx;
    int type;
    double x;
    double y;
    double z;
    double r;
    int parent;
};

struct Marker {
    double x;
    double y;
    double z;
    double r = 0.0;
    int shape = 0;
    std::string name = '';
    std::string comment = '';
    int red = 255;
    int green = 0;
    int blue = 0;
};



//在其他源文件中，包含 config.h 头文件后，可以直接使用全局变量

#endif
