


std::string reconstruction(const std::string& imgfile = "", const std::string& outdir = "", bool rf = false,
                    const std::string& vaa3dpath = "/home/lyx/software/v3d_external/bin/vaa3d",
                    const std::string& vn2path = "/home/lyx/software/v3d_external/bin/plugins/neuron_tracing/Vaa3D_Neuron2/libvn2.so") {
    std::vector<std::vector<std::vector<unsigned char>>> image;

    if (!imgfile.empty()) {
        // Load the image
        image = load_image(imgfile, rf);
    } else {
        // Handle the case when imgfile is not provided
        // Replace this with your desired behavior
        throw std::runtime_error("Image file is not provided.");
    }

    if (image[0].size() == 2) {
        image = { image };
    }

    // Find the coordinates of the maximum value in the image
    int z, y, x;
    unsigned char maxVal = 0;
    for (size_t i = 0; i < image.size(); ++i) {
        for (size_t j = 0; j < image[i].size(); ++j) {
            for (size_t k = 0; k < image[i][j].size(); ++k) {
                if (image[i][j][k] > maxVal) {
                    maxVal = image[i][j][k];
                    z = i;
                    y = j;
                    x = k;
                }
            }
        }
    }

    std::cout << "Coordinates: " << y << ", " << x << ", " << z << std::endl;
    
    std::string markerfile = outdir + "/" + imgfile.substr(imgfile.find_last_of('/') + 1) + "_y" +
                             std::to_string(y) + "_x" + std::to_string(x) + "_z" + std::to_string(z) + ".marker";
    write_marker({{ x, y, z }}, markerfile);

    if (!imgfile.empty()) {
        std::string newfile = outdir + "/" + imgfile.substr(imgfile.find_last_of('/') + 1) + ".swc";
        std::string vn2_str = vaa3dpath + " -x " + vn2path + 
    " -f app2 -i " + imgfile + " -o " + newfile + " -p " + markerfile +
    " 0 AUTO 0 0 0 1 5 0 0 0";
    // Execute the Vaa3D command
    int status = system(vn2_str.c_str());
    if (status != 0) {
    // Handle the case when the command execution fails
    throw std::runtime_error("Failed to execute Vaa3D command.");
    }
    return newfile;
    } else {
    // Handle the case when imgfile is not provided
    // Replace this with your desired behavior
    throw std::runtime_error("Image file is not provided.");
    }
}

int main() {
try {
std::string imgfile = ""; // Provide the image file path here
std::string outdir = ""; // Provide the output directory path here
bool rf = false; // Specify the flip_tif value here
std::string vaa3dpath = "/home/lyx/software/v3d_external/bin/vaa3d";
std::string vn2path = "/home/lyx/software/v3d_external/bin/plugins/neuron_tracing/Vaa3D_Neuron2/libvn2.so";
    std::string swcFile = get_swc(imgfile, outdir, rf, vaa3dpath, vn2path);
    std::cout << "SWC file generated: " << swcFile << std::endl;
} catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
}

return 0;


        std::string imgfile_name = imgfile.substr(imgfile.find_last_of('/') + 1);
        std::string imgfile_extension = imgfile.substr(imgfile.find_last_of('.') + 1);

        newfile = outdir + "/" + imgfile_name + "_gsdt" + std::to_string(u8) + "." + imgfile_extension;

}
