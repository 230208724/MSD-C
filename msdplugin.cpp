
/*unsettled part

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
