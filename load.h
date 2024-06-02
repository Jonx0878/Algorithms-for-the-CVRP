#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <curl/curl.h>

struct InstanceData {
    int dimension = -1;
    int capacity = -1;
    std::vector<int> demand;
    std::vector<std::vector<int>> distances; // Upper triangular matrix without the diagonal
    std::vector<std::vector<int>> routes;
    std::vector<std::pair<double, double>> coords;
    int upper_bound = -1;
};

// Gets the letters of the URL for downloading files
std::string get_letters_before_separator(const std::string& filename) {
    std::string letters;
    for (char ch : filename) {
        if (ch == '_' || ch == '-' || (ch >= '0' && ch <= '9')) {
            break;
        }
        letters.push_back(ch);
    }
    // Catch Exception cases
    if (letters == "Loggi" || letters == "ORTEC") {
        return "D";
    }
    else if (letters == "Antwerp" || letters == "Brussels" || letters == "Flanders" || letters == "Ghent" || letters == "Leuven") {
        return "XXL";
    }

    return letters;
}

// Callback function to write data from the HTTP response to a file
size_t write_data(void* ptr, size_t size, size_t nmemb, FILE* stream) {
    return fwrite(ptr, size, nmemb, stream);
}

// Function to download file from URL and save it locally
bool download_file(const std::string& url, const std::string& filepath) {
    CURL* curl;
    FILE* fp;
    CURLcode res;

    curl = curl_easy_init();
    if (curl) {
        // Open file for writing
        errno_t err = fopen_s(&fp, filepath.c_str(), "wb");
        if (err != 0 || fp == nullptr) {
            std::cerr << "Error opening file for writing." << std::endl;
            return false;
        }

        // Set up cURL options
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);

        // Perform the request
        res = curl_easy_perform(curl);

        // Cleanup
        curl_easy_cleanup(curl);
        fclose(fp);

        if (res != CURLE_OK) {
            remove(filepath.c_str()); // Delete the file if download failed
            std::cerr << "Error downloading file from URL." << std::endl;
            return false;
        }
        return true;
    }

    std::cerr << "Error initializing cURL." << std::endl;
    return false;
}

// Function to download instance file (.vrp) and solution file (.sol)
bool download_instance_and_solution(const std::string& filename) {
    const std::string instance_download_location = "Instance Data/Instances/";
    const std::string solution_download_location = "Instance Data/Solutions/";

    // Construct URLs for downloading instance (.vrp) and solution (.sol) files
    std::string base_url = "http://vrp.galgos.inf.puc-rio.br/media/com_vrp/instances/";
    std::string letters = get_letters_before_separator(filename);
    base_url += letters + "/" + filename;
    std::string instance_url = base_url + ".vrp";
    std::string solution_url = base_url + ".sol";

    // Download instance file
    if (!download_file(instance_url, instance_download_location + filename + ".vrp")) {
        return false;
    }

    // Download solution file
    if (!download_file(solution_url, solution_download_location + filename + ".sol")) {
        // If solution file download fails, delete the downloaded instance file
        remove((instance_download_location + filename + ".vrp").c_str());
        return false;
    }

    return true;
}


// Reads a .vrp file
bool load_instance_data(const std::string& filename, InstanceData& data) {
    std::string instance_path = "Instance Data/Instances/" + filename + ".vrp";
    std::ifstream file(instance_path);
    if (!file.is_open()) {
        std::cout << "Instance file not found locally. Attempting to download...\n";
        if (!download_instance_and_solution(filename)) {
            std::cerr << "Error: Unable to download instance file from URL.\n";
            exit(1);
        }
        file.open(instance_path);
    }

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << '\n';
        exit(1);
    }

    std::string line;
    bool calculate_distances = true;
    while (std::getline(file, line)) {
        if (line.find("DIMENSION") != std::string::npos) {
            size_t pos = line.find(":");
            data.dimension = std::stoi(line.substr(pos + 1));
            data.distances.resize(data.dimension);
            for (int i = 0; i < data.dimension; ++i) {
                data.distances[i].resize(data.dimension);
            }
        }
        else if (line.find("CAPACITY") != std::string::npos) {
            size_t pos = line.find(":");
            data.capacity = std::stoi(line.substr(pos + 1));
        }
        else if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos) {
            if (data.dimension == -1) {
                std::cerr << "Error: DIMENSION should be defined before EDGE_WEIGHT_SECTION.\n";
                exit(1);
            }
            for (int j = 1; j < data.dimension; ++j) {
                for (int i = 0; i < j; ++i) {
                    file >> data.distances[i][j]; // Store distances in upper triangular matrix
                }
            }
            calculate_distances = false;
        }
        else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
            if (data.dimension == -1) {
                std::cerr << "Error: DIMENSION should be defined before NODE_COORD_SECTION.\n";
                exit(1);
            }
            std::vector<std::pair<double, double>> coords(data.dimension);
            for (int i = 0; i < data.dimension; ++i) {
                int node;
                double x, y;
                file >> node >> x >> y;
                coords[node - 1] = { x, y };
            }
            data.coords = coords;
            if (calculate_distances) {
                for (int i = 0; i < data.dimension - 1; ++i) {
                    for (int j = i + 1; j < data.dimension; ++j) {
                        double dx = coords[i].first - coords[j].first;
                        double dy = coords[i].second - coords[j].second;
                        int dist = int(std::round(std::sqrt(dx * dx + dy * dy)));
                        data.distances[i][j] = dist;
                    }
                }
            }
        }
        else if (line.find("DEMAND_SECTION") != std::string::npos) {
            if (data.dimension == -1) {
                std::cerr << "Error: DIMENSION should be defined before DEMAND_SECTION.\n";
                exit(1);
            }
            data.demand.resize(data.dimension);
            for (int i = 0; i < data.dimension; ++i) {
                std::getline(file, line);
                int node, demand;
                sscanf_s(line.c_str(), "%d %d", &node, &demand);
                data.demand[node - 1] = demand;
            }
            break; // No need to continue parsing
        }
    }

    file.close();

    if (data.dimension == -1) {
        std::cout << "Something went wrong when downloading the files :(\n";
        std::remove(instance_path.c_str());
        std::remove(("Instance Data/Solutions/" + filename + ".sol").c_str());
        return false;
    }
    return true;
}

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

// Loads solution data
void load_solution_data(const std::string& filename, InstanceData& vrp_data) {
    std::string solution_path = "Instance Data/Solutions/" + filename + ".sol";

    std::ifstream file(solution_path);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open solution file " << filename << std::endl;
        return; // Return early if file opening fails
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("Route #") != std::string::npos) {
            std::vector<int> route;
            std::stringstream ss(line.substr(line.find(":") + 2));
            int node;
            while (ss >> node) {
                route.push_back(node);
            }
            vrp_data.routes.push_back(route);
        }
        else if (line.find("Cost") != std::string::npos) {
            vrp_data.upper_bound = std::stoi(line.substr(line.find(" ") + 1));
        }
    }

    file.close();
}


// Loads instance data
InstanceData load_instance(const std::string& filename) {
    InstanceData data;
    std::cout << filename << "\n";

    auto start = std::chrono::high_resolution_clock::now(); // Start measuring loading time
    bool inst_loaded = load_instance_data(filename, data);
    if (inst_loaded) {
        load_solution_data(filename, data);
    }

    auto stop = std::chrono::high_resolution_clock::now(); // Stop measuring loading time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // Calculate loading time in milliseconds
    std::cout << "Loading time: " << duration.count() << " milliseconds\n";
    return data;
}
