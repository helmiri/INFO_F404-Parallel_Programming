#include "point.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <stdlib.h>

void generate_colors(int num_points, int *colors, int range,
                     std::mt19937 &rng) {
    std::uniform_int_distribution<int> cols(0, range);
    for (size_t i = 0; i < num_points; i++) {
        colors[i] = cols(rng);
    }
}

void write_points(const char *path, int num_points, float *points, int *colors,
                  int dimensions) {
    std::ofstream file;
    file.open(path);
    for (size_t j = 0; j < num_points; j++) {
        if (colors != nullptr) {
            file << colors[j] << " ";
        }
        for (size_t i = j * dimensions; i < (j + 1) * dimensions; i++) {
            file << points[i] << " ";
        }
        file << "\n";
    }
    file.close();
}

void generate_points(int num_points, float *points, int dimensions,
                     std::mt19937 &rng) {
    std::uniform_real_distribution<float> coords(-100, 100);

    const float round = std::pow(10.0, 3);
    for (size_t i = 0; i < num_points * dimensions; i++) {
        points[i] = std::ceil(coords(rng) * round) / round;
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr
            << "Usage: \n - Environment: \n\t./generator [color_range] "
               "[Num_points] [dimensionality] [output_file]\n - Challenge "
               "file:\n\t./generator[Num_points] [dimensionality] [output_file]"
            << std::endl;
        return EXIT_FAILURE;
    }

    std::random_device rd;
    std::mt19937 rng(rd());

    bool is_ch = (argc == 4); // Is it a challenge set?
    int num_points, dimensions, range;
    const char *output;

    range = atoi(argv[1 - is_ch]);
    num_points = atoi(argv[2 - is_ch]);
    dimensions = atoi(argv[3 - is_ch]);
    output = argv[4 - is_ch];

    int colors[num_points];
    float points[num_points * dimensions];

    generate_points(num_points, points, dimensions, rng);
    if (!is_ch) {
        generate_colors(num_points, colors, range, rng);
        write_points(output, num_points, points, colors, dimensions);
    } else {
        write_points(output, num_points, points, nullptr, dimensions);
    }

    return EXIT_SUCCESS;
}