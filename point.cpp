#include "point.hpp"

Point::Point(int dimensions, std::vector<float> &coordinates) {
    // Element at index 0 is the color of the point if dimensions < coordinates.size()
    for (size_t i = coordinates.size() - dimensions; i < coordinates.size();
         i++) {
        this->coordinates.push_back(coordinates[i]);
    }
    this->color = (dimensions == coordinates.size()) ? -1 : coordinates[0];
}

float Point::get_distance(Point &other) {
    float dist = 0;
    std::vector<float> other_coordinates = other.get_coordinates();
    for (size_t i = 0; i < coordinates.size(); i++) {
        dist += std::pow((coordinates[i] - other_coordinates[i]), 2);
    }

    return std::sqrt(dist);
}

std::vector<float> &Point::get_coordinates() { return coordinates; }

int Point::get_color() { return color; }

int Point::get_dimensions() { return coordinates.size(); }
