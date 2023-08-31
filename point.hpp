#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>
#include <math.h>
#include <vector>

class Point {
  private:
    std::vector<float> coordinates;
    int color;

  public:
    /**
     * @brief Construct a new Point object
     *
     * @param dimensions Size of a point's coordinates
     * @param coordinates The coordinates
     */
    Point(int dimensions, std::vector<float> &coordinates);

    /**
     * @brief Computes the distance between this Point and another
     *
     * @param other Target point
     * @return float The distance
     */
    float get_distance(Point &other);

    int get_color();
    int get_dimensions();
    std::vector<float> &get_coordinates();
};

#endif