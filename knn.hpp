#ifndef KNN_HPP
#define KNN_HPP

#include "point.hpp"
#include <cstdlib>
#include <fstream>
#include <mpi.h>
#include <stdio.h>
#include <string>
#include <vector>

/**
 * @brief Associates a distance with the color of the point to keep track of
 *        which distance belongs to what color after sorting the array of
 * Distances
 *
 */
typedef struct {
    float distance;
    int color;
} Distance;

/**
 * @brief Read the contents of the file specified by path
 *
 * @param path Path to the file to be read
 * @param array Vector where the values will be stored
 * @return int The number of elements on each line in the file (size of
 *             coordinates + color)
 */
int load_file(char *path, std::vector<float> &array);

/**
 * @brief Split a line by whitespace, converst the tokens to float and stores
 *        them in the vector
 *
 * @param line The line to be split
 * @param coordinates Vector where the values will be stored
 * @return int The number of elements on each line in the file (size of
 *             coordinates + color)
 */
int split_line(std::string line, std::vector<float> &coordinates);

/**
 * @brief Evenly distributes the points in the environment between all processes
 *
 * @param coordinates Vector containing the points in the environment in the
 *                    following format: (c_0,
 * p_00,...,p_n0,...,ci,p_0i,...,p_ni) where c_0-c_i are the colors of points
 * p_0-p_i and p_0i-p_ni are the coordinates of dimension n of points p_0-p_i
 * @param world_size Number of processes
 * @param dimensions Size of coordinates + color. Needed to properly segment the
 *                   vector.
 * @param points Vector where the points assigned to root process will be stored
 * @param counts Vector where counts[i] contains the elements sent to process i
 * @param displacements Vector where displacements[i] contains the index k in
 *                      coordinates such that counts[i] elements will be sent to
 * process i starting from coordinates[k]
 * @return int The result of the MPI functions called
 */
int scatter_points(std::vector<float> &coordinates, int world_size,
                   int dimensions, std::vector<Point> &points,
                   std::vector<int> &counts, std::vector<int> &displacements);

/**
 * @brief Receives environment points from root process
 *
 * @param array Vector where the points will be stored
 * @return int The result of the MPI functions called
 */
int receive_points(std::vector<Point> &array);

/**
 * @brief Comparison function to be used with std::qsort to sort the values
 *        returned by the processes
 *
 * @param a Distance
 * @param b Distance
 * @return int -1 if a < b, 1 if a > b, 0 if a == b
 */
int compare_distances(const void *a, const void *b);

/**
 * @brief Get the majority color
 *
 * @param k The number of closest points to be considered
 * @param distances Vector of Distance
 * @return int The majority color
 */
int get_majority_color(int k, std::vector<Distance> &distances);

/**
 * @brief Calculate distances from point to all points in vector
 *
 * @param points Environment points
 * @param point Point to be classified
 * @return std::vector<float> Vector containing the distances and colors of the
 *                            points with format (d_0, c_0,...,d_i, c_i) where
 * d_i and c_i are respectively the distance and color of point i
 */
std::vector<float> calculate_distances(std::vector<Point> &points,
                                       Point &point);

/**
 * @brief Calculates the number of environment points that each process will
 * receive
 *
 * @param num_elements Total number of elements (coordinates + colors)
 * @param world_size Number of processes
 * @param dimensions Size of each point (coordinates + color)
 * @return std::vector<int> Vector where element i is the number of elements
 *                          that process i will receive
 */
std::vector<int> assign_elements(int num_elements, int world_size,
                                 int dimensions);

/**
 * @brief Calculates the displacements from where to count elements
 *
 * @param world_size The number of processes
 * @param counts Vector where element i is the number of elements that process i
 *               will receive
 * @return std::vector<int> Vector where element i is the the index in the
 *                          vector of points from where to count[i] elements
 */
std::vector<int> assign_displacements(int world_size, std::vector<int> &counts);

/**
 * @brief Write the result of 1 computation to a file
 *
 * @param output Destination file
 * @param color Color of the point
 * @param challenger Point to be classified
 */
void write_result(std::ofstream &output, int color, Point &challenger);

/**
 * @brief Root process
 *
 * @param world_size Number of processes
 * @param env_file Path to environment file
 * @param challenge_file Path to challenge file
 * @param k Number of points to be considred for classification
 * @param output_file Path to file where the results will be saved
 */
void root(int world_size, char *env_file, char *challenge_file, int k,
          char *output_file);

/**
 * @brief Non-root process
 *
 */
void process();

#endif