#include "knn.hpp"
#include <chrono>

int main(int argc, char **argv) {
    int rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (rank == 0) {
        root(world_size, argv[1], argv[2], atoi(argv[3]), argv[4]);
    } else {
        process();
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

int load_file(char *path, std::vector<float> &array) {
    std::ifstream env_file(path);
    int dimensions;
    std::string line;
    while (std::getline(env_file, line)) {
        dimensions = split_line(line, array);
    }
    env_file.close();
    return dimensions;
}

int split_line(std::string line, std::vector<float> &coordinates) {
    std::string delimiter = " ";
    int dimensions = 0;
    char *token = strtok((char *)line.c_str(), " ");
    while (token != NULL) {
        dimensions += 1;
        coordinates.push_back(atof(token));
        token = strtok(NULL, " ");
    }
    return dimensions;
}

int receive_points(std::vector<Point> &array) {
    int result, num_elements, dimensions;
    // Receive number of elements to be received from root process
    result = MPI_Recv(&num_elements, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                      MPI_STATUS_IGNORE);
    if (result != MPI_SUCCESS) {
        return result;
    }

    // Receive the size of the points to properly interpret the array of floats
    result = MPI_Bcast(&dimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        return result;
    }

    // Receive environment points
    std::vector<float> elements(num_elements);
    result = MPI_Scatterv(NULL, NULL, NULL, NULL, elements.data(), num_elements,
                          MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Load points into vector for easier manipulation
    for (size_t i = 0; i < num_elements; i += dimensions) {
        std::vector<float> temp = std::vector<float>(
            elements.begin() + i, elements.begin() + i + dimensions);
        array.push_back(*new Point(dimensions - 1, temp));
    }

    return result;
}

std::vector<int> assign_elements(int num_elements, int world_size,
                                 int dimensions) {

    int num_points = num_elements / dimensions;
    int elements_per_process = num_points / world_size;

    // Process i receives counts[i] elements
    std::vector<int> counts(world_size);
    for (size_t i = 0; i < world_size; i++) {
        counts[i] = dimensions * elements_per_process;
    }

    // Evenly distributes the remainder, if any.
    // (Total number of elements - number of elements already assigned) / Size
    // of a point = Number of points remaining to be assigned
    int remainder =
        (num_elements - (dimensions * elements_per_process * world_size)) /
        dimensions;
    for (size_t i = 0; i < world_size && remainder > 0; i++) {
        counts[i] += dimensions; // Add the size of a point to the elements
                                 // received by process i
        remainder -= 1;
    }
    return counts;
}

std::vector<int> assign_displacements(int world_size,
                                      std::vector<int> &counts) {
    std::vector<int> displacements(world_size);
    displacements[0] = 0;
    int total_displacement = 0;
    // Displacement i is the sum of counts 0-i
    for (size_t i = 1; i < world_size; i++) {
        total_displacement += counts[i - 1];
        displacements[i] = total_displacement;
    }
    return displacements;
}

int scatter_points(std::vector<float> &coordinates, int world_size,
                   int dimensions, std::vector<Point> &points,
                   std::vector<int> &counts, std::vector<int> &displacements) {
    counts = assign_elements(coordinates.size(), world_size, dimensions);
    displacements = assign_displacements(world_size, counts);

    // Inform other processes on the size of the buffers
    int error;
    for (size_t i = 1; i < world_size; i++) {
        error = MPI_Send(&counts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        if (error != MPI_SUCCESS) {
            return error;
        }
    }

    int result;
    // Inform other processes on the size of the point
    result = MPI_Bcast(&dimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        return result;
    }
    // Distribute points
    std::vector<float> elements(counts[0]);
    result = MPI_Scatterv(coordinates.data(), counts.data(),
                          displacements.data(), MPI_FLOAT, elements.data(),
                          counts[0], MPI_FLOAT, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        return result;
    }
    // Convert to Point for ease of use
    for (size_t i = 0; i < elements.size(); i += dimensions) {
        std::vector<float> temp = std::vector<float>(
            elements.begin() + i, elements.begin() + i + dimensions);
        points.push_back(*new Point(dimensions - 1, temp));
    }
    return result;
}

std::vector<float> calculate_distances(std::vector<Point> &points,
                                       Point &point) {
    std::vector<float> distances;
    for (size_t i = 0; i < points.size(); i++) {
        distances.push_back(point.get_distance(points[i]));
        distances.push_back(points[i].get_color());
    }
    return distances;
}

int compare_distances(const void *a, const void *b) {
    Distance *first = (Distance *)a;
    Distance *second = (Distance *)b;
    if (first->distance < second->distance) {
        return -1;
    } else if (first->distance > second->distance) {
        return 1;
    }
    return 0;
}

int get_majority_color(int k, std::vector<Distance> &distances) {
    // Vectors such that counters[i] is the number of distinct_colors[i] found
    // in the k first elements of distances
    std::vector<int> distinct_colors;
    std::vector<int> counters;

    bool check;
    int index;
    int color;
    for (size_t i = 0; i < k; i++) {
        color = distances[i].color;
        // Get counters index of color
        for (size_t i = 0; i < distinct_colors.size(); i++) {
            check = distinct_colors[i] == color;
            if (check) {
                index = i;
                break;
            }
        }

        if (!check) {
            // Create it
            counters.push_back(1);
            distinct_colors.push_back(color);
        } else {
            counters[index] += 1;
        }
    }
    // Get index of max
    int current_index = 0;
    int current_max = 0;
    for (size_t i = 0; i < counters.size(); i++) {
        if (counters[i] > current_max) {
            current_max = counters[i];
            current_index = i;
        }
    }
    return distinct_colors[current_index];
}

void root(int world_size, char *env_file, char *challenge_file, int k,
          char *output_file) {
    std::vector<float> array;
    std::vector<Point> points;

    // Read and scatter environment file
    int dimensions = load_file(env_file, array);
    std::vector<int> counts, displacements;

    if (scatter_points(array, world_size, dimensions, points, counts,
                       displacements) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    array.clear();

    // Read challenge file
    dimensions = load_file(challenge_file, array);
    int num_points = array.size() / dimensions;

    // Broadcast stop condition: Number of points in challenge file
    if (MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Calculate counts and displacements of the vector where the results of the
    // computation will be gathered
    std::vector<int> counts_buffer;
    std::vector<int> displs_buffer;
    for (size_t i = 0; i < counts.size(); i++) {
        counts_buffer.push_back((counts[i] / (dimensions + 1)) * 2);
        displs_buffer.push_back((displacements[i] / (dimensions + 1)) * 2);
    }

    std::ofstream output;
    output.open(output_file);

    // Size of the results vector
    int buffer_size = 0;
    for (size_t i = 0; i < counts.size(); i++) {
        buffer_size += (counts[i] / (dimensions + 1)) * 2;
    }

    for (size_t i = 0; i < array.size(); i += dimensions) {
        // For each point in challenge file, broadcast it
        std::vector<float> temp = std::vector<float>(
            array.begin() + i, array.begin() + i + dimensions);
        if (MPI_Bcast(temp.data(), dimensions, MPI_FLOAT, 0, MPI_COMM_WORLD) !=
            MPI_SUCCESS) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Calculate own points
        Point challenger = Point(temp.size(), temp);
        std::vector<float> distances = calculate_distances(points, challenger);
        distances.resize(buffer_size);

        // Gather results from from other processes and add them to own
        if (MPI_Gatherv(MPI_IN_PLACE, distances.size(), MPI_FLOAT,
                        distances.data(), counts_buffer.data(),
                        displs_buffer.data(), MPI_FLOAT, 0,
                        MPI_COMM_WORLD) != MPI_SUCCESS) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Convert received data to Distances to preserve distance-color
        // association during sort
        std::vector<Distance> col_distance;
        for (size_t i = 0; i < distances.size() - 1; i += 2) {
            Distance dist = {distances[i], distances[i + 1]};
            col_distance.push_back(dist);
        }

        std::qsort(col_distance.data(), col_distance.size(), sizeof(Distance),
                   compare_distances);

        write_result(output, get_majority_color(k, col_distance), challenger);
    }
    output.close();
}

void write_result(std::ofstream &output, int color, Point &challenger) {
    output << color << " ";
    std::vector<float> coordinates = challenger.get_coordinates();
    for (size_t i = 0; i < challenger.get_dimensions(); i++) {
        output << coordinates[i] << " ";
    }
    output << "\n";
}

void process() {
    std::vector<Point> points;
    if (receive_points(points) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int num_points;
    if (MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    std::vector<float> coordinates(points[0].get_dimensions());
    for (size_t i = 0; i < num_points; i++) {
        if (MPI_Bcast(coordinates.data(), points[0].get_dimensions(), MPI_FLOAT,
                      0, MPI_COMM_WORLD) != MPI_SUCCESS) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        Point challenger = Point(coordinates.size(), coordinates);
        std::vector<float> distances = calculate_distances(points, challenger);
        if (MPI_Gatherv(distances.data(), distances.size(), MPI_FLOAT, NULL,
                        NULL, NULL, MPI_FLOAT, 0,
                        MPI_COMM_WORLD) != MPI_SUCCESS) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
}