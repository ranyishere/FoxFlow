#ifndef DGGML_FUNCTIONS_HPP
#define DGGML_FUNCTIONS_HPP

namespace particle_functions {
    double distance(const std::array<double, 2> &pos1, const  std::array<double, 2> &pos2) {
        double dx = pos1[0] - pos2[0];
        double dy = pos1[1] - pos2[1];
        return std::sqrt(dx * dx + dy * dy);
    }
}

#endif //DGGML_FUNCTIONS_HPP
