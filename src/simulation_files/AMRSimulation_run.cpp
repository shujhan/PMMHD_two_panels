#include "AMRSimulation.hpp"

int AMRSimulation::run() {

    while (iter_num < num_steps) {
        step();
    }
    return 0;
}