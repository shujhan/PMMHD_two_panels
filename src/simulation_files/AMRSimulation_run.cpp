#include "AMRSimulation.hpp"

int AMRSimulation::run() {
    while (iter_num < num_steps) {
        step();
    }

    
    write_to_file();
    {
        MHDDiagnostics d = compute_diagnostics();
        write_diagnostics(d);
    }

    return 0;
}