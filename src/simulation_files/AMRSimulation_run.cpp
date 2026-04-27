#include "AMRSimulation.hpp"

int AMRSimulation::run() {
    // ---- Conservation diagnostics (every step) ----
    {
        MHDDiagnostics d = compute_diagnostics();
        write_diagnostics(d);
        // std::cout << "[diag] iter=" << d.iter
        //           << " t="    << d.t
        //           << " E="    << d.E_tot
        //           << " H_C="  << d.H_C << std::endl;
    }

    while (iter_num < num_steps) {
        step();
    }
    return 0;
}