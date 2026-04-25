#include "AMRSimulation.hpp"
#include <fstream>
#include <iomanip>

MHDDiagnostics AMRSimulation::compute_diagnostics() {
    MHDDiagnostics d{};
    d.iter = iter_num;
    d.t    = t;

    // All quantities are already cached and consistent at this point:
    //   - constructor populated u,B at t=0
    //   - step() refreshes u,B at post-remesh positions before calling this
    std::vector<double>& W   = general_list[0]->get_weights();
    std::vector<double>& u1  = general_list[0]->get_u1s();
    std::vector<double>& u2  = general_list[0]->get_u2s();
    std::vector<double>& B1  = general_list[0]->get_b1s();
    std::vector<double>& B2  = general_list[0]->get_b2s();

    const size_t N = W.size();
    double E_kin = 0.0, E_mag = 0.0, H_C = 0.0;

    #pragma omp parallel for reduction(+:E_kin, E_mag, H_C)
    for (size_t i = 0; i < N; ++i) {
        const double wi = W[i];
        E_kin += 0.5 * wi * (u1[i]*u1[i] + u2[i]*u2[i]);
        E_mag += 0.5 * wi * (B1[i]*B1[i] + B2[i]*B2[i]);
        H_C   +=       wi * (u1[i]*B1[i] + u2[i]*B2[i]);
    }

    d.E_kin = E_kin;
    d.E_mag = E_mag;
    d.E_tot = E_kin + E_mag;
    d.H_C   = H_C;
    return d;
}

int AMRSimulation::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
        f << "iter,t,E_kin,E_mag,E_tot,H_C\n";
        header_written = true;
    } else {
        f.open(path, std::ios::out | std::ios::app);
    }

    if (!f.is_open()) {
        std::cerr << "[diag] failed to open " << path << std::endl;
        return 1;
    }

    f << std::setprecision(16);
    f << d.iter << "," << d.t << ","
      << d.E_kin << "," << d.E_mag << "," << d.E_tot << ","
      << d.H_C << "\n";
    return 0;
}