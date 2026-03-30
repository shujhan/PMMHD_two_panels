#include "AMRSimulation.hpp"

int AMRSimulation::evaluate_u1_field(std::vector<double>& es_local, std::vector<double>& xs_local, std::vector<double>& q_ws_local, double t) {

    std::vector<double> xtemp_cpy(xs_local), xtemp_cpy2(xs_local);
    (*calculate_e)(es_local.data(), xtemp_cpy.data(), es_local.size(),
                    xtemp_cpy2.data(), q_ws_local.data(), xtemp_cpy.size());
    if (use_external_field) {
        (*calculate_e_external)(es_local.data(), xs_local.data(), es.size(), t);
    }
    return 0;
}


