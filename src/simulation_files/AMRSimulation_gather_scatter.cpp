#include "AMRSimulation.hpp"

int AMRSimulation::gather() {
    
    xs = std::vector<double> ();
    ys = std::vector<double> ();
    ws = std::vector<double> ();
    es = std::vector<double> ();
    species_start = std::vector<size_t> ();
    species_end = std::vector<size_t> ();
    int Nx = 0;

    for (auto &species : species_list) {
        species_start.push_back(Nx);
        Nx += species->xs.size();
        species_end.push_back(Nx);
    }
    xs.reserve(Nx);
    ys.reserve(Nx);
    ws.reserve(Nx);
    es.reserve(Nx);
    
    for (auto &species : species_list) {
        xs.insert( xs.end(), species->xs.begin(), species->xs.end());
        ys.insert( ys.end(), species->ys.begin(), species->ys.end());
        ws.insert( ws.end(), species->ws.begin(), species->ws.end());
        es.insert( es.end(), species->es.begin(), species->es.end());
    }

    need_gather = false;
    return 0;
}

int AMRSimulation::scatter(bool send_e) {
    // distribute to each species
    size_t start_ind = 0;
    for (auto &species : species_list) {
        for (int ii = 0; ii < species->xs.size(); ++ii) {
            species->xs[ii] = xs[start_ind + ii];
            species->ys[ii] = ys[start_ind + ii];
            species->ws[ii] = ws[start_ind + ii];
            if (send_e) {
                species->es[ii] = es[start_ind + ii];
            }
        }
        start_ind += species->xs.size();
    }

    need_scatter = false;
    return 0;
}