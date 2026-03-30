#include "AMRSimulation.hpp"

AMRSimulation::AMRSimulation() {}

AMRSimulation::AMRSimulation(std::string sim_dir, std::string deck_address)
{
    this->sim_dir = sim_dir;
    this->deck_address = deck_address;

    // create ptree
    pt::ptree deck;
    load_deck(deck_address, deck);
    cout << "deck address: " << deck_address << endl;

    get_box_t_params(deck);

    // create e solver
    calculate_e = make_field_return_ptr(deck);
    // make_external_field(deck);

    // load vorticity and current density
    try {
        pt::ptree &initial_list_deck = deck.get_child("initial_list");
        auto it = initial_list_deck.begin();

        // get vorticity
        pt::ptree &vorticity = it->second;
        distribution* w0 = make_f0_return_ptr(vorticity);
        ic_list.push_back(w0);

        // get current density 
        ++it; 
        pt::ptree &current_density = it->second;
        distribution* j0 = make_f0_return_ptr(current_density);
        ic_list.push_back(j0);

        distribution* q_plus  = new combined_distribution(w0, j0, +1);
        distribution* q_minus = new combined_distribution(w0, j0, -1);

        ic_list.push_back(q_plus);
        ic_list.push_back(q_minus);

        // each AMR structure has both w0 and j0, q+ is for w+j, q- is for w-j
        general_list.push_back(make_species_return_ptr(vorticity, w0, j0, q_plus));
        general_list.push_back(make_species_return_ptr(current_density, w0, j0, q_minus));

        size_t xs_n = general_list[0]->get_xs_size();

        xs.assign(xs_n,0);
        ys.assign(xs_n,0);
        u1s.assign(xs_n, 0.0);
        u2s.assign(xs_n, 0.0);
        b1s.assign(xs_n, 0.0);
        b2s.assign(xs_n, 0.0);

        u1s_grad_x.assign(xs_n, 0.0);
        u1s_grad_y.assign(xs_n, 0.0);
        u2s_grad_x.assign(xs_n, 0.0);
        u2s_grad_y.assign(xs_n, 0.0);
        b1s_grad_x.assign(xs_n, 0.0);
        b1s_grad_y.assign(xs_n, 0.0);
        b2s_grad_x.assign(xs_n, 0.0);
        b2s_grad_y.assign(xs_n, 0.0);
    } catch(std::exception& e) {
        cout << "Invalid deck format.  No vorticity or current density" << endl;
        return;
    }

    // N_sp = species_list.size();

    iter_num = 0;
    t = 0;

    // print AMR description
    print_sim_setup();
    //write to file
    write_to_file();
}


//destructor
AMRSimulation::~AMRSimulation() {
    for (size_t ii = 0; ii < ic_list.size(); ++ii) {
        delete ic_list[ii];
    }
    for (size_t ii = 0; ii < general_list.size(); ++ii) {
        delete general_list[ii];
    }
    delete calculate_e;
}


