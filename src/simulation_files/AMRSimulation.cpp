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

    if (bcs == periodic_bcs) {
        periodizer = new Periodizer(x_min, x_max, y_min, y_max, greens_epsilon, calculate_e,
                            80,  // M
                            30,  // m
                            1.4);
        periodizer->precompute_Q();
        // inject into the two species
        general_list[0]->set_periodizer(periodizer);
        general_list[1]->set_periodizer(periodizer);
    }

    iter_num = 0;
    t = 0;

    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();
    std::vector<double> u_ws_step0(xs.size(), 0.0);
    std::vector<double> b_ws_step0(xs.size(), 0.0);
    u_ws_step0 = general_list[0]->get_u_weights();
    b_ws_step0 = general_list[0]->get_b_weights();
    u1s.assign(xs.size(), 0.0); u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0); b2s.assign(xs.size(), 0.0);
    general_list[0]->evaluate_u_field(u1s, u2s, xs, ys, u_ws_step0, t);
    general_list[0]->evaluate_b_field(b1s, b2s, xs, ys, b_ws_step0, t);
    //external field for polarized_alfven wave
    // for (size_t i = 0; i < b1s.size(); ++i) {
    //     b1s[i] +=  2.0/sqrt(5);
    //     b2s[i] +=  1.0/sqrt(5);
    // }

    general_list[0]->set_u1s(u1s);
    general_list[0]->set_u2s(u2s);
    general_list[0]->set_b1s(b1s);
    general_list[0]->set_b2s(b2s);
    general_list[1]->set_u1s(u1s);
    general_list[1]->set_u2s(u2s);
    general_list[1]->set_b1s(b1s);
    general_list[1]->set_b2s(b2s);



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
    delete periodizer;
}


