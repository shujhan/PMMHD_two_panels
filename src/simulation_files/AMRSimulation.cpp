#include "AMRSimulation.hpp"

AMRSimulation::AMRSimulation() {}

AMRSimulation::AMRSimulation(std::string sim_dir, std::string deck_address) 
    : need_scatter (false)
    // , outFile("debug.txt")
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
        general_list.push_back(make_species_return_ptr(vorticity, w0));

        // get current density 
        ++it; 
        pt::ptree &current_density = it->second;
        distribution* j0 = make_f0_return_ptr(current_density);
        ic_list.push_back(j0);
        general_list.push_back(make_species_return_ptr(current_density, j0));


    } catch(std::exception& e) {
        cout << "Invalid deck format.  No vorticity or current density" << endl;
        return;
    }

    N_sp = species_list.size();

    
    // initialize e
    iter_num = 0;
    t = 0;
    evaluate_u1_field();
    evaluate_u2_field();
    evaluate_b1_field();
    evaluate_b2_field();
    

    // print AMR description
    print_sim_setup();
    //write to file
    write_to_file();
}


//destructor
AMRSimulation::~AMRSimulation() {
    for (int ii = 0; ii < ic_list.size(); ++ii) {
        delete ic_list[ii];
        delete general_list[ii];
    }
    delete calculate_e;
}


