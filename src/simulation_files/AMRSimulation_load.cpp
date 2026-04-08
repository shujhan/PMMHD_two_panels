#include "AMRSimulation.hpp"

int AMRSimulation::load_deck(std::string &deck_address, pt::ptree &deck) {

    try {
        pt::read_json(deck_address, deck);
    } catch(std::exception& err) {
        cout << "unable to open input deck" << endl;
        return 1;
    }
    return 0;
}

int AMRSimulation::get_box_t_params(pt::ptree &deck) {

    std::string project_name = deck.get<std::string>("sim_name", "no_name_found");
    x_min = deck.get<double>("xmin", 0.0), x_max = deck.get<double>("xmax", 0.0);
    y_min = deck.get<double>("ymin", 0.0), y_max = deck.get<double>("ymax",0.0);
    int bcs_int = deck.get<int>("bcs",0);
    
    if (bcs_int < 0 || bcs_int >= 2) {
        cout << "Invalid boundary condition provided in input deck. Valid BCs are: " << endl;
        cout << "0: periodic in x and y, 1: periodic in x, open in y" << endl;
        return 1;
    }
    bcs = static_cast<BoundaryConditions> (bcs_int);
    
    Lx = x_max - x_min;

    int which_quad = deck.get<int>("quadrature",0);
    quad = static_cast<Quadrature>(which_quad);
    method = deck.get<int>("method",0);

    num_steps = deck.get<int>("num_steps", 1);
    n_steps_remesh = deck.get<int>("remesh_period", 1); 
    n_steps_diag = deck.get<int>("diag_period", 1); 
    dt = deck.get<double> ("dt", 0.5); 

    return 0;
}

Field* AMRSimulation::make_field_return_ptr(pt::ptree &deck) {
    
    Field* calculate_field;

    double greens_epsilon = deck.get<double>("greens_epsilon",0.1);
    int use_treecode = deck.get<int>("use_treecode", 0);  // 0 no treecode; 1 with treecode
    // double beta = deck.get<double>("beta", -1.0);
    double mac = deck.get<double>("mac", -1.0); 
    int degree = deck.get<int>("degree", -1); 
    int max_source = deck.get<int>("max_source", 3000);
    int max_target = deck.get<int>("max_target", 3000); 
    
    if (use_treecode > 0) {
        if (bcs!=periodic_bcs) { // open_bcs in y 
            calculate_field = new U_Treecode(Lx, greens_epsilon, mac, degree, max_source, max_target);
            cout << "using treecode: periodic in x and open in y" << endl;
        }
        else { // periodic in y, suing free space kernel 
            KernelMode m = periodic_y;
            calculate_field = new U_Treecode(Lx, greens_epsilon, mac, degree, max_source, max_target);
            calculate_field->set_mode(m);
            cout << "using treecode: periodic in x and y" << endl;
        }
    }
    else {
        if (bcs!=periodic_bcs) { // open_bcs in y 
            calculate_field = new U_DirectSum(Lx, greens_epsilon);
            cout << "using direct sum: : periodic in x and open in y" << endl;
        }
        else {
            KernelMode m = periodic_y;
            calculate_field = new U_DirectSum(Lx, greens_epsilon);
            calculate_field->set_mode(m);
            cout << "using treecode: periodic in x and y" << endl;
        }
    }
    return calculate_field;
}



distribution* AMRSimulation::make_f0_return_ptr(pt::ptree &species_deck_portion) {
    pt::ptree deck = species_deck_portion;

    double kx = 2.0 * M_PI / Lx * deck.get<double>("normalized_wavenumber",1.0);
    double amp = deck.get<double>("amp", 0.0);
    int ics_type = deck.get<int>("ics_type", 1);
    distribution* f0;
    switch (ics_type)
    {
        case 1: 
            f0 = new w0_uniform();
            break;
        case 2: 
            f0 = new j0_uniform();
            break;

        // for current sheet 
        case 3: 
            f0 = new w0_zero();
            break;
        case 4: 
            f0 = new j0_current_sheet(kx,amp);
            break;

        // for alfven
        case 5: 
            f0 = new w0_alfven(kx, amp);
            break;
        case 6: 
            f0 = new j0_alfven(kx,amp);
            break;

        default:
            cout << "Using default initial conditions, all 1s" << endl;
            f0 = new w0_uniform();
            break;
    }
    return f0;
}




AMRStructure* AMRSimulation::make_species_return_ptr(pt::ptree &species_deck_portion, distribution* w0, distribution* j0, distribution* q0) {

    pt::ptree deck = species_deck_portion;
    // get parameters
    std::string sp_name = deck.get_child("name").get_value<std::string>();

    double y_min = deck.get<double>("ymin", -1.0);
    double y_max = deck.get<double>("ymax",1.0);

    int initial_height = deck.get<int>("initial_height",6);
    int y_height = deck.get<int>("y_height",0);
    int max_height = deck.get<int>("max_height", initial_height);
    bool do_adaptively_refine = deck.get<bool> ("do_adaptively_refine", false);
    double amr_epsilons = deck.get<double>("amr_epsilons",0.1);
    

    AMRStructure *species = new AMRStructure{sim_dir, sp_name,
                w0, j0, q0, initial_height, y_height,max_height,
                x_min, x_max, y_min, y_max, bcs, quad, this->calculate_e,
                do_adaptively_refine, amr_epsilons};
    return species;
}
