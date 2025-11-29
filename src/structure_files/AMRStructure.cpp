#include "AMRStructure.hpp"


AMRStructure::AMRStructure() {}

AMRStructure::AMRStructure(std::string sim_dir, std::string species_name,
                            distribution* f0,
                            int initial_height, int y_height, int max_height, 
                            double x_min, double x_max, double y_min, double y_max, 
                            Quadrature quad,
                            bool do_adaptively_refine, double amr_epsilons)
                           : f0(f0),
                           initial_height(initial_height), y_height(y_height),
                           max_height(max_height), 
                           x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), 
                           quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine)
{
    this->sim_dir = sim_dir;
    this->species_name = species_name;
    Lx = x_max - x_min;
    Ly = y_max - y_min;
    npanels_x = int(pow(2, initial_height));
    npanels_y = int(pow(2, initial_height + y_height));
    initial_dx = Lx / npanels_x;
    initial_dy = Ly / npanels_y;
    this->amr_epsilons = amr_epsilons;
    bcs = 0;

    bool is_initial_step = true;
    generate_mesh([&](double x, double y) { return (*f0)(x,y); }, do_adaptively_refine, is_initial_step);





}


//destructor
AMRStructure::~AMRStructure() = default;


