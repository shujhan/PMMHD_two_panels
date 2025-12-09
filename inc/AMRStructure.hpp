#ifndef AMRSTRUCTURE_HPP
#define AMRSTRUCTURE_HPP

#include <string> 
#include <vector>
#include <math.h>
#include <functional> 
#include <fstream>
#include <iostream>             // std::cout, std::endl
#include <assert.h> 
#include <algorithm>

#include "initial_distributions.hpp"
#include "Panel.hpp"

enum Quadrature {trap, simpsons, last_quad};


class AMRStructure {
    std::string sim_dir;
    std::string species_name;
    double x_min, x_max;
    double y_min, y_max;
    distribution* f0;

    int initial_height;
    int y_height;
    int max_height;
    double Lx, Ly;
    int npanels_x, npanels_y;
    double initial_dx, initial_dy;


    // panels
    std::vector <Panel> panels;
    std::vector <int> leaf_inds;

    std::vector<double> xs, ys, fs, ws;
    std::vector<double> u1s, u2s; // u1: velocity in x; u2: velocity in y. 


    bool is_initial_mesh_set;
    int minimum_unrefined_index;
    bool need_further_refinement;
    bool do_adaptively_refine; 
    double amr_epsilons;

    Quadrature quad;
    int bcs; // 0 for periodic bc 



    // interpolation parameters
    double f_beyond_boundary = 0;





// constructor
    public:
        AMRStructure(); 
        AMRStructure(std::string sim_dir, std::string species_name, distribution* f0, 
                    int initial_height, int y_height, int max_height,
                    double x_min, double x_max, double y_min, double y_max,
                    Quadrature quad,
                    bool do_adaptively_refine, double amr_epsilons);
    // destructor
        ~AMRStructure();



    // generate uniform panel mesh 
        void generate_mesh(std::function<double (double,double)> f,
                    bool do_adaptive_refine, bool is_initial_step);

        int create_prerefined_mesh(); 

        void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine);

        void set_leaves_weights();
        void recursively_set_leaves_weights(int panel_ind);

        // TODO 
        void test_panel(int panel_ind, bool verbose);




    // interpolation
        void AMRStructure::interpolate_to_initial_xys(std::vector<double>& fs, std::vector<double>& xs, 
                                                std::vector<double>& ys, int nx, int ny);








        int write_particles_to_file(int iter_num);
        int write_panels_to_file(int iter_num);
        int write_particles_to_file(bool pre_remesh, int iter_num);
        int write_panels_to_file(bool pre_remesh, int iter_num);

        // io
        friend std::ostream& operator<<(std::ostream& os, const AMRStructure& amr);
        void print_amr();
        int write_to_file(int iter_num);
        int write_to_file(bool pre_remesh, int iter_num);
        void print_panel_points();





};

#endif
