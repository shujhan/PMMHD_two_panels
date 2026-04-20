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
#include <set>   
#include <Eigen/Dense>
using namespace Eigen;

#include "initial_distributions.hpp"
#include "Panel.hpp"
#include "FieldStructure.hpp"
#include "Periodizer.hpp"

enum Quadrature {trap, simpsons, last_quad};
enum BoundaryConditions {periodic_bcs, open_bcs};


class AMRStructure {
    std::string sim_dir;
    std::string species_name;
    double x_min, x_max;
    double y_min, y_max;
    // distribution* f0;
    distribution* w0;
    distribution* j0;
    distribution* q0; // q+ or q-

    int initial_height;
    int y_height;
    int max_height;
    double Lx, Ly;
    int npanels_x, npanels_y;
    double initial_dx, initial_dy;


    // panels
    std::vector <Panel> panels;
    std::vector <int> leaf_inds;

    std::vector<double> xs, ys, w0s,j0s,q0s;
    std::vector<double> u1s, u2s, b1s, b2s; // u1: velocity in x; u2: velocity in y. 

    std::vector<double> weights;
    std::vector<double> u_weights, b_weights;

    bool is_initial_mesh_set;
    int minimum_unrefined_index;
    bool need_further_refinement;
    bool do_adaptively_refine; 
    double amr_epsilons;

    // time stepping parameters
    
    // double dt;
    // double t;
    // int n_steps_remesh;
    // int n_steps_diag;
    // int method;


    Quadrature quad;
    // int bcs; // 0 for periodic bc 
    BoundaryConditions bcs;

    Field* calculate_e;
    double greens_epsilon;
    Periodizer* periodizer = nullptr;
   



    std::vector <Panel> old_panels;
    std::vector<double> old_xs, old_ys, old_q0s;

    // interpolation parameters
    double f_beyond_boundary = 0;
    bool use_limiter = false;
    double limit_val = 0.0;
    bool allow_boundary_extrapolation = false;
    bool do_unshear = false;
    bool sqrt_f = false;




// constructor
    public:
    int iter_num = 0;
        AMRStructure(); 
        AMRStructure(std::string sim_dir, std::string species_name, distribution* w0, distribution* j0, distribution* q0,
                    int initial_height, int y_height, int max_height,  double greens_epsilon,
                    double x_min, double x_max, double y_min, double y_max, BoundaryConditions bcs,
                    Quadrature quad, Field* calculate_e, Periodizer* periodizer,
                    bool do_adaptively_refine, double amr_epsilons);
    // destructor
        ~AMRStructure();


    // upadte w and j using new q 
        void recover_wj_from_q(const AMRStructure* q_minus_struct);
        void copy_wj_from(const AMRStructure* other);
        void update_panel_positions_and_q(const std::vector<double>& vx,
                                  const std::vector<double>& vy,
                                  const std::vector<double>& source,
                                  double dt,
                                  int sign);




        
    // generate uniform panel mesh 
        void generate_mesh(std::function<double (double,double)> w0, std::function<double (double,double)> j0,
                            std::function<double (double,double)> q0,
                    bool do_adaptive_refine, bool is_initial_step);

        int create_prerefined_mesh(); 

        void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine);

        void set_leaves_weights();
        void recursively_set_leaves_weights(int panel_ind);

        // TODO 
        void test_panel(int panel_ind, bool verbose);



    // field evaluation 
        // u1 and u2, use u_weights
        int evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // b1 and b2, use b_weights 
        int evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        int compute_rhs_state(std::vector<double>& xs_in, std::vector<double>& ys_in, std::vector<double>& u1s_in, std::vector<double>& u2s_in,
                            std::vector<double>& b1s_in, std::vector<double>& b2s_in, double t_in,
                            std::vector<double>& u1s_grad_x, std::vector<double>& u1s_grad_y, std::vector<double>& u2s_grad_x, std::vector<double>& u2s_grad_y,
                            std::vector<double>& b1s_grad_x, std::vector<double>& b1s_grad_y, std::vector<double>& b2s_grad_x,std::vector<double>& b2s_grad_y,
                            std::vector<double>& source_term
                        );


    // remesh
    void remesh();


    // interpolation
        void interpolate_to_initial_xys(std::vector<double>& fs, std::vector<double>& xs, 
                                                std::vector<double>& ys, int nx, int ny);
        void shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& ys);
        int find_leaf_containing_xy_recursively(double &x, double &y, bool& beyond_boundary, int panel_ind);
        int find_leaf_containing_point_from_neighbor(double& tx, double& ty, bool& beyond_boundary, int leaf_ind, std::set<int>& history);
        void interpolate_from_panel_to_points(std::vector<double>& values_q0s, std::vector<double>& xs, std::vector<double>& ys,
                                                std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val);






        // io
        friend std::ostream& operator<<(std::ostream& os, const AMRStructure& amr);
        void print_amr();
        int write_to_file();
        int write_to_file(bool pre_remesh);

        int write_particles_to_file();
        int write_panels_to_file();
        int write_particles_to_file(bool pre_remesh);
        int write_panels_to_file(bool pre_remesh);


        std::vector<double>& get_xs() { return xs; }
        std::vector<double>& get_ys() { return ys; }
        std::vector<double>& get_q0s() { return q0s; }
        std::vector<double>& get_u_weights() { return u_weights; }
        std::vector<double>& get_b_weights() { return b_weights; }
        std::vector<double>& get_u1s() { return u1s; }
        std::vector<double>& get_u2s() { return u2s; }
        std::vector<double>& get_b1s() { return b1s; }
        std::vector<double>& get_b2s() { return b2s; }
        
        size_t get_xs_size() { return xs.size(); }
        void set_u1s(const std::vector<double>& vals);
        void set_u2s(const std::vector<double>& vals);
        void set_b1s(const std::vector<double>& vals);
        void set_b2s(const std::vector<double>& vals);

        void set_xs(const std::vector<double>& vals);
        void set_ys(const std::vector<double>& vals);
        void set_q0s(const std::vector<double>& vals);

        void set_periodizer(Periodizer* p) { periodizer = p; }

};

#endif
