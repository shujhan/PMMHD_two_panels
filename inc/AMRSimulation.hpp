#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <string> 
#include <vector>               // std::vector
#include <fstream>
#include <iostream> 

#include "/sw/pkgs/arc/stacks/gcc/10.3.0/boost/1.78.0/include/boost/property_tree/ptree.hpp"        //property_tree::ptree, property_tree::read_json
#include "/sw/pkgs/arc/stacks/gcc/10.3.0/boost/1.78.0/include/boost/property_tree/json_parser.hpp"
namespace pt = boost::property_tree;
using namespace std;

#include "FieldStructure.hpp"
#include "AMRStructure.hpp"
#include "Periodizer.hpp"



struct MHDDiagnostics {
    int    iter;
    double t;
    double E_kin;
    double E_mag;
    double E_tot;
    double H_C;
};

class AMRSimulation {
    ofstream outFile;
    string sim_dir;
    string deck_address;

    // for simulation box
    double x_min, x_max, y_min, y_max;
    double Lx, Ly;

    // load 
    vector<distribution*> ic_list;
    vector<AMRStructure*> general_list;
    // int N_sp;
    // std::vector<size_t> species_start, species_end;

    std::vector<double> xs, ys;
    std::vector<double> u1s, u2s, b1s, b2s;
    std::vector<double> u1s_grad_x, u1s_grad_y;
    std::vector<double> u2s_grad_x, u2s_grad_y;
    std::vector<double> b1s_grad_x, b1s_grad_y;
    std::vector<double> b2s_grad_x, b2s_grad_y;

    // // field parameters
    Quadrature quad;
    Field* calculate_e;
    BoundaryConditions bcs; // 0 periodic in y, 1 open in y 
    Periodizer* periodizer = nullptr;
    double greens_epsilon;

    // External uniform field for polarized-Alfvén test
    // Set to 0 for tests without external field
    double B0x = 0.0;
    double B0y = 0.0;

    int method; // 0: rk4, otherwise euler

    // pusher parameters
    int num_steps;
    int n_steps_remesh;
    int n_steps_diag;
    double dt;

    int iter_num;
    double t;







    public:
        // constructor
        AMRSimulation();
        AMRSimulation(std::string sim_dir, std::string deck_address);
        // destructor
        ~AMRSimulation();

        // functions for building panels for each species
        int load_deck(std::string &deck_address, pt::ptree &deck);
        int get_box_t_params(pt::ptree &deck);
        Field* make_field_return_ptr(pt::ptree &deck);
        // void make_external_field(pt::ptree &deck);
        distribution* make_f0_return_ptr(pt::ptree &species_deck_portion);
        AMRStructure* make_species_return_ptr(pt::ptree &species_deck_portion, distribution* w0, distribution* j0, distribution* q0);


        // functions for pushing particles
        // int gather();
        // int scatter(bool send_e);
        int run();
        int step();
        int euler();
        int rk4();


        void compute_stage_rhs(double t_stage,
                       std::vector<double>& stage_xs,    
                       std::vector<double>& stage_ys,   
                       std::vector<double>& stage_u1s,
                       std::vector<double>& stage_u2s,
                       std::vector<double>& stage_b1s,
                       std::vector<double>& stage_b2s,
                       std::vector<double>& stage_source);


        MHDDiagnostics compute_diagnostics();   
        int write_diagnostics(const MHDDiagnostics& d); 


        // function for remeshing 


        
        

        int write_to_file();
        int write_to_file(bool pre_remesh);
        void print_sim_setup();








};














#endif 
