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
    int N_sp;
    std::vector<size_t> species_start, species_end;

    std::vector<double> xs, ys, ws, es, rho_ws;
    // std::vector<double> species_qs, species_ms, species_qms;

    // field parameters
    Quadrature quad;
    ElectricField* calculate_e;


    // pusher parameters
    int relativistic; // 0 non-relativistic; 1 relativistic
    int num_steps;
    int n_steps_remesh;
    int n_steps_diag;
    double dt;
    bool need_scatter;
    bool need_gather;

    int bcs; // 0 for periodic 
    int iter_num, t;



    public:
        // constructor
        AMRSimulation();
        AMRSimulation(std::string sim_dir, std::string deck_address);
        // destructor
        ~AMRSimulation();

        // functions for building panels for each species
        int load_deck(std::string &deck_address, pt::ptree &deck);
        int get_box_t_params(pt::ptree &deck);
        ElectricField* make_field_return_ptr(pt::ptree &deck);
        // void make_external_field(pt::ptree &deck);
        distribution* make_f0_return_ptr(pt::ptree &species_deck_portion);
        AMRStructure* make_species_return_ptr(pt::ptree &species_deck_portion, distribution* f0);




        // functions for field evaluation
        int evaluate_field_uniform_grid(double t);
        // int evaluate_field(std::vector<double>& es_local, std::vector<double>& xs_local, std::vector<double>& q_ws_local, double t);





        // functions for pushing particles
        int gather();
        int scatter(bool send_e);
        int run();
        int step();
        int euler();
        int rk4(bool get_4th_e);


        // function for remeshing 


        
        

        int write_to_file();
        int write_to_file(bool pre_remesh);
        void print_sim_setup();








};














#endif 
