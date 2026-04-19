#include "AMRSimulation.hpp"

// #define DEBUG

int AMRSimulation::step() {
    iter_num += 1;
    general_list[0]->iter_num = iter_num;
    general_list[1]->iter_num = iter_num;
    std::cout << "step " << iter_num << std::endl;
    // if (need_gather) {
    //     // gather from species (need at first and after every remesh)
    //     // gather vorticity and current density xs and ys
    //     gather();
    // }

    if (method > 0) {
        euler();
    }
    else {
        rk4();
    }
    
    t += dt;

    // must remesh to calculate w0 and j0 based on q+ and q-
    general_list[0]->remesh();
    general_list[1]->remesh();
    // if (iter_num % n_steps_remesh == 0) {
    //     remesh();
    // }

    // after remesh, calculate w0 and j0, also reset leave weights
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();

    // Then in future for amr, do amr here for new w0 and j0. 



    // if dump : write to file
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }

    return 0;
}


int AMRSimulation::euler() {
    cout << "enter euler" << endl;
    // get information from first panel 
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    if(iter_num > 1) {
        xs = general_list[0]->get_xs();
        ys = general_list[0]->get_ys();
        std::vector<double> u_ws(xs.size(), 0.0);
        std::vector<double> b_ws(xs.size(), 0.0);
        u_ws = general_list[0]->get_u_weights();
        b_ws = general_list[0]->get_b_weights();
        general_list[0]->evaluate_u_field(u1s, u2s, xs, ys, u_ws, t);
        general_list[0]->evaluate_b_field(b1s, b2s, xs, ys, b_ws, t);

        //external field for alfven wave
        // for (size_t i = 0; i < b1s.size(); ++i) {
        //     b1s[i] +=  1.0;
        // }

        //external field for polarized_alfven wave
        // for (size_t i = 0; i < b1s.size(); ++i) {
        //     b1s[i] +=  2/sqrt(5);
        //     b2s[i] +=  1/sqrt(5);
        // }
        
        general_list[0]->set_u1s(u1s);
        general_list[0]->set_u2s(u2s);
        general_list[0]->set_b1s(b1s);
        general_list[0]->set_b2s(b2s);

        general_list[1]->set_u1s(u1s);
        general_list[1]->set_u2s(u2s);
        general_list[1]->set_b1s(b1s);
        general_list[1]->set_b2s(b2s);
    }
    else{
        xs = general_list[0]->get_xs();
        ys = general_list[0]->get_ys();
        u1s = general_list[0]->get_u1s();
        u2s = general_list[0]->get_u2s();
        b1s = general_list[0]->get_b1s();
        b2s = general_list[0]->get_b2s();
    }


    // source terms using finte difference 
    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    std::vector<double> source_term(xs.size(), 0.0);
    general_list[0]->compute_rhs_state(xs, ys, u1s, u2s, b1s, b2s, t, 
                u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
                b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, source_term);

    // build advection velocities
    std::vector<double> vx_plus(xs.size()), vy_plus(xs.size());
    std::vector<double> vx_minus(xs.size()), vy_minus(xs.size());

    for (size_t ii = 0; ii < xs.size(); ++ii) {
        vx_plus[ii]  = u1s[ii] - b1s[ii];
        vy_plus[ii]  = u2s[ii] - b2s[ii];
        vx_minus[ii] = u1s[ii] + b1s[ii];
        vy_minus[ii] = u2s[ii] + b2s[ii];
    }
    // q+ panel
    general_list[0]->update_panel_positions_and_q(vx_plus, vy_plus, source_term, dt, +1);
    // q- panel
    general_list[1]->update_panel_positions_and_q(vx_minus, vy_minus, source_term, dt, -1);

    return 0;
}


int AMRSimulation::rk4() {
    cout << "enter rk4" << endl;

    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();

    // copy xs,ysq0s from two panels so that even we change panel, we will push based on their original values
    std::vector<double> original_q_plus_xs(xs.size()), original_q_plus_ys(xs.size()), original_q_plus_q0s;
    std::vector<double> original_q_minus_xs(xs.size()), original_q_minus_ys(xs.size()), original_q_minus_q0s;
    original_q_plus_xs = general_list[0]->get_xs();
    original_q_plus_ys = general_list[0]->get_ys();
    original_q_plus_q0s = general_list[0]->get_q0s();

    original_q_minus_xs = general_list[1]->get_xs();
    original_q_minus_ys = general_list[1]->get_ys();
    original_q_minus_q0s = general_list[1]->get_q0s();




    const int xs_size = (int)xs.size();

    // Stage storage: k1, k2, k3, k4 for x, y, w and j
    std::vector<double> k1_u1s(xs_size, 0.0), k1_u2s(xs_size, 0.0), k1_b1s(xs_size, 0.0), k1_b2s(xs_size, 0.0), k1_source_term(xs.size(), 0.0);
    std::vector<double> k2_u1s(xs_size, 0.0), k2_u2s(xs_size, 0.0), k2_b1s(xs_size, 0.0), k2_b2s(xs_size, 0.0), k2_source_term(xs.size(), 0.0);
    std::vector<double> k3_u1s(xs_size, 0.0), k3_u2s(xs_size, 0.0), k3_b1s(xs_size, 0.0), k3_b2s(xs_size, 0.0), k3_source_term(xs.size(), 0.0);
    std::vector<double> k4_u1s(xs_size, 0.0), k4_u2s(xs_size, 0.0), k4_b1s(xs_size, 0.0), k4_b2s(xs_size, 0.0), k4_source_term(xs.size(), 0.0);

    // Temporary positions
    std::vector<double> k2_xs_plus(xs_size), k2_ys_plus(xs_size), k2_q_plus(xs_size), k2_xs_minus(xs_size), k2_ys_minus(xs_size), k2_q_minus(xs_size);
    std::vector<double> k3_xs_plus(xs_size), k3_ys_plus(xs_size), k3_q_plus(xs_size), k3_xs_minus(xs_size), k3_ys_minus(xs_size), k3_q_minus(xs_size);
    std::vector<double> k4_xs_plus(xs_size), k4_ys_plus(xs_size), k4_q_plus(xs_size), k4_xs_minus(xs_size), k4_ys_minus(xs_size), k4_q_minus(xs_size);


    std::vector<double> k1_u_ws(xs.size(), 0.0), k2_u_ws(xs.size(), 0.0), k3_u_ws(xs.size(), 0.0), k4_u_ws(xs.size(), 0.0);
    std::vector<double> k1_b_ws(xs.size(), 0.0), k2_b_ws(xs.size(), 0.0), k3_b_ws(xs.size(), 0.0), k4_b_ws(xs.size(), 0.0);

    if(iter_num > 1) {
        k1_u_ws = general_list[0]->get_u_weights();
        k1_b_ws = general_list[0]->get_b_weights();
        general_list[0]->evaluate_u_field(k1_u1s, k1_u2s, xs, ys, k1_u_ws, t);
        general_list[0]->evaluate_b_field(k1_b1s, k1_b2s, xs, ys, k1_b_ws, t);

        //external field for alfven wave
        // for (size_t i = 0; i < k1_b1s.size(); ++i) {
        //     k1_b1s[i] +=  1.0;
        // }

        //external field for polarized_alfven wave
        // for (size_t i = 0; i < k1_b1s.size(); ++i) {
        //     k1_b1s[i] +=  2/sqrt(5);
        //     k1_b1s[i] +=  1/sqrt(5);
        // }


        general_list[0]->set_u1s(k1_u1s);
        general_list[0]->set_u2s(k1_u2s);
        general_list[0]->set_b1s(k1_b1s);
        general_list[0]->set_b2s(k1_b2s);

        general_list[1]->set_u1s(k1_u1s);
        general_list[1]->set_u2s(k1_u2s);
        general_list[1]->set_b1s(k1_b1s);
        general_list[1]->set_b2s(k1_b2s);
    }
    else{
        k1_u1s = general_list[0]->get_u1s();
        k1_u2s = general_list[0]->get_u2s();
        k1_b1s = general_list[0]->get_b1s();
        k1_b2s = general_list[0]->get_b2s();
    }

    // ------------------------------------------------------------
    // Stage 1 : k1 = (k1_u1s, k1_u2s, k1_b1s, k1_b2s, k1_source_term)
    // ------------------------------------------------------------

    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    // source term to update q0
    general_list[0]->compute_rhs_state(xs, ys, k1_u1s, k1_u2s, k1_b1s, k1_b2s, t, 
                u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
                b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k1_source_term);


    // build advection velocities
    std::vector<double> k1_vx_plus(xs.size()), k1_vy_plus(xs.size());
    std::vector<double> k1_vx_minus(xs.size()), k1_vy_minus(xs.size());

    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k1_vx_plus[ii]  = k1_u1s[ii] - k1_b1s[ii];
        k1_vy_plus[ii]  = k1_u2s[ii] - k1_b2s[ii];
        k1_vx_minus[ii] = k1_u1s[ii] + k1_b1s[ii];
        k1_vy_minus[ii] = k1_u2s[ii] + k1_b2s[ii];
    }

    // q+ panel 
    // double k1_dt = 0.5 * dt;
    // general_list[0]->update_panel_positions_and_q(k1_vx_plus, k1_vy_plus, k1_source_term, k1_dt, +1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k2_xs_plus[ii] = original_q_plus_xs[ii] + 0.5 * dt * k1_vx_plus[ii]; 
        k2_ys_plus[ii] = original_q_plus_ys[ii] + 0.5 * dt * k1_vy_plus[ii]; 
        k2_q_plus[ii] = original_q_plus_q0s[ii] + 0.5 * dt * k1_source_term[ii]; 
    }

    // q- panel
    // general_list[1]->update_panel_positions_and_q(k1_vx_minus, k1_vy_minus, k1_source_term, k1_dt, -1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k2_xs_minus[ii] = original_q_minus_xs[ii] + 0.5 * dt * k1_vx_minus[ii]; 
        k2_ys_minus[ii] = original_q_minus_ys[ii] + 0.5 * dt * k1_vy_minus[ii]; 
        k2_q_minus[ii] = original_q_minus_q0s[ii] - 0.5 * dt * k1_source_term[ii]; 
    }

    general_list[0]->set_xs(k2_xs_plus);
    general_list[0]->set_ys(k2_ys_plus);
    general_list[0]->set_q0s(k2_q_plus);

    general_list[1]->set_xs(k2_xs_minus);
    general_list[1]->set_ys(k2_ys_minus);
    general_list[1]->set_q0s(k2_q_minus);



    // remesh first to calculate w0 and j0, then reset leave weights for future u and b evaluations
    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();


    // ------------------------------------------------------------
    // Stage 2 : k2 = (k2_u1s, k2_u2s, k2_b1s, k2_b2s, k2_source_term)
    // ------------------------------------------------------------
    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();
    k2_u_ws = general_list[0]->get_u_weights();
    k2_b_ws = general_list[0]->get_b_weights();
    general_list[0]->evaluate_u_field(k2_u1s, k2_u2s, xs, ys, k2_u_ws, t + 0.5 * dt);
    general_list[0]->evaluate_b_field(k2_b1s, k2_b2s, xs, ys, k2_b_ws, t + 0.5 * dt);

    //external field for alfven wave
    // for (size_t i = 0; i < k2_b1s.size(); ++i) {
    //     k2_b1s[i] +=  1.0;
    // }

    //external field for polarized_alfven wave
    // for (size_t i = 0; i < k2_b1s.size(); ++i) {
    //     k2_b1s[i] +=  2/sqrt(5);
    //     k2_b1s[i] +=  1/sqrt(5);
    // }

    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    // source term to update q0
    general_list[0]->compute_rhs_state(xs, ys, k2_u1s, k2_u2s, k2_b1s, k2_b2s, t + 0.5 * dt, 
                u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
                b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k2_source_term);

    // build advection velocities
    std::vector<double> k2_vx_plus(xs.size()), k2_vy_plus(xs.size());
    std::vector<double> k2_vx_minus(xs.size()), k2_vy_minus(xs.size());

    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k2_vx_plus[ii]  = k2_u1s[ii] - k2_b1s[ii];
        k2_vy_plus[ii]  = k2_u2s[ii] - k2_b2s[ii];
        k2_vx_minus[ii] = k2_u1s[ii] + k2_b1s[ii];
        k2_vy_minus[ii] = k2_u2s[ii] + k2_b2s[ii];
    }
    // double k2_dt = 0.5 * dt;
    // q+ panel 
    // general_list[0]->update_panel_positions_and_q(k2_vx_plus, k2_vy_plus, k2_source_term, k2_dt, +1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k3_xs_plus[ii] = original_q_plus_xs[ii] + 0.5 * dt * k2_vx_plus[ii]; 
        k3_ys_plus[ii] = original_q_plus_ys[ii] + 0.5 * dt * k2_vy_plus[ii]; 
        k3_q_plus[ii] = original_q_plus_q0s[ii] + 0.5 * dt * k2_source_term[ii]; 
    }
    // q- panel
    // general_list[1]->update_panel_positions_and_q(k2_vx_minus, k2_vy_minus, k2_source_term, k2_dt, -1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k3_xs_minus[ii] = original_q_minus_xs[ii] + 0.5 * dt * k2_vx_minus[ii]; 
        k3_ys_minus[ii] = original_q_minus_ys[ii] + 0.5 * dt * k2_vy_minus[ii]; 
        k3_q_minus[ii] = original_q_minus_q0s[ii] - 0.5 * dt * k2_source_term[ii]; 
    }

    general_list[0]->set_xs(k3_xs_plus);
    general_list[0]->set_ys(k3_ys_plus);
    general_list[0]->set_q0s(k3_q_plus);

    general_list[1]->set_xs(k3_xs_minus);
    general_list[1]->set_ys(k3_ys_minus);
    general_list[1]->set_q0s(k3_q_minus);

    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();



    // ------------------------------------------------------------
    // Stage 3 : k3 = (k3_u1s, k3_u2s, k3_b1s, k3_b2s, k3_source_term)
    // ------------------------------------------------------------
    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();
    k3_u_ws = general_list[0]->get_u_weights();
    k3_b_ws = general_list[0]->get_b_weights();
    general_list[0]->evaluate_u_field(k3_u1s, k3_u2s, xs, ys, k3_u_ws, t + 0.5 * dt);
    general_list[0]->evaluate_b_field(k3_b1s, k3_b2s, xs, ys, k3_b_ws, t + 0.5 * dt);

    //external field for alfven wave
    // for (size_t i = 0; i < k3_b1s.size(); ++i) {
    //     k3_b1s[i] +=  1.0;
    // }

    //external field for polarized_alfven wave
    // for (size_t i = 0; i < k3_b1s.size(); ++i) {
    //     k3_b1s[i] +=  2/sqrt(5);
    //     k3_b1s[i] +=  1/sqrt(5);
    // }

    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    // source term to update q0
    general_list[0]->compute_rhs_state(xs, ys, k3_u1s, k3_u2s, k3_b1s, k3_b2s, t + 0.5 * dt, 
                u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
                b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k3_source_term);

    // build advection velocities
    std::vector<double> k3_vx_plus(xs.size()), k3_vy_plus(xs.size());
    std::vector<double> k3_vx_minus(xs.size()), k3_vy_minus(xs.size());

    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k3_vx_plus[ii]  = k3_u1s[ii] - k3_b1s[ii];
        k3_vy_plus[ii]  = k3_u2s[ii] - k3_b2s[ii];
        k3_vx_minus[ii] = k3_u1s[ii] + k3_b1s[ii];
        k3_vy_minus[ii] = k3_u2s[ii] + k3_b2s[ii];
    }
    // double k3_dt = dt;
    // q+ panel 
    // general_list[0]->update_panel_positions_and_q(k3_vx_plus, k3_vy_plus, k3_source_term, k3_dt, +1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k4_xs_plus[ii] = original_q_plus_xs[ii] + dt * k3_vx_plus[ii]; 
        k4_ys_plus[ii] = original_q_plus_ys[ii] + dt * k3_vy_plus[ii]; 
        k4_q_plus[ii] = original_q_plus_q0s[ii] + dt * k3_source_term[ii]; 
    }
    // q- panel
    // general_list[1]->update_panel_positions_and_q(k3_vx_minus, k3_vy_minus, k3_source_term, k3_dt, -1);
    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k4_xs_minus[ii] = original_q_minus_xs[ii] + dt * k3_vx_minus[ii]; 
        k4_ys_minus[ii] = original_q_minus_ys[ii] + dt * k3_vy_minus[ii]; 
        k4_q_minus[ii] = original_q_minus_q0s[ii] - dt * k3_source_term[ii]; 
    }

    general_list[0]->set_xs(k4_xs_plus);
    general_list[0]->set_ys(k4_ys_plus);
    general_list[0]->set_q0s(k4_q_plus);

    general_list[1]->set_xs(k4_xs_minus);
    general_list[1]->set_ys(k4_ys_minus);
    general_list[1]->set_q0s(k4_q_minus);

    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();


    // ------------------------------------------------------------
    // Stage 4 : k4 = (k4_u1s, k4_u2s, k4_b1s, k4_b2s, k4_source_term)
    // ------------------------------------------------------------

    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();
    k4_u_ws = general_list[0]->get_u_weights();
    k4_b_ws = general_list[0]->get_b_weights();
    general_list[0]->evaluate_u_field(k4_u1s, k4_u2s, xs, ys, k4_u_ws, t + dt);
    general_list[0]->evaluate_b_field(k4_b1s, k4_b2s, xs, ys, k4_b_ws, t + dt);

    //external field for alfven wave
    // for (size_t i = 0; i < k4_b1s.size(); ++i) {
    //     k4_b1s[i] +=  1.0;
    // }

    //external field for polarized_alfven wave
    // for (size_t i = 0; i < k4_b1s.size(); ++i) {
    //     k4_b1s[i] +=  2/sqrt(5);
    //     k4_b1s[i] +=  1/sqrt(5);
    // }

    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    // source term to update q0
    general_list[0]->compute_rhs_state(xs, ys, k4_u1s, k4_u2s, k4_b1s, k4_b2s, t + dt, 
                u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
                b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k4_source_term);

    // build advection velocities
    std::vector<double> k4_vx_plus(xs.size()), k4_vy_plus(xs.size());
    std::vector<double> k4_vx_minus(xs.size()), k4_vy_minus(xs.size());

    for (size_t ii = 0; ii < xs.size(); ++ii) {
        k4_vx_plus[ii]  = k4_u1s[ii] - k4_b1s[ii];
        k4_vy_plus[ii]  = k4_u2s[ii] - k4_b2s[ii];
        k4_vx_minus[ii] = k4_u1s[ii] + k4_b1s[ii];
        k4_vy_minus[ii] = k4_u2s[ii] + k4_b2s[ii];
    }




    for (int i = 0; i < xs_size; ++i) {
        original_q_plus_xs[i] += (dt / 6.0) * (k1_vx_plus[i] + 2.0 * k2_vx_plus[i] + 2.0 * k3_vx_plus[i] + k4_vx_plus[i]);
        original_q_plus_ys[i] += (dt / 6.0) * (k1_vy_plus[i] + 2.0 * k2_vy_plus[i] + 2.0 * k3_vy_plus[i] + k4_vy_plus[i]);
        original_q_plus_q0s[i] += (dt / 6.0) * (k1_source_term[i] + 2.0 * k2_source_term[i] + 2.0 * k3_source_term[i] + k4_source_term[i]);
        
        original_q_minus_xs[i] += (dt / 6.0) * (k1_vx_minus[i] + 2.0 * k2_vx_minus[i] + 2.0 * k3_vx_minus[i] + k4_vx_minus[i]);
        original_q_minus_ys[i] += (dt / 6.0) * (k1_vy_minus[i] + 2.0 * k2_vy_minus[i] + 2.0 * k3_vy_minus[i] + k4_vy_minus[i]);
        original_q_minus_q0s[i] -= (dt / 6.0) * (k1_source_term[i] + 2.0 * k2_source_term[i] + 2.0 * k3_source_term[i] + k4_source_term[i]);
    }

    // copy back to 2 panels and then reconstruct based on this updation
    general_list[0]->set_xs(original_q_plus_xs);
    general_list[0]->set_ys(original_q_plus_ys);
    general_list[0]->set_q0s(original_q_plus_q0s);

    general_list[1]->set_xs(original_q_minus_xs);
    general_list[1]->set_ys(original_q_minus_ys);
    general_list[1]->set_q0s(original_q_minus_q0s);

    // out of rk4, in step:
    // general_list[0]->remesh();
    // general_list[1]->remesh();
    // general_list[0]->recover_wj_from_q(general_list[1]);
    // general_list[1]->copy_wj_from(general_list[0]);
    // general_list[0]->set_leaves_weights();
    // general_list[1]->set_leaves_weights();
    return 0;
}
















// int AMRSimulation::rk4_step(bool get_4th_e) {
// // initialize rk4 vectors
//     std::vector<double> xk = xs;
//     std::vector<double> vs(xs.size()), v1(xs.size()), v2(xs.size()), v3(xs.size()), v4(xs.size() );
//     std::vector<double> p1 = ps, p2(xs.size()), p3(xs.size()), p4(xs.size());;
//     std::vector<double> ef1 = es, ef2(xs.size()), ef3(xs.size()), ef4(xs.size());
//     int N = xs.size();


// //   math_vector p1, p2, p3, p4;
// //   math_vector f1(total_num_points), f2(total_num_points), f3(total_num_points), f4(total_num_points);
// //   math_vector tempx;
  
//     // k1 = (v1,f1) = F(un) = F(xn,pn) = (vn, q/m E(xn) )
//   // v1 = vs = vn
//     // calculate_E_mq(a1, xs, xs, q_ws, L, epsilon);
//     // v1 = vs;

//     // get vs and v1 from p1
//     for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
//         for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
//             vs[xi] = ps[xi]/ species_ms[sp_i];
//             v1[xi] = vs[xi];
//         }
//     }

//     for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
// #if TESTFLAG
//     outFile << "In RK4, the " << sp_i << "th species_qms: " << species_qms[sp_i] << std::endl;
// #endif
//         for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {

//     // k2 = (v2,a2) = F(un + h/2 k1) = F(xn + delt/2 v1, vn + delt/2 a1)
//     //              = ( vn + delt a1 / 2, q/m E(xn + delt v1 /2) )
//             ef1[xi] *= species_qms[sp_i];

//             // k->2
//             // p2.push_back(ps[xi] + 0.5 * dt * ef1[xi]);
//             p2[xi] = ps[xi] + 0.5 * dt * ef1[xi];
//             v2[xi] = p2[xi]/ species_ms[sp_i];
//             xk[xi] = xs[xi] + 0.5 * dt * p1[xi];
//         }
//     }
    
//     // get x2 = xn + delt v1 / 2
    

//     // auto start = high_resolution_clock::now();
//     double t2 = t + dt * 0.5;
//     evaluate_field(ef2,xk,q_ws,t2);
//     // std::vector<double> xk_cpy (xk);
//     // std::vector<double> xk_cpy2 (xk);
//     // (*calculate_e)(ef2.data(), xk_cpy.data(), ef2.size(),
//     //                 xk_cpy2.data(), q_ws.data(), xk.size());
//     // auto stop = high_resolution_clock::now();
//     // add_time(field_time, duration_cast<duration<double>>(stop - start) );

//     for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
//         for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {

//             ef2[xi] *= species_qms[sp_i];

//             // k->3
//             p3[xi] = ps[xi] + 0.5 * dt * ef2[xi];
//             v3[xi] = p3[xi]/ species_ms[sp_i];
//             xk[xi] = xs[xi] + 0.5 * dt * p2[xi];
//         }
//     }




//     // k3 = (v3,a3) = F(un + h/2 k2) = F(xn + delt/2 v2, vn + delt/2 a2)
//     //              = ( vn + delt a2 / 2, q/m E(xn + delt v2 /2) )
//     // for (int ii = 0; ii < N; ++ii) {
//     //     xk[ii] = xs[ii] + 0.5 * dt * p2[ii];
//     // }

//     // start = high_resolution_clock::now();
//     double t3 = t2;
//     evaluate_field(ef3, xk, q_ws, t3);
//     // xk_cpy = xk;
//     // xk_cpy2 = xk;
//     // (*calculate_e)(ef3.data(), xk_cpy.data(), ef3.size(),
//     //                 xk_cpy2.data(), q_ws.data(), xk.size());
//     // stop = high_resolution_clock::now();
//     // add_time(field_time, duration_cast<duration<double>>(stop - start) );
//     // for (int ii = 0; ii < N; ++ii) {
//     //     a3[ii] *= qm;
//     // }
//     for (size_t sp_i = 0; sp_i < N_sp; ++ sp_i) {
//         for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
//             ef3[xi] *= species_qms[sp_i];
//             // k : 3-> 4
//             p4[xi] = ps[xi] + dt * ef3[xi];
//             v4[xi] = p4[xi]/ species_ms[sp_i];
//             xk[xi] = xs[xi] + dt * p3[xi];
//         }
//     }
//     // k4 = (v4,a4) = F(un + h k3) = F(xn + delt v3, pn + delt f3)
//     //              = ( (pn + delt f3) / m, q E(xn + delt v3) )

//     double t4 = t + dt;
//     evaluate_field(ef4, xk, q_ws, t4);
//     // xk_cpy = xk;
//     // xk_cpy2 = xk;
//     // (*calculate_e)(ef4.data(), xk_cpy.data(), ef4.size(), 
//     //                 xk_cpy2.data(), q_ws.data(), xk.size());
//     // stop = high_resolution_clock::now();
//     // add_time(field_time, duration_cast<duration<double>>(stop - start) );

//     // for (int ii = 0; ii < N; ++ii) {
//     //     a4[ii] *= qm;
//     // }
//     for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
//         for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
//             ef4[xi] *= species_qms[sp_i];
//         }
//     }
//     // un+1 = (xn+1,pn+1) = un + h/6 (k1 + 2k2 + 2k3 + k4)
//     //                    = (xn + delt/6 (v1 + 2 v2 + 2 v3 + v4),
//     //                    = pn + delt/6 (a1 + 2 a2 + 2 a3 + a4) )
//     // store xn+1, pn+1 in xn,pn
//     for (int ii = 0; ii < N; ++ii) {
//         xs[ii] += dt / 6.0 * (p1[ii] + 2 * p2[ii] + 2 * p3[ii] + p4[ii]);
//         ps[ii] += dt / 6.0 * (ef1[ii] + 2 * ef2[ii] + 2 * ef3[ii] + ef4[ii]);
//     }
//     if (get_4th_e) {
//         // start = high_resolution_clock::now();
//         double tk = t + dt;
//         evaluate_field(es, xs, q_ws, tk);
//         // xk_cpy = xs;
//         // xk_cpy2 = xs;
//         // (*calculate_e)(es.data(), xk_cpy.data(), xs.size(),
//         //                 xk_cpy2.data(), q_ws.data(), xs.size());
//         // stop = high_resolution_clock::now();
//         // add_time(field_time, duration_cast<duration<double>>(stop - start) );
//     }

//     need_scatter = true;
//     return 0;
// }







