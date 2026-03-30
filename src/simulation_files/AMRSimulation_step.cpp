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

    // rk4 step
    // rk4_step(false);

    // euler 
    euler();
    
    t += dt;

    // must remesh to calculate w0 and j0 based on q+ and q-
    general_list[0]->remesh();
    general_list[1]->remesh();
    // if (iter_num % n_steps_remesh == 0) {
    //     remesh();
    // }

    // after remesh, calculate w0 and j0
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);

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
    xs = general_list[0]->get_xs();
    ys = general_list[0]->get_ys();
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    std::vector<double> u_ws(xs.size(), 0.0);
    std::vector<double> b_ws(xs.size(), 0.0);
    u_ws = general_list[0]->get_u_weights();
    b_ws = general_list[0]->get_b_weights();
    general_list[0]->evaluate_u_field(u1s, u2s, xs, ys, u_ws, t);
    general_list[0]->evaluate_b_field(b1s, b2s, xs, ys, b_ws, t);

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







