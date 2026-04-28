#include "AMRSimulation.hpp"

// #define DEBUG

int AMRSimulation::step() {
    // ---- Diagnostics & periodic dumps at the START of the step (state at t_n) ----
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }
    {
        MHDDiagnostics d = compute_diagnostics();
        write_diagnostics(d);
    }

    // ---- Advance counters ----
    iter_num += 1;
    general_list[0]->iter_num = iter_num;
    general_list[1]->iter_num = iter_num;
    std::cout << "step " << iter_num << std::endl;

    // ---- Time integration: method==0 -> RK4, else Euler ----
    if (method > 0) {
        euler();
    } else {
        rk4();
    }

    t += dt;

    // ---- Post-step bookkeeping: rebuild tree + recover w0,j0 from q+- + refresh weights ----
    // Both euler() and rk4() leave panels with updated (xs, ys, q0s) but
    // stale w0,j0 and leaf weights. This block makes the panel state consistent
    // for the NEXT step's stage-1 RHS evaluation.
    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();

    // in future amr implementation

    
    return 0;
}


// int AMRSimulation::euler() {
//     cout << "enter euler" << endl;
//     // get information from first panel 
//     u1s.assign(xs.size(), 0.0);
//     u2s.assign(xs.size(), 0.0);
//     b1s.assign(xs.size(), 0.0);
//     b2s.assign(xs.size(), 0.0);
//     std::vector<double> u_ws(xs.size(), 0.0);
//     std::vector<double> b_ws(xs.size(), 0.0);
//     u_ws = general_list[0]->get_u_weights();
//     b_ws = general_list[0]->get_b_weights();
//     if(iter_num >= 1) {
//         xs = general_list[0]->get_xs();
//         ys = general_list[0]->get_ys();

//         general_list[0]->evaluate_u_field(u1s, u2s, xs, ys, u_ws, t);
//         general_list[0]->evaluate_b_field(b1s, b2s, xs, ys, b_ws, t);

//         //external field for alfven wave
//         // for (size_t i = 0; i < b1s.size(); ++i) {
//         //     b1s[i] +=  1.0;
//         // }

//         // external field for polarized_alfven wave
//         // for (size_t i = 0; i < b1s.size(); ++i) {
//         //     b1s[i] +=  2.0/sqrt(5);
//         //     b2s[i] +=  1.0/sqrt(5);
//         // }
        
//         general_list[0]->set_u1s(u1s);
//         general_list[0]->set_u2s(u2s);
//         general_list[0]->set_b1s(b1s);
//         general_list[0]->set_b2s(b2s);

//         general_list[1]->set_u1s(u1s);
//         general_list[1]->set_u2s(u2s);
//         general_list[1]->set_b1s(b1s);
//         general_list[1]->set_b2s(b2s);
//     }
//     else{
//         xs = general_list[0]->get_xs();
//         ys = general_list[0]->get_ys();
//         u1s = general_list[0]->get_u1s();
//         u2s = general_list[0]->get_u2s();
//         b1s = general_list[0]->get_b1s();
//         b2s = general_list[0]->get_b2s();
//     }


//     // source terms using finte difference 
//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     std::vector<double> source_term(xs.size(), 0.0);
//     // general_list[0]->compute_rhs_state(xs, ys, u1s, u2s, b1s, b2s, t, 
//     //             u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
//     //             b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, source_term);


//     // compute source term using integral 
//     if (bcs == periodic_bcs) {
//         general_list[0]->evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs, ys, u_ws, t);
//         general_list[0]->evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs, ys, u_ws, t);
//         general_list[0]->evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs, ys, b_ws, t);
//         general_list[0]->evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs, ys, b_ws, t);
//     }
//     else { // future change to open y 
//         general_list[0]->evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs, ys, u_ws, t);
//         general_list[0]->evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs, ys, u_ws, t);
//         general_list[0]->evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs, ys, b_ws, t);
//         general_list[0]->evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs, ys, b_ws, t);
//     }

//     for (int i = 0; i < xs.size(); ++i) {
//         source_term[i] = 2*(b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i])
//                         -2* (b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i]);
//     }



//     // build advection velocities
//     std::vector<double> vx_plus(xs.size()), vy_plus(xs.size());
//     std::vector<double> vx_minus(xs.size()), vy_minus(xs.size());

//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         vx_plus[ii]  = u1s[ii] - b1s[ii];
//         vy_plus[ii]  = u2s[ii] - b2s[ii];
//         vx_minus[ii] = u1s[ii] + b1s[ii];
//         vy_minus[ii] = u2s[ii] + b2s[ii];
//     }
//     // q+ panel
//     general_list[0]->update_panel_positions_and_q(vx_plus, vy_plus, source_term, dt, +1);
//     // q- panel
//     general_list[1]->update_panel_positions_and_q(vx_minus, vy_minus, source_term, dt, -1);

//     return 0;
// }


int AMRSimulation::euler() {
    cout << "enter euler" << endl;

    // Cache originals (so we push from the state at t_n, consistent with rk4())
    std::vector<double> orig_xp = general_list[0]->get_xs();
    std::vector<double> orig_yp = general_list[0]->get_ys();
    std::vector<double> orig_qp = general_list[0]->get_q0s();

    std::vector<double> orig_xm = general_list[1]->get_xs();
    std::vector<double> orig_ym = general_list[1]->get_ys();
    std::vector<double> orig_qm = general_list[1]->get_q0s();

    const size_t N = orig_xp.size();

    // ============================================================
    // STAGE 1 (the only stage for forward Euler) : RHS at (X_n, t_n)
    // Panels currently hold X_n with valid weights from the previous step()'s
    // post-block (or from the initial constructor for iter 0).
    // ============================================================
    std::vector<double> u1, u2, b1, b2, S;
    compute_stage_rhs(t, orig_xp, orig_yp, u1, u2, b1, b2, S);

    // Cache field values on the panels (needed by compute_diagnostics() in step())
    general_list[0]->set_u1s(u1);
    general_list[0]->set_u2s(u2);
    general_list[0]->set_b1s(b1);
    general_list[0]->set_b2s(b2);
    general_list[1]->set_u1s(u1);
    general_list[1]->set_u2s(u2);
    general_list[1]->set_b1s(b1);
    general_list[1]->set_b2s(b2);

    // Build advection velocities for q+ and q- panels
    std::vector<double> vxp(N), vyp(N), vxm(N), vym(N);
    for (size_t i = 0; i < N; ++i) {
        vxp[i] = u1[i] - b1[i];
        vyp[i] = u2[i] - b2[i];
        vxm[i] = u1[i] + b1[i];
        vym[i] = u2[i] + b2[i];
    }

    // Forward-Euler update from the cached originals (NOT in-place)
    for (size_t i = 0; i < N; ++i) {
        orig_xp[i] += dt * vxp[i];
        orig_yp[i] += dt * vyp[i];
        orig_qp[i] += dt * S[i];

        orig_xm[i] += dt * vxm[i];
        orig_ym[i] += dt * vym[i];
        orig_qm[i] -= dt * S[i];
    }

    general_list[0]->set_xs(orig_xp);
    general_list[0]->set_ys(orig_yp);
    general_list[0]->set_q0s(orig_qp);

    general_list[1]->set_xs(orig_xm);
    general_list[1]->set_ys(orig_ym);
    general_list[1]->set_q0s(orig_qm);

    // remesh + recover + leaves_weights happens in step() (same as rk4())
    return 0;
}






void AMRSimulation::compute_stage_rhs(double t_stage,
        std::vector<double>& stage_xs,
        std::vector<double>& stage_ys,
        std::vector<double>& stage_u1s,
        std::vector<double>& stage_u2s,
        std::vector<double>& stage_b1s,
        std::vector<double>& stage_b2s,
        std::vector<double>& stage_source)
{
    const size_t N = stage_xs.size();
    stage_u1s.assign(N, 0.0);
    stage_u2s.assign(N, 0.0);
    stage_b1s.assign(N, 0.0);
    stage_b2s.assign(N, 0.0);
    stage_source.assign(N, 0.0);

    // Weights from the CURRENT panel state.
    // Caller is responsible for ensuring the panels have been
    // set + remeshed + recovered + set_leaves_weights() before this call.
    std::vector<double> u_ws = general_list[0]->get_u_weights();
    std::vector<double> b_ws = general_list[0]->get_b_weights();

    // (1) Fields at stage points
    general_list[0]->evaluate_u_field(stage_u1s, stage_u2s,
                                      stage_xs, stage_ys, u_ws, t_stage);
    general_list[0]->evaluate_b_field(stage_b1s, stage_b2s,
                                      stage_xs, stage_ys, b_ws, t_stage);

    for (size_t i = 0; i < N; ++i) {
        stage_b1s[i] += B0x;
        stage_b2s[i] += B0y;
    }

    general_list[0]->set_u1s(stage_u1s);
    general_list[0]->set_u2s(stage_u2s);
    general_list[0]->set_b1s(stage_b1s);
    general_list[0]->set_b2s(stage_b2s);


    // (Optional external field — keep identical to euler() for consistency)
    // for (size_t i = 0; i < stage_b1s.size(); ++i) stage_b1s[i] += 1.0;            // Alfvén
    // for (size_t i = 0; i < stage_b1s.size(); ++i) {                                // polarized Alfvén
    //     stage_b1s[i] += 2.0/std::sqrt(5.0);
    //     stage_b2s[i] += 1.0/std::sqrt(5.0);
    // }

    // (2) Gradients via INTEGRAL formulation
    std::vector<double> du1dx(N, 0.0), du1dy(N, 0.0);
    std::vector<double> du2dx(N, 0.0), du2dy(N, 0.0);
    std::vector<double> db1dx(N, 0.0), db1dy(N, 0.0);
    std::vector<double> db2dx(N, 0.0), db2dy(N, 0.0);

    if (bcs == periodic_bcs) {
        general_list[0]->evaluate_u1s_grad(du1dx, du1dy, stage_xs, stage_ys, u_ws, t_stage);
        general_list[0]->evaluate_u2s_grad(du2dx, du2dy, stage_xs, stage_ys, u_ws, t_stage);
        general_list[0]->evaluate_b1s_grad(db1dx, db1dy, stage_xs, stage_ys, b_ws, t_stage);
        general_list[0]->evaluate_b2s_grad(db2dx, db2dy, stage_xs, stage_ys, b_ws, t_stage);
    } else {
        general_list[0]->evaluate_u1s_grad(du1dx, du1dy, stage_xs, stage_ys, u_ws, t_stage);
        general_list[0]->evaluate_u2s_grad(du2dx, du2dy, stage_xs, stage_ys, u_ws, t_stage);
        general_list[0]->evaluate_b1s_grad(db1dx, db1dy, stage_xs, stage_ys, b_ws, t_stage);
        general_list[0]->evaluate_b2s_grad(db2dx, db2dy, stage_xs, stage_ys, b_ws, t_stage);
    }

    // (3) Source term — exact formula from your euler()
    for (size_t i = 0; i < N; ++i) {
        stage_source[i] = 2.0 * (db1dx[i]*du2dx[i] + db2dx[i]*du2dy[i])
                        - 2.0 * (db1dy[i]*du1dx[i] + db2dy[i]*du1dy[i]);
    }
}





int AMRSimulation::rk4() {
    cout << "enter rk4" << endl;

    // ---------- Cache the original state of both panels at t_n ----------
    std::vector<double> orig_xp = general_list[0]->get_xs();
    std::vector<double> orig_yp = general_list[0]->get_ys();
    std::vector<double> orig_qp = general_list[0]->get_q0s();

    std::vector<double> orig_xm = general_list[1]->get_xs();
    std::vector<double> orig_ym = general_list[1]->get_ys();
    std::vector<double> orig_qm = general_list[1]->get_q0s();

    const size_t N = orig_xp.size();

    // ============================================================
    // STAGE 1 : RHS at (X_n, t_n)
    // Panels currently hold X_n with valid weights, so just evaluate.
    // ============================================================
    std::vector<double> k1_u1, k1_u2, k1_b1, k1_b2, k1_S;
    compute_stage_rhs(t, orig_xp, orig_yp,
                      k1_u1, k1_u2, k1_b1, k1_b2, k1_S);

    std::vector<double> k1_vxp(N), k1_vyp(N), k1_vxm(N), k1_vym(N);
    for (size_t i = 0; i < N; ++i) {
        k1_vxp[i] = k1_u1[i] - k1_b1[i];
        k1_vyp[i] = k1_u2[i] - k1_b2[i];
        k1_vxm[i] = k1_u1[i] + k1_b1[i];
        k1_vym[i] = k1_u2[i] + k1_b2[i];
    }

    // ---------- Build stage-2 probe state X_n + (dt/2) k1 ----------
    std::vector<double> x2p(N), y2p(N), q2p(N);
    std::vector<double> x2m(N), y2m(N), q2m(N);
    for (size_t i = 0; i < N; ++i) {
        x2p[i] = orig_xp[i] + 0.5 * dt * k1_vxp[i];
        y2p[i] = orig_yp[i] + 0.5 * dt * k1_vyp[i];
        q2p[i] = orig_qp[i] + 0.5 * dt * k1_S[i];

        x2m[i] = orig_xm[i] + 0.5 * dt * k1_vxm[i];
        y2m[i] = orig_ym[i] + 0.5 * dt * k1_vym[i];
        q2m[i] = orig_qm[i] - 0.5 * dt * k1_S[i];
    }
    general_list[0]->set_xs(x2p);
    general_list[0]->set_ys(y2p);
    general_list[0]->set_q0s(q2p);
    general_list[1]->set_xs(x2m);
    general_list[1]->set_ys(y2m);
    general_list[1]->set_q0s(q2m);

    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();

    // ============================================================
    // STAGE 2 : RHS at (X_n + dt/2 k1, t_n + dt/2)
    // ============================================================
    std::vector<double> k2_u1, k2_u2, k2_b1, k2_b2, k2_S;
    {
        std::vector<double> sxp = general_list[0]->get_xs();
        std::vector<double> syp = general_list[0]->get_ys();
        compute_stage_rhs(t + 0.5 * dt, sxp, syp,
                          k2_u1, k2_u2, k2_b1, k2_b2, k2_S);
    }

    std::vector<double> k2_vxp(N), k2_vyp(N), k2_vxm(N), k2_vym(N);
    for (size_t i = 0; i < N; ++i) {
        k2_vxp[i] = k2_u1[i] - k2_b1[i];
        k2_vyp[i] = k2_u2[i] - k2_b2[i];
        k2_vxm[i] = k2_u1[i] + k2_b1[i];
        k2_vym[i] = k2_u2[i] + k2_b2[i];
    }

    // ---------- Build stage-3 probe state X_n + (dt/2) k2 ----------
    std::vector<double> x3p(N), y3p(N), q3p(N);
    std::vector<double> x3m(N), y3m(N), q3m(N);
    for (size_t i = 0; i < N; ++i) {
        x3p[i] = orig_xp[i] + 0.5 * dt * k2_vxp[i];
        y3p[i] = orig_yp[i] + 0.5 * dt * k2_vyp[i];
        q3p[i] = orig_qp[i] + 0.5 * dt * k2_S[i];

        x3m[i] = orig_xm[i] + 0.5 * dt * k2_vxm[i];
        y3m[i] = orig_ym[i] + 0.5 * dt * k2_vym[i];
        q3m[i] = orig_qm[i] - 0.5 * dt * k2_S[i];
    }
    general_list[0]->set_xs(x3p);
    general_list[0]->set_ys(y3p);
    general_list[0]->set_q0s(q3p);
    general_list[1]->set_xs(x3m);
    general_list[1]->set_ys(y3m);
    general_list[1]->set_q0s(q3m);

    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();

    // ============================================================
    // STAGE 3 : RHS at (X_n + dt/2 k2, t_n + dt/2)
    // ============================================================
    std::vector<double> k3_u1, k3_u2, k3_b1, k3_b2, k3_S;
    {
        std::vector<double> sxp = general_list[0]->get_xs();
        std::vector<double> syp = general_list[0]->get_ys();
        compute_stage_rhs(t + 0.5 * dt, sxp, syp,
                          k3_u1, k3_u2, k3_b1, k3_b2, k3_S);
    }

    std::vector<double> k3_vxp(N), k3_vyp(N), k3_vxm(N), k3_vym(N);
    for (size_t i = 0; i < N; ++i) {
        k3_vxp[i] = k3_u1[i] - k3_b1[i];
        k3_vyp[i] = k3_u2[i] - k3_b2[i];
        k3_vxm[i] = k3_u1[i] + k3_b1[i];
        k3_vym[i] = k3_u2[i] + k3_b2[i];
    }

    // ---------- Build stage-4 probe state X_n + dt k3 (full step!) ----------
    std::vector<double> x4p(N), y4p(N), q4p(N);
    std::vector<double> x4m(N), y4m(N), q4m(N);
    for (size_t i = 0; i < N; ++i) {
        x4p[i] = orig_xp[i] + dt * k3_vxp[i];
        y4p[i] = orig_yp[i] + dt * k3_vyp[i];
        q4p[i] = orig_qp[i] + dt * k3_S[i];

        x4m[i] = orig_xm[i] + dt * k3_vxm[i];
        y4m[i] = orig_ym[i] + dt * k3_vym[i];
        q4m[i] = orig_qm[i] - dt * k3_S[i];
    }
    general_list[0]->set_xs(x4p);
    general_list[0]->set_ys(y4p);
    general_list[0]->set_q0s(q4p);
    general_list[1]->set_xs(x4m);
    general_list[1]->set_ys(y4m);
    general_list[1]->set_q0s(q4m);

    general_list[0]->remesh();
    general_list[1]->remesh();
    general_list[0]->recover_wj_from_q(general_list[1]);
    general_list[1]->copy_wj_from(general_list[0]);
    general_list[0]->set_leaves_weights();
    general_list[1]->set_leaves_weights();

    // ============================================================
    // STAGE 4 : RHS at (X_n + dt k3, t_n + dt)
    // ============================================================
    std::vector<double> k4_u1, k4_u2, k4_b1, k4_b2, k4_S;
    {
        std::vector<double> sxp = general_list[0]->get_xs();
        std::vector<double> syp = general_list[0]->get_ys();
        compute_stage_rhs(t + dt, sxp, syp,
                          k4_u1, k4_u2, k4_b1, k4_b2, k4_S);
    }

    std::vector<double> k4_vxp(N), k4_vyp(N), k4_vxm(N), k4_vym(N);
    for (size_t i = 0; i < N; ++i) {
        k4_vxp[i] = k4_u1[i] - k4_b1[i];
        k4_vyp[i] = k4_u2[i] - k4_b2[i];
        k4_vxm[i] = k4_u1[i] + k4_b1[i];
        k4_vym[i] = k4_u2[i] + k4_b2[i];
    }

    // ============================================================
    // FINAL COMBINE : X_{n+1} = X_n + (dt/6)(k1 + 2 k2 + 2 k3 + k4)
    // ============================================================
    for (size_t i = 0; i < N; ++i) {
        orig_xp[i] += (dt / 6.0) * (k1_vxp[i] + 2.0*k2_vxp[i] + 2.0*k3_vxp[i] + k4_vxp[i]);
        orig_yp[i] += (dt / 6.0) * (k1_vyp[i] + 2.0*k2_vyp[i] + 2.0*k3_vyp[i] + k4_vyp[i]);
        orig_qp[i] += (dt / 6.0) * (k1_S[i]   + 2.0*k2_S[i]   + 2.0*k3_S[i]   + k4_S[i]);

        orig_xm[i] += (dt / 6.0) * (k1_vxm[i] + 2.0*k2_vxm[i] + 2.0*k3_vxm[i] + k4_vxm[i]);
        orig_ym[i] += (dt / 6.0) * (k1_vym[i] + 2.0*k2_vym[i] + 2.0*k3_vym[i] + k4_vym[i]);
        orig_qm[i] -= (dt / 6.0) * (k1_S[i]   + 2.0*k2_S[i]   + 2.0*k3_S[i]   + k4_S[i]);
    }

    general_list[0]->set_xs(orig_xp);
    general_list[0]->set_ys(orig_yp);
    general_list[0]->set_q0s(orig_qp);
    general_list[1]->set_xs(orig_xm);
    general_list[1]->set_ys(orig_ym);
    general_list[1]->set_q0s(orig_qm);

    // remesh + recover_wj_from_q + set_leaves_weights happens in step()
    return 0;
}














//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------- Finite Difference ----------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//

// int AMRSimulation::rk4() {
//     cout << "enter rk4" << endl;

//     xs = general_list[0]->get_xs();
//     ys = general_list[0]->get_ys();

//     // copy xs,ysq0s from two panels so that even we change panel, we will push based on their original values
//     std::vector<double> original_q_plus_xs(xs.size()), original_q_plus_ys(xs.size()), original_q_plus_q0s;
//     std::vector<double> original_q_minus_xs(xs.size()), original_q_minus_ys(xs.size()), original_q_minus_q0s;
//     original_q_plus_xs = general_list[0]->get_xs();
//     original_q_plus_ys = general_list[0]->get_ys();
//     original_q_plus_q0s = general_list[0]->get_q0s();

//     original_q_minus_xs = general_list[1]->get_xs();
//     original_q_minus_ys = general_list[1]->get_ys();
//     original_q_minus_q0s = general_list[1]->get_q0s();




//     const int xs_size = (int)xs.size();

//     // Stage storage: k1, k2, k3, k4 for x, y, w and j
//     std::vector<double> k1_u1s(xs_size, 0.0), k1_u2s(xs_size, 0.0), k1_b1s(xs_size, 0.0), k1_b2s(xs_size, 0.0), k1_source_term(xs.size(), 0.0);
//     std::vector<double> k2_u1s(xs_size, 0.0), k2_u2s(xs_size, 0.0), k2_b1s(xs_size, 0.0), k2_b2s(xs_size, 0.0), k2_source_term(xs.size(), 0.0);
//     std::vector<double> k3_u1s(xs_size, 0.0), k3_u2s(xs_size, 0.0), k3_b1s(xs_size, 0.0), k3_b2s(xs_size, 0.0), k3_source_term(xs.size(), 0.0);
//     std::vector<double> k4_u1s(xs_size, 0.0), k4_u2s(xs_size, 0.0), k4_b1s(xs_size, 0.0), k4_b2s(xs_size, 0.0), k4_source_term(xs.size(), 0.0);

//     // Temporary positions
//     std::vector<double> k2_xs_plus(xs_size), k2_ys_plus(xs_size), k2_q_plus(xs_size), k2_xs_minus(xs_size), k2_ys_minus(xs_size), k2_q_minus(xs_size);
//     std::vector<double> k3_xs_plus(xs_size), k3_ys_plus(xs_size), k3_q_plus(xs_size), k3_xs_minus(xs_size), k3_ys_minus(xs_size), k3_q_minus(xs_size);
//     std::vector<double> k4_xs_plus(xs_size), k4_ys_plus(xs_size), k4_q_plus(xs_size), k4_xs_minus(xs_size), k4_ys_minus(xs_size), k4_q_minus(xs_size);


//     std::vector<double> k1_u_ws(xs.size(), 0.0), k2_u_ws(xs.size(), 0.0), k3_u_ws(xs.size(), 0.0), k4_u_ws(xs.size(), 0.0);
//     std::vector<double> k1_b_ws(xs.size(), 0.0), k2_b_ws(xs.size(), 0.0), k3_b_ws(xs.size(), 0.0), k4_b_ws(xs.size(), 0.0);

//     if(iter_num > 1) {
//         k1_u_ws = general_list[0]->get_u_weights();
//         k1_b_ws = general_list[0]->get_b_weights();
//         general_list[0]->evaluate_u_field(k1_u1s, k1_u2s, xs, ys, k1_u_ws, t);
//         general_list[0]->evaluate_b_field(k1_b1s, k1_b2s, xs, ys, k1_b_ws, t);

//         //external field for alfven wave
//         // for (size_t i = 0; i < k1_b1s.size(); ++i) {
//         //     k1_b1s[i] +=  1.0;
//         // }

//         // external field for polarized_alfven wave
//         // for (size_t i = 0; i < k1_b1s.size(); ++i) {
//         //     k1_b1s[i] +=  2.0/sqrt(5);
//         //     k1_b2s[i] +=  1.0/sqrt(5);
//         // }


//         general_list[0]->set_u1s(k1_u1s);
//         general_list[0]->set_u2s(k1_u2s);
//         general_list[0]->set_b1s(k1_b1s);
//         general_list[0]->set_b2s(k1_b2s);

//         general_list[1]->set_u1s(k1_u1s);
//         general_list[1]->set_u2s(k1_u2s);
//         general_list[1]->set_b1s(k1_b1s);
//         general_list[1]->set_b2s(k1_b2s);
//     }
//     else{
//         k1_u1s = general_list[0]->get_u1s();
//         k1_u2s = general_list[0]->get_u2s();
//         k1_b1s = general_list[0]->get_b1s();
//         k1_b2s = general_list[0]->get_b2s();
//     }

//     // ------------------------------------------------------------
//     // Stage 1 : k1 = (k1_u1s, k1_u2s, k1_b1s, k1_b2s, k1_source_term)
//     // ------------------------------------------------------------

//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     // source term to update q0
//     general_list[0]->compute_rhs_state(xs, ys, k1_u1s, k1_u2s, k1_b1s, k1_b2s, t, 
//                 u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
//                 b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k1_source_term);


//     // build advection velocities
//     std::vector<double> k1_vx_plus(xs.size()), k1_vy_plus(xs.size());
//     std::vector<double> k1_vx_minus(xs.size()), k1_vy_minus(xs.size());

//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k1_vx_plus[ii]  = k1_u1s[ii] - k1_b1s[ii];
//         k1_vy_plus[ii]  = k1_u2s[ii] - k1_b2s[ii];
//         k1_vx_minus[ii] = k1_u1s[ii] + k1_b1s[ii];
//         k1_vy_minus[ii] = k1_u2s[ii] + k1_b2s[ii];
//     }

//     // q+ panel 
//     // double k1_dt = 0.5 * dt;
//     // general_list[0]->update_panel_positions_and_q(k1_vx_plus, k1_vy_plus, k1_source_term, k1_dt, +1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k2_xs_plus[ii] = original_q_plus_xs[ii] + 0.5 * dt * k1_vx_plus[ii]; 
//         k2_ys_plus[ii] = original_q_plus_ys[ii] + 0.5 * dt * k1_vy_plus[ii]; 
//         k2_q_plus[ii] = original_q_plus_q0s[ii] + 0.5 * dt * k1_source_term[ii]; 
//     }

//     // q- panel
//     // general_list[1]->update_panel_positions_and_q(k1_vx_minus, k1_vy_minus, k1_source_term, k1_dt, -1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k2_xs_minus[ii] = original_q_minus_xs[ii] + 0.5 * dt * k1_vx_minus[ii]; 
//         k2_ys_minus[ii] = original_q_minus_ys[ii] + 0.5 * dt * k1_vy_minus[ii]; 
//         k2_q_minus[ii] = original_q_minus_q0s[ii] - 0.5 * dt * k1_source_term[ii]; 
//     }

//     general_list[0]->set_xs(k2_xs_plus);
//     general_list[0]->set_ys(k2_ys_plus);
//     general_list[0]->set_q0s(k2_q_plus);

//     general_list[1]->set_xs(k2_xs_minus);
//     general_list[1]->set_ys(k2_ys_minus);
//     general_list[1]->set_q0s(k2_q_minus);



//     // remesh first to calculate w0 and j0, then reset leave weights for future u and b evaluations
//     general_list[0]->remesh();
//     general_list[1]->remesh();
//     general_list[0]->recover_wj_from_q(general_list[1]);
//     general_list[1]->copy_wj_from(general_list[0]);
//     general_list[0]->set_leaves_weights();
//     general_list[1]->set_leaves_weights();


//     // ------------------------------------------------------------
//     // Stage 2 : k2 = (k2_u1s, k2_u2s, k2_b1s, k2_b2s, k2_source_term)
//     // ------------------------------------------------------------
//     xs = general_list[0]->get_xs();
//     ys = general_list[0]->get_ys();
//     k2_u_ws = general_list[0]->get_u_weights();
//     k2_b_ws = general_list[0]->get_b_weights();
//     general_list[0]->evaluate_u_field(k2_u1s, k2_u2s, xs, ys, k2_u_ws, t + 0.5 * dt);
//     general_list[0]->evaluate_b_field(k2_b1s, k2_b2s, xs, ys, k2_b_ws, t + 0.5 * dt);

//     //external field for alfven wave
//     // for (size_t i = 0; i < k2_b1s.size(); ++i) {
//     //     k2_b1s[i] +=  1.0;
//     // }

//     //external field for polarized_alfven wave
//     // for (size_t i = 0; i < k2_b1s.size(); ++i) {
//     //     k2_b1s[i] +=  2.0/sqrt(5);
//     //     k2_b2s[i] +=  1.0/sqrt(5);
//     // }

//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     // source term to update q0
//     general_list[0]->compute_rhs_state(xs, ys, k2_u1s, k2_u2s, k2_b1s, k2_b2s, t + 0.5 * dt, 
//                 u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
//                 b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k2_source_term);

//     // build advection velocities
//     std::vector<double> k2_vx_plus(xs.size()), k2_vy_plus(xs.size());
//     std::vector<double> k2_vx_minus(xs.size()), k2_vy_minus(xs.size());

//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k2_vx_plus[ii]  = k2_u1s[ii] - k2_b1s[ii];
//         k2_vy_plus[ii]  = k2_u2s[ii] - k2_b2s[ii];
//         k2_vx_minus[ii] = k2_u1s[ii] + k2_b1s[ii];
//         k2_vy_minus[ii] = k2_u2s[ii] + k2_b2s[ii];
//     }
//     // double k2_dt = 0.5 * dt;
//     // q+ panel 
//     // general_list[0]->update_panel_positions_and_q(k2_vx_plus, k2_vy_plus, k2_source_term, k2_dt, +1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k3_xs_plus[ii] = original_q_plus_xs[ii] + 0.5 * dt * k2_vx_plus[ii]; 
//         k3_ys_plus[ii] = original_q_plus_ys[ii] + 0.5 * dt * k2_vy_plus[ii]; 
//         k3_q_plus[ii] = original_q_plus_q0s[ii] + 0.5 * dt * k2_source_term[ii]; 
//     }
//     // q- panel
//     // general_list[1]->update_panel_positions_and_q(k2_vx_minus, k2_vy_minus, k2_source_term, k2_dt, -1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k3_xs_minus[ii] = original_q_minus_xs[ii] + 0.5 * dt * k2_vx_minus[ii]; 
//         k3_ys_minus[ii] = original_q_minus_ys[ii] + 0.5 * dt * k2_vy_minus[ii]; 
//         k3_q_minus[ii] = original_q_minus_q0s[ii] - 0.5 * dt * k2_source_term[ii]; 
//     }

//     general_list[0]->set_xs(k3_xs_plus);
//     general_list[0]->set_ys(k3_ys_plus);
//     general_list[0]->set_q0s(k3_q_plus);

//     general_list[1]->set_xs(k3_xs_minus);
//     general_list[1]->set_ys(k3_ys_minus);
//     general_list[1]->set_q0s(k3_q_minus);

//     general_list[0]->remesh();
//     general_list[1]->remesh();
//     general_list[0]->recover_wj_from_q(general_list[1]);
//     general_list[1]->copy_wj_from(general_list[0]);
//     general_list[0]->set_leaves_weights();
//     general_list[1]->set_leaves_weights();



//     // ------------------------------------------------------------
//     // Stage 3 : k3 = (k3_u1s, k3_u2s, k3_b1s, k3_b2s, k3_source_term)
//     // ------------------------------------------------------------
//     xs = general_list[0]->get_xs();
//     ys = general_list[0]->get_ys();
//     k3_u_ws = general_list[0]->get_u_weights();
//     k3_b_ws = general_list[0]->get_b_weights();
//     general_list[0]->evaluate_u_field(k3_u1s, k3_u2s, xs, ys, k3_u_ws, t + 0.5 * dt);
//     general_list[0]->evaluate_b_field(k3_b1s, k3_b2s, xs, ys, k3_b_ws, t + 0.5 * dt);

//     //external field for alfven wave
//     // for (size_t i = 0; i < k3_b1s.size(); ++i) {
//     //     k3_b1s[i] +=  1.0;
//     // }

//     //external field for polarized_alfven wave
//     // for (size_t i = 0; i < k3_b1s.size(); ++i) {
//     //     k3_b1s[i] +=  2/sqrt(5);
//     //     k3_b2s[i] +=  1/sqrt(5);
//     // }

//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     // source term to update q0
//     general_list[0]->compute_rhs_state(xs, ys, k3_u1s, k3_u2s, k3_b1s, k3_b2s, t + 0.5 * dt, 
//                 u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
//                 b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k3_source_term);

//     // build advection velocities
//     std::vector<double> k3_vx_plus(xs.size()), k3_vy_plus(xs.size());
//     std::vector<double> k3_vx_minus(xs.size()), k3_vy_minus(xs.size());

//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k3_vx_plus[ii]  = k3_u1s[ii] - k3_b1s[ii];
//         k3_vy_plus[ii]  = k3_u2s[ii] - k3_b2s[ii];
//         k3_vx_minus[ii] = k3_u1s[ii] + k3_b1s[ii];
//         k3_vy_minus[ii] = k3_u2s[ii] + k3_b2s[ii];
//     }
//     // double k3_dt = dt;
//     // q+ panel 
//     // general_list[0]->update_panel_positions_and_q(k3_vx_plus, k3_vy_plus, k3_source_term, k3_dt, +1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k4_xs_plus[ii] = original_q_plus_xs[ii] + dt * k3_vx_plus[ii]; 
//         k4_ys_plus[ii] = original_q_plus_ys[ii] + dt * k3_vy_plus[ii]; 
//         k4_q_plus[ii] = original_q_plus_q0s[ii] + dt * k3_source_term[ii]; 
//     }
//     // q- panel
//     // general_list[1]->update_panel_positions_and_q(k3_vx_minus, k3_vy_minus, k3_source_term, k3_dt, -1);
//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k4_xs_minus[ii] = original_q_minus_xs[ii] + dt * k3_vx_minus[ii]; 
//         k4_ys_minus[ii] = original_q_minus_ys[ii] + dt * k3_vy_minus[ii]; 
//         k4_q_minus[ii] = original_q_minus_q0s[ii] - dt * k3_source_term[ii]; 
//     }

//     general_list[0]->set_xs(k4_xs_plus);
//     general_list[0]->set_ys(k4_ys_plus);
//     general_list[0]->set_q0s(k4_q_plus);

//     general_list[1]->set_xs(k4_xs_minus);
//     general_list[1]->set_ys(k4_ys_minus);
//     general_list[1]->set_q0s(k4_q_minus);

//     general_list[0]->remesh();
//     general_list[1]->remesh();
//     general_list[0]->recover_wj_from_q(general_list[1]);
//     general_list[1]->copy_wj_from(general_list[0]);
//     general_list[0]->set_leaves_weights();
//     general_list[1]->set_leaves_weights();


//     // ------------------------------------------------------------
//     // Stage 4 : k4 = (k4_u1s, k4_u2s, k4_b1s, k4_b2s, k4_source_term)
//     // ------------------------------------------------------------

//     xs = general_list[0]->get_xs();
//     ys = general_list[0]->get_ys();
//     k4_u_ws = general_list[0]->get_u_weights();
//     k4_b_ws = general_list[0]->get_b_weights();
//     general_list[0]->evaluate_u_field(k4_u1s, k4_u2s, xs, ys, k4_u_ws, t + dt);
//     general_list[0]->evaluate_b_field(k4_b1s, k4_b2s, xs, ys, k4_b_ws, t + dt);

//     //external field for alfven wave
//     // for (size_t i = 0; i < k4_b1s.size(); ++i) {
//     //     k4_b1s[i] +=  1.0;
//     // }

//     //external field for polarized_alfven wave
//     // for (size_t i = 0; i < k4_b1s.size(); ++i) {
//     //     k4_b1s[i] +=  2.0/sqrt(5);
//     //     k4_b2s[i] +=  1.0/sqrt(5);
//     // }

//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     // source term to update q0
//     general_list[0]->compute_rhs_state(xs, ys, k4_u1s, k4_u2s, k4_b1s, k4_b2s, t + dt, 
//                 u1s_grad_x, u1s_grad_y, u2s_grad_x, u2s_grad_y, 
//                 b1s_grad_x, b1s_grad_y, b2s_grad_x, b2s_grad_y, k4_source_term);

//     // build advection velocities
//     std::vector<double> k4_vx_plus(xs.size()), k4_vy_plus(xs.size());
//     std::vector<double> k4_vx_minus(xs.size()), k4_vy_minus(xs.size());

//     for (size_t ii = 0; ii < xs.size(); ++ii) {
//         k4_vx_plus[ii]  = k4_u1s[ii] - k4_b1s[ii];
//         k4_vy_plus[ii]  = k4_u2s[ii] - k4_b2s[ii];
//         k4_vx_minus[ii] = k4_u1s[ii] + k4_b1s[ii];
//         k4_vy_minus[ii] = k4_u2s[ii] + k4_b2s[ii];
//     }




//     for (int i = 0; i < xs_size; ++i) {
//         original_q_plus_xs[i] += (dt / 6.0) * (k1_vx_plus[i] + 2.0 * k2_vx_plus[i] + 2.0 * k3_vx_plus[i] + k4_vx_plus[i]);
//         original_q_plus_ys[i] += (dt / 6.0) * (k1_vy_plus[i] + 2.0 * k2_vy_plus[i] + 2.0 * k3_vy_plus[i] + k4_vy_plus[i]);
//         original_q_plus_q0s[i] += (dt / 6.0) * (k1_source_term[i] + 2.0 * k2_source_term[i] + 2.0 * k3_source_term[i] + k4_source_term[i]);
        
//         original_q_minus_xs[i] += (dt / 6.0) * (k1_vx_minus[i] + 2.0 * k2_vx_minus[i] + 2.0 * k3_vx_minus[i] + k4_vx_minus[i]);
//         original_q_minus_ys[i] += (dt / 6.0) * (k1_vy_minus[i] + 2.0 * k2_vy_minus[i] + 2.0 * k3_vy_minus[i] + k4_vy_minus[i]);
//         original_q_minus_q0s[i] -= (dt / 6.0) * (k1_source_term[i] + 2.0 * k2_source_term[i] + 2.0 * k3_source_term[i] + k4_source_term[i]);
//     }

//     // copy back to 2 panels and then reconstruct based on this updation
//     general_list[0]->set_xs(original_q_plus_xs);
//     general_list[0]->set_ys(original_q_plus_ys);
//     general_list[0]->set_q0s(original_q_plus_q0s);

//     general_list[1]->set_xs(original_q_minus_xs);
//     general_list[1]->set_ys(original_q_minus_ys);
//     general_list[1]->set_q0s(original_q_minus_q0s);

//     // out of rk4, in step:
//     // general_list[0]->remesh();
//     // general_list[1]->remesh();
//     // general_list[0]->recover_wj_from_q(general_list[1]);
//     // general_list[1]->copy_wj_from(general_list[0]);
//     // general_list[0]->set_leaves_weights();
//     // general_list[1]->set_leaves_weights();
//     return 0;
// }








