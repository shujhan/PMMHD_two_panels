#include "AMRStructure.hpp"

// int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, 
//                      std::vector<double>& xs_local,std::vector<double>& ys_local,
//                      std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = original;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u1s_local.data(), u2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }



// int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = original;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b1s_local.data(), b2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, 
//                      std::vector<double>& xs_local, std::vector<double>& ys_local,
//                      std::vector<double>& ws_local, double t) 
// {
//     const int n_local = (int)xs_local.size();

//     // build particle list: local domain first, then image copies
//     std::vector<double> xtemp_cpy;
//     std::vector<double> ytemp_cpy;
//     std::vector<double> ws_temp_cpy;
//     std::vector<double> u1s_temp;
//     std::vector<double> u2s_temp;

//     // x is always periodic
//     // y is periodic only when bcs == periodic_bcs
//     if (bcs == periodic_bcs) {
//         xtemp_cpy.reserve(9 * n_local);
//         ytemp_cpy.reserve(9 * n_local);
//         ws_temp_cpy.reserve(9 * n_local);
//         const double shift_x[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
//         const double shift_y[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
//         for (int img = 0; img < 9; ++img) {
//             for (int i = 0; i < n_local; ++i) {
//                 xtemp_cpy.push_back(xs_local[i] + shift_x[img]);
//                 ytemp_cpy.push_back(ys_local[i] + shift_y[img]);
//                 ws_temp_cpy.push_back(ws_local[i]);
//             }
//         }
        
//         // xtemp_cpy.reserve(25 * n_local);
//         // ytemp_cpy.reserve(25 * n_local);
//         // ws_temp_cpy.reserve(25 * n_local);

//         // // center block first
//         // for (int k = 0; k < n_local; ++k) {
//         //     xtemp_cpy.push_back(xs_local[k]);
//         //     ytemp_cpy.push_back(ys_local[k]);
//         //     ws_temp_cpy.push_back(ws_local[k]);
//         // }

//         // // then the surrounding image blocks
//         // for (int j = -2; j <= 2; ++j) {
//         //     for (int i = -2; i <= 2; ++i) {
//         //         if (i == 0 && j == 0) continue;
//         //         for (int k = 0; k < n_local; ++k) {
//         //             xtemp_cpy.push_back(xs_local[k] + i * Lx);
//         //             ytemp_cpy.push_back(ys_local[k] + j * Ly);
//         //             ws_temp_cpy.push_back(ws_local[k]);
//         //         }
//         //     }
//         // }

//         KernelMode m = periodic_y;
//         calculate_e->set_mode(m);
//         u1s_temp.assign(xtemp_cpy.size(), 0.0);
//         u2s_temp.assign(xtemp_cpy.size(), 0.0);
//         (*calculate_e)(u1s_temp.data(), u2s_temp.data(),
//                     xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());

//         // keep only the local-domain field values (first block)
//         for (int i = 0; i < n_local; ++i) {
//             u1s_local[i] = u1s_temp[i];
//             u2s_local[i] = u2s_temp[i];
//         }
//     }
//     else {
//         xtemp_cpy = xs_local;
//         ytemp_cpy = ys_local;
//         ws_temp_cpy = ws_local;
//         u1s_temp.assign(xtemp_cpy.size(), 0.0);
//         u2s_temp.assign(xtemp_cpy.size(), 0.0);
//         KernelMode m = original;
//         calculate_e->set_mode(m);
//         (*calculate_e)(u1s_temp.data(), u2s_temp.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                         ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//         for (int i = 0; i < n_local; ++i) {
//             u1s_local[i] = u1s_temp[i];
//             u2s_local[i] = u2s_temp[i];
//         }
//     }
//     return 0;
// }




// int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, 
//                      std::vector<double>& xs_local, std::vector<double>& ys_local,
//                      std::vector<double>& ws_local, double t) 
// {
//     const int n_local = (int)xs_local.size();

//     // build particle list: local domain first, then image copies
//     std::vector<double> xtemp_cpy;
//     std::vector<double> ytemp_cpy;
//     std::vector<double> ws_temp_cpy;
//     std::vector<double> b1s_temp;
//     std::vector<double> b2s_temp;

//     // x is always periodic
//     // y is periodic only when bcs == periodic_bcs
//     if (bcs == periodic_bcs) {
//         xtemp_cpy.reserve(9 * n_local);
//         ytemp_cpy.reserve(9 * n_local);
//         ws_temp_cpy.reserve(9 * n_local);
//         const double shift_x[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
//         const double shift_y[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
//         for (int img = 0; img < 9; ++img) {
//             for (int i = 0; i < n_local; ++i) {
//                 xtemp_cpy.push_back(xs_local[i] + shift_x[img]);
//                 ytemp_cpy.push_back(ys_local[i] + shift_y[img]);
//                 ws_temp_cpy.push_back(ws_local[i]);
//             }
//         }
//         KernelMode m = periodic_y;
//         calculate_e->set_mode(m);
//         b1s_temp.assign(xtemp_cpy.size(), 0.0);
//         b2s_temp.assign(xtemp_cpy.size(), 0.0);
//         (*calculate_e)(b1s_temp.data(), b2s_temp.data(),
//                     xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());

//         // keep only the local-domain field values (first block)
//         for (int i = 0; i < n_local; ++i) {
//             b1s_local[i] = b1s_temp[i];
//             b2s_local[i] = b2s_temp[i];
//         }
//     }
//     else {
//         xtemp_cpy = xs_local;
//         ytemp_cpy = ys_local;
//         ws_temp_cpy = ws_local;
//         b1s_temp.assign(xtemp_cpy.size(), 0.0);
//         b2s_temp.assign(xtemp_cpy.size(), 0.0);
//         KernelMode m = original;
//         calculate_e->set_mode(m);
//         (*calculate_e)(b1s_temp.data(), b2s_temp.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                         ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//         for (int i = 0; i < n_local; ++i) {
//             b1s_local[i] = b1s_temp[i];
//             b2s_local[i] = b2s_temp[i];
//         }
//     }
//     return 0;
// }



int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local,
                                   std::vector<double>& xs_local, std::vector<double>& ys_local,
                                   std::vector<double>& ws_local, double t)
{
    const int n_local = (int)xs_local.size();

    if (bcs == periodic_bcs) {
        // ---------- Step A: build 9x image list, 3x3 near sum at PARTICLES ----------
        const int n9 = 9 * n_local;
        std::vector<double> xtmp(n9), ytmp(n9), wtmp(n9);
        const double sx[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
        const double sy[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
        for (int img = 0; img < 9; ++img) {
            for (int i = 0; i < n_local; ++i) {
                xtmp[img * n_local + i] = xs_local[i] + sx[img];
                ytmp[img * n_local + i] = ys_local[i] + sy[img];
                wtmp[img * n_local + i] = ws_local[i];
            }
        }

        std::vector<double> u1tmp(n9, 0.0), u2tmp(n9, 0.0);
        KernelMode km = periodic_y;   // free-space 2D Biot-Savart kernel branch
        calculate_e->set_mode(km);
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xtmp.data(), n9,
                       ytmp.data(), wtmp.data(), n9);

        for (int i = 0; i < n_local; ++i) {
            u1s_local[i] = u1tmp[i];
            u2s_local[i] = u2tmp[i];
        }

        // ---------- Step B: same near sum at the 4m wall collocation targets ----------
        const std::vector<double>& wxs = periodizer->get_wall_xs();
        const std::vector<double>& wys = periodizer->get_wall_ys();
        const int mm = periodizer->get_m();
        const int nw = 4 * mm;

        std::vector<double> near_walls_u1(nw, 0.0);
        std::vector<double> near_walls_u2(nw, 0.0);

        const double pi = 3.14159265358979323846;
        const double eps = greens_epsilon;

        #pragma omp parallel for
        for (int i = 0; i < nw; ++i) {
            double tx = wxs[i];
            double ty = wys[i];
            double s1 = 0.0, s2 = 0.0;
            for (int k = 0; k < n9; ++k) {
                double dx = tx - xtmp[k];
                double dy = ty - ytmp[k];
                double r2 = dx * dx + dy * dy + eps * eps;
                s1 -= (1.0 / (2.0 * pi)) * dy / r2 * wtmp[k];
                s2 += (1.0 / (2.0 * pi)) * dx / r2 * wtmp[k];
            }
            near_walls_u1[i] = s1;
            near_walls_u2[i] = s2;
        }

        // ---------- Step C: pack 8m vector for Periodizer ----------
        std::vector<double> near_u_on_walls(8 * mm);
        for (int i = 0; i < mm; ++i) {
            near_u_on_walls[          i] = near_walls_u1[         i];  // L u1
            near_u_on_walls[   mm  + i] = near_walls_u2[         i];  // L u2
            near_u_on_walls[ 2*mm  + i] = near_walls_u1[  mm   + i];  // D u1
            near_u_on_walls[ 3*mm  + i] = near_walls_u2[  mm   + i];  // D u2
            near_u_on_walls[ 4*mm  + i] = near_walls_u1[2*mm   + i];  // R u1
            near_u_on_walls[ 5*mm  + i] = near_walls_u2[2*mm   + i];  // R u2
            near_u_on_walls[ 6*mm  + i] = near_walls_u1[3*mm   + i];  // U u1
            near_u_on_walls[ 7*mm  + i] = near_walls_u2[3*mm   + i];  // U u2
        }

        // ---------- Step D: solve Q xi = g for proxy coefficients ----------
        periodizer->compute_xi(near_u_on_walls);

        // ---------- Step E: add proxy correction at the particles ----------
        periodizer->add_correction(xs_local, ys_local, u1s_local, u2s_local);
    }
    else {
        // channel case (periodic x, open y): unchanged, use analytic channel kernel
        std::vector<double> u1tmp(n_local, 0.0), u2tmp(n_local, 0.0);
        KernelMode km = original;
        calculate_e->set_mode(km);
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xs_local.data(), n_local,
                       ys_local.data(), ws_local.data(), n_local);
        u1s_local = u1tmp;
        u2s_local = u2tmp;
    }
    return 0;
}


int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local,
                                   std::vector<double>& xs_local, std::vector<double>& ys_local,
                                   std::vector<double>& ws_local, double t)
{
    const int n_local = (int)xs_local.size();

    if (bcs == periodic_bcs) {
        // ---------- Step A: build 9x image list, 3x3 near sum at PARTICLES ----------
        const int n9 = 9 * n_local;
        std::vector<double> xtmp(n9), ytmp(n9), wtmp(n9);
        const double sx[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
        const double sy[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
        for (int img = 0; img < 9; ++img) {
            for (int i = 0; i < n_local; ++i) {
                xtmp[img * n_local + i] = xs_local[i] + sx[img];
                ytmp[img * n_local + i] = ys_local[i] + sy[img];
                wtmp[img * n_local + i] = ws_local[i];
            }
        }

        std::vector<double> b1tmp(n9, 0.0), b2tmp(n9, 0.0);
        KernelMode km = periodic_y;
        calculate_e->set_mode(km);
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xtmp.data(), n9,
                       ytmp.data(), wtmp.data(), n9);

        for (int i = 0; i < n_local; ++i) {
            b1s_local[i] = b1tmp[i];
            b2s_local[i] = b2tmp[i];
        }

        // ---------- Step B: same near sum at the 4m wall collocation targets ----------
        const std::vector<double>& wxs = periodizer->get_wall_xs();
        const std::vector<double>& wys = periodizer->get_wall_ys();
        const int mm = periodizer->get_m();
        const int nw = 4 * mm;

        std::vector<double> near_walls_b1(nw, 0.0);
        std::vector<double> near_walls_b2(nw, 0.0);

        const double pi = 3.14159265358979323846;
        const double eps = greens_epsilon;

        #pragma omp parallel for
        for (int i = 0; i < nw; ++i) {
            double tx = wxs[i];
            double ty = wys[i];
            double s1 = 0.0, s2 = 0.0;
            for (int k = 0; k < n9; ++k) {
                double dx = tx - xtmp[k];
                double dy = ty - ytmp[k];
                double r2 = dx * dx + dy * dy + eps * eps;
                s1 -= (1.0 / (2.0 * pi)) * dy / r2 * wtmp[k];
                s2 += (1.0 / (2.0 * pi)) * dx / r2 * wtmp[k];
            }
            near_walls_b1[i] = s1;
            near_walls_b2[i] = s2;
        }

        // ---------- Step C: pack 8m vector for Periodizer ----------
        std::vector<double> near_b_on_walls(8 * mm);
        for (int i = 0; i < mm; ++i) {
            near_b_on_walls[          i] = near_walls_b1[         i];  // L b1
            near_b_on_walls[   mm  + i] = near_walls_b2[         i];  // L b2
            near_b_on_walls[ 2*mm  + i] = near_walls_b1[  mm   + i];  // D b1
            near_b_on_walls[ 3*mm  + i] = near_walls_b2[  mm   + i];  // D b2
            near_b_on_walls[ 4*mm  + i] = near_walls_b1[2*mm   + i];  // R b1
            near_b_on_walls[ 5*mm  + i] = near_walls_b2[2*mm   + i];  // R b2
            near_b_on_walls[ 6*mm  + i] = near_walls_b1[3*mm   + i];  // U b1
            near_b_on_walls[ 7*mm  + i] = near_walls_b2[3*mm   + i];  // U b2
        }

        // ---------- Step D: solve Q xi = g for proxy coefficients (for B) ----------
        periodizer->compute_xi(near_b_on_walls);

        // ---------- Step E: add proxy correction at the particles ----------
        periodizer->add_correction(xs_local, ys_local, b1s_local, b2s_local);
    }
    else {
        // channel case (periodic x, open y): unchanged
        std::vector<double> b1tmp(n_local, 0.0), b2tmp(n_local, 0.0);
        KernelMode km = original;
        calculate_e->set_mode(km);
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xs_local.data(), n_local,
                       ys_local.data(), ws_local.data(), n_local);
        b1s_local = b1tmp;
        b2s_local = b2tmp;
    }
    return 0;
}






















































// int AMRStructure::evaluate_u1s_grad(std::vector<double>& u1s_grad_x_local, std::vector<double>& u1s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u1_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u1s_grad_x_local.data(), u1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_u2s_grad(std::vector<double>& u2s_grad_x_local, std::vector<double>& u2s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u2_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u2s_grad_x_local.data(), u2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_b1s_grad(std::vector<double>& b1s_grad_x_local, std::vector<double>& b1s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u1_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b1s_grad_x_local.data(), b1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_b2s_grad(std::vector<double>& b2s_grad_x_local, std::vector<double>& b2s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u2_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b2s_grad_x_local.data(), b2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_vorticity_grad(std::vector<double>& vorticity_grad_x_local, std::vector<double>& vorticity_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = vorticity_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(vorticity_grad_x_local.data(), vorticity_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_j_grad(std::vector<double>& j_grad_x_local, std::vector<double>& j_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = vorticity_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(j_grad_x_local.data(), j_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_vorticity_laplacian(std::vector<double>& vorticity_laplacian_local, std::vector<double>& vorticity_none_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = laplacian;
//     calculate_e->set_mode(m);
//     (*calculate_e)(vorticity_laplacian_local.data(), vorticity_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_j_laplacian(std::vector<double>& j_laplacian_local, std::vector<double>& j_none_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = laplacian;
//     calculate_e->set_mode(m);
//     (*calculate_e)(j_laplacian_local.data(),j_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }



int AMRStructure::compute_rhs_state(
    std::vector<double>& xs_in,
    std::vector<double>& ys_in,
    std::vector<double>& u1s_in,
    std::vector<double>& u2s_in,
    std::vector<double>& b1s_in,
    std::vector<double>& b2s_in,
    double t_in,
    std::vector<double>& u1s_grad_x,
    std::vector<double>& u1s_grad_y,
    std::vector<double>& u2s_grad_x,
    std::vector<double>& u2s_grad_y,
    std::vector<double>& b1s_grad_x,
    std::vector<double>& b1s_grad_y,
    std::vector<double>& b2s_grad_x,
    std::vector<double>& b2s_grad_y,
    std::vector<double>& source_term
) {
    cout << "enter compute_rhs_state" << endl;

    // for each leaf panel, do centered finite difference
    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
        Panel* panel = &(panels[panel_ind]);

        // only use leaf panels
        if (panel->child_inds_start > -1) {
            continue;
        }

        const int* panel_point_inds = panel->point_inds;
        double panel_xs[9], panel_ys[9];
        double panel_u1s[9], panel_u2s[9];
        double panel_b1s[9], panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            panel_xs[ii]  = xs[pind];
            panel_ys[ii]  = ys[pind];
            panel_u1s[ii] = u1s[pind];
            panel_u2s[ii] = u2s[pind];
            panel_b1s[ii] = b1s[pind];
            panel_b2s[ii] = b2s[pind];
        }

        Panel* left_panel = &(panels[panel->left_nbr_ind]);
        const int* left_panel_point_inds = left_panel->point_inds;
        double left_panel_u1s[9], left_panel_u2s[9];
        double left_panel_b1s[9], left_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = left_panel_point_inds[ii];
            left_panel_u1s[ii] = u1s[pind];
            left_panel_u2s[ii] = u2s[pind];
            left_panel_b1s[ii] = b1s[pind];
            left_panel_b2s[ii] = b2s[pind];
        }

        Panel* right_panel = &(panels[panel->right_nbr_ind]);
        const int* right_panel_point_inds = right_panel->point_inds;
        double right_panel_u1s[9], right_panel_u2s[9];
        double right_panel_b1s[9], right_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = right_panel_point_inds[ii];
            right_panel_u1s[ii] = u1s[pind];
            right_panel_u2s[ii] = u2s[pind];
            right_panel_b1s[ii] = b1s[pind];
            right_panel_b2s[ii] = b2s[pind];
        }

        double top_panel_u1s[9], top_panel_u2s[9];
        double top_panel_b1s[9], top_panel_b2s[9];
        if (panel->top_nbr_ind > -2) {
            Panel* top_panel = &(panels[panel->top_nbr_ind]);
            const int* top_panel_point_inds = top_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = top_panel_point_inds[ii];
                top_panel_u1s[ii] = u1s[pind];
                top_panel_u2s[ii] = u2s[pind];
                top_panel_b1s[ii] = b1s[pind];
                top_panel_b2s[ii] = b2s[pind];
            }
        }

        double bottom_panel_u1s[9], bottom_panel_u2s[9];
        double bottom_panel_b1s[9], bottom_panel_b2s[9];
        if (panel->bottom_nbr_ind > -2) {
            Panel* bottom_panel = &(panels[panel->bottom_nbr_ind]);
            const int* bottom_panel_point_inds = bottom_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = bottom_panel_point_inds[ii];
                bottom_panel_u1s[ii] = u1s[pind];
                bottom_panel_u2s[ii] = u2s[pind];
                bottom_panel_b1s[ii] = b1s[pind];
                bottom_panel_b2s[ii] = b2s[pind];
            }
        }

        std::vector<double> dx_u1s(9, 0.0);
        std::vector<double> dx_u2s(9, 0.0);
        std::vector<double> dx_b1s(9, 0.0);
        std::vector<double> dx_b2s(9, 0.0);
        std::vector<double> dy_u1s(9, 0.0);
        std::vector<double> dy_u2s(9, 0.0);
        std::vector<double> dy_b1s(9, 0.0);
        std::vector<double> dy_b2s(9, 0.0);

        double hx = panel_xs[3] - panel_xs[0];
        double hy = panel_ys[1] - panel_ys[0];

        //////////// u1s ///////////////
        dx_u1s[0] = (panel_u1s[3] - left_panel_u1s[3]) / (2 * hx);
        dx_u1s[3] = (panel_u1s[6] - panel_u1s[0]) / (2 * hx);
        dx_u1s[6] = (right_panel_u1s[3] - panel_u1s[3]) / (2 * hx);

        dx_u1s[1] = (panel_u1s[4] - left_panel_u1s[4]) / (2 * hx);
        dx_u1s[4] = (panel_u1s[7] - panel_u1s[1]) / (2 * hx);
        dx_u1s[7] = (right_panel_u1s[4] - panel_u1s[4]) / (2 * hx);

        dx_u1s[2] = (panel_u1s[5] - left_panel_u1s[5]) / (2 * hx);
        dx_u1s[5] = (panel_u1s[8] - panel_u1s[2]) / (2 * hx);
        dx_u1s[8] = (right_panel_u1s[5] - panel_u1s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_u1s[0] = (panel_u1s[1] - bottom_panel_u1s[1]) / (2 * hy);
            dy_u1s[3] = (panel_u1s[4] - bottom_panel_u1s[4]) / (2 * hy);
            dy_u1s[6] = (panel_u1s[7] - bottom_panel_u1s[7]) / (2 * hy);
        }
        else {
            dy_u1s[0] = (-3 * (panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2])) / (2 * hy);
            dy_u1s[3] = (-3 * (panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5])) / (2 * hy);
            dy_u1s[6] = (-3 * (panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8])) / (2 * hy);
        }

        dy_u1s[1] = (panel_u1s[2] - panel_u1s[0]) / (2 * hy);
        dy_u1s[4] = (panel_u1s[5] - panel_u1s[3]) / (2 * hy);
        dy_u1s[7] = (panel_u1s[8] - panel_u1s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_u1s[2] = (top_panel_u1s[1] - panel_u1s[1]) / (2 * hy);
            dy_u1s[5] = (top_panel_u1s[4] - panel_u1s[4]) / (2 * hy);
            dy_u1s[8] = (top_panel_u1s[7] - panel_u1s[7]) / (2 * hy);
        }
        else {
            dy_u1s[2] = (-3 * (panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2])) / (2 * hy);
            dy_u1s[5] = (-3 * (panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5])) / (2 * hy);
            dy_u1s[8] = (-3 * (panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8])) / (2 * hy);
        }

        //////////// u2s ///////////////
        dx_u2s[0] = (panel_u2s[3] - left_panel_u2s[3]) / (2 * hx);
        dx_u2s[3] = (panel_u2s[6] - panel_u2s[0]) / (2 * hx);
        dx_u2s[6] = (right_panel_u2s[3] - panel_u2s[3]) / (2 * hx);

        dx_u2s[1] = (panel_u2s[4] - left_panel_u2s[4]) / (2 * hx);
        dx_u2s[4] = (panel_u2s[7] - panel_u2s[1]) / (2 * hx);
        dx_u2s[7] = (right_panel_u2s[4] - panel_u2s[4]) / (2 * hx);

        dx_u2s[2] = (panel_u2s[5] - left_panel_u2s[5]) / (2 * hx);
        dx_u2s[5] = (panel_u2s[8] - panel_u2s[2]) / (2 * hx);
        dx_u2s[8] = (right_panel_u2s[5] - panel_u2s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_u2s[0] = (panel_u2s[1] - bottom_panel_u2s[1]) / (2 * hy);
            dy_u2s[3] = (panel_u2s[4] - bottom_panel_u2s[4]) / (2 * hy);
            dy_u2s[6] = (panel_u2s[7] - bottom_panel_u2s[7]) / (2 * hy);
        }
        else {
            dy_u2s[0] = (-3 * (panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2])) / (2 * hy);
            dy_u2s[3] = (-3 * (panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5])) / (2 * hy);
            dy_u2s[6] = (-3 * (panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8])) / (2 * hy);
        }

        dy_u2s[1] = (panel_u2s[2] - panel_u2s[0]) / (2 * hy);
        dy_u2s[4] = (panel_u2s[5] - panel_u2s[3]) / (2 * hy);
        dy_u2s[7] = (panel_u2s[8] - panel_u2s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_u2s[2] = (top_panel_u2s[1] - panel_u2s[1]) / (2 * hy);
            dy_u2s[5] = (top_panel_u2s[4] - panel_u2s[4]) / (2 * hy);
            dy_u2s[8] = (top_panel_u2s[7] - panel_u2s[7]) / (2 * hy);
        }
        else {
            dy_u2s[2] = (-3 * (panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2])) / (2 * hy);
            dy_u2s[5] = (-3 * (panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5])) / (2 * hy);
            dy_u2s[8] = (-3 * (panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8])) / (2 * hy);
        }

        //////////// b1s ///////////////
        dx_b1s[0] = (panel_b1s[3] - left_panel_b1s[3]) / (2 * hx);
        dx_b1s[3] = (panel_b1s[6] - panel_b1s[0]) / (2 * hx);
        dx_b1s[6] = (right_panel_b1s[3] - panel_b1s[3]) / (2 * hx);

        dx_b1s[1] = (panel_b1s[4] - left_panel_b1s[4]) / (2 * hx);
        dx_b1s[4] = (panel_b1s[7] - panel_b1s[1]) / (2 * hx);
        dx_b1s[7] = (right_panel_b1s[4] - panel_b1s[4]) / (2 * hx);

        dx_b1s[2] = (panel_b1s[5] - left_panel_b1s[5]) / (2 * hx);
        dx_b1s[5] = (panel_b1s[8] - panel_b1s[2]) / (2 * hx);
        dx_b1s[8] = (right_panel_b1s[5] - panel_b1s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_b1s[0] = (panel_b1s[1] - bottom_panel_b1s[1]) / (2 * hy);
            dy_b1s[3] = (panel_b1s[4] - bottom_panel_b1s[4]) / (2 * hy);
            dy_b1s[6] = (panel_b1s[7] - bottom_panel_b1s[7]) / (2 * hy);
        }
        else {
            dy_b1s[0] = (-3 * (panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2])) / (2 * hy);
            dy_b1s[3] = (-3 * (panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5])) / (2 * hy);
            dy_b1s[6] = (-3 * (panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8])) / (2 * hy);
        }

        dy_b1s[1] = (panel_b1s[2] - panel_b1s[0]) / (2 * hy);
        dy_b1s[4] = (panel_b1s[5] - panel_b1s[3]) / (2 * hy);
        dy_b1s[7] = (panel_b1s[8] - panel_b1s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_b1s[2] = (top_panel_b1s[1] - panel_b1s[1]) / (2 * hy);
            dy_b1s[5] = (top_panel_b1s[4] - panel_b1s[4]) / (2 * hy);
            dy_b1s[8] = (top_panel_b1s[7] - panel_b1s[7]) / (2 * hy);
        }
        else {
            dy_b1s[2] = (-3 * (panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2])) / (2 * hy);
            dy_b1s[5] = (-3 * (panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5])) / (2 * hy);
            dy_b1s[8] = (-3 * (panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8])) / (2 * hy);
        }

        //////////// b2s ///////////////
        dx_b2s[0] = (panel_b2s[3] - left_panel_b2s[3]) / (2 * hx);
        dx_b2s[3] = (panel_b2s[6] - panel_b2s[0]) / (2 * hx);
        dx_b2s[6] = (right_panel_b2s[3] - panel_b2s[3]) / (2 * hx);

        dx_b2s[1] = (panel_b2s[4] - left_panel_b2s[4]) / (2 * hx);
        dx_b2s[4] = (panel_b2s[7] - panel_b2s[1]) / (2 * hx);
        dx_b2s[7] = (right_panel_b2s[4] - panel_b2s[4]) / (2 * hx);

        dx_b2s[2] = (panel_b2s[5] - left_panel_b2s[5]) / (2 * hx);
        dx_b2s[5] = (panel_b2s[8] - panel_b2s[2]) / (2 * hx);
        dx_b2s[8] = (right_panel_b2s[5] - panel_b2s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_b2s[0] = (panel_b2s[1] - bottom_panel_b2s[1]) / (2 * hy);
            dy_b2s[3] = (panel_b2s[4] - bottom_panel_b2s[4]) / (2 * hy);
            dy_b2s[6] = (panel_b2s[7] - bottom_panel_b2s[7]) / (2 * hy);
        }
        else {
            dy_b2s[0] = (-3 * (panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2])) / (2 * hy);
            dy_b2s[3] = (-3 * (panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5])) / (2 * hy);
            dy_b2s[6] = (-3 * (panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8])) / (2 * hy);
        }

        dy_b2s[1] = (panel_b2s[2] - panel_b2s[0]) / (2 * hy);
        dy_b2s[4] = (panel_b2s[5] - panel_b2s[3]) / (2 * hy);
        dy_b2s[7] = (panel_b2s[8] - panel_b2s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_b2s[2] = (top_panel_b2s[1] - panel_b2s[1]) / (2 * hy);
            dy_b2s[5] = (top_panel_b2s[4] - panel_b2s[4]) / (2 * hy);
            dy_b2s[8] = (top_panel_b2s[7] - panel_b2s[7]) / (2 * hy);
        }
        else {
            dy_b2s[2] = (-3 * (panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2])) / (2 * hy);
            dy_b2s[5] = (-3 * (panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5])) / (2 * hy);
            dy_b2s[8] = (-3 * (panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8])) / (2 * hy);
        }

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            u1s_grad_x[pind]       = dx_u1s[ii];
            u2s_grad_x[pind]       = dx_u2s[ii];
            b1s_grad_x[pind]       = dx_b1s[ii];
            b2s_grad_x[pind]       = dx_b2s[ii];

            u1s_grad_y[pind]       = dy_u1s[ii];
            u2s_grad_y[pind]       = dy_u2s[ii];
            b1s_grad_y[pind]       = dy_b1s[ii];
            b2s_grad_y[pind]       = dy_b2s[ii];
        }
    }

    for (int i = 0; i < xs.size(); ++i) {
        source_term[i] = 2*(b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i])
                        -2* (b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i]);
    }
    return 0;
}

