#include "AMRStructure.hpp"


// void AMRStructure::interpolate_to_initial_xys(
//     std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& ps, 
//     int nx, int np, bool verbose) 
// {
//     std::vector<double> shifted_xs(xs.size());
//     if (bcs == 0) {
//         shift_xs(shifted_xs, xs, ps);
//         if (verbose) {
//             cout << "shift xs" << endl;
//             std::copy(shifted_xs.begin(), shifted_xs.end(), std::ostream_iterator<double>(cout, ", "));
//             cout << endl;
//         }
//     }

// #ifdef DEBUG
// cout << "sorting " << endl;
// #endif
//     std::vector<int> sort_indices(xs.size());
//     for (int ii = 0; ii < xs.size(); ii++ ) { sort_indices[ii] = ii; }
//     // std::iota(sort_indices.begin(), sort_indices.end(), 0);
//     bool do_sort = true;
//     double sort_threshold = initial_dp / 10.0;
//     if (do_sort) {
//         std::sort(sort_indices.begin(), sort_indices.end(),
//             [&] (int a, int b) 
//             { 
//                 if (fabs(ps[a] - ps[b]) >= sort_threshold) { return ps[a] < ps[b]; }
//                 else {
//                     return shifted_xs[a] < shifted_xs[b];
//                 }
//             });
//         // std::cout << "sorted indices " << endl;
//         // std::copy(sort_indices.begin(), sort_indices.end(), std::ostream_iterator<int>(cout, " "));
//         // std::cout << endl;
//         // cout << "xs.size " << xs.size() << endl;
//         // cout << "sort_indices.size " << sort_indices.size() << endl;
//         // cout << "shifted_xs.size " << shifted_xs.size() << endl;
//     }
// #ifdef DEBUG
// cout << "Done sorting" << endl;
// #endif
//     std::vector<double> sortxs(shifted_xs.size()), sortps(ps.size());
//     for (int ii = 0; ii < xs.size(); ii++) {
//         sortxs[ii] = shifted_xs[sort_indices[ii]];
//         sortps[ii] = ps[sort_indices[ii]];
//     }
//     // if (verbose) {
//     #ifdef DEBUG_L3
//         cout << "sorted xs, size = " << sortxs.size() << endl;
//         std::copy(sortxs.begin(), sortxs.end(), std::ostream_iterator<double>(cout, ", "));
//         cout << endl << "sorted ps, size = "  << sortps.size() << endl;
//         std::copy(sortps.begin(), sortps.end(), std::ostream_iterator<double>(cout, ", "));
//         cout << endl;
//     #endif /* DEBUG_L3 */
//     // }

//     stop = high_resolution_clock::now();
//     // duration = duration_cast<microseconds>(stop - start);
//     // if (print_profile) {
//     //     cout << "sort time " << duration.count() << " microseconds" << endl;
//     // }
    
//     auto search_start = high_resolution_clock::now();
//     std::vector<double> sortfs(xs.size() );

//     std::vector<int> leaf_panel_of_points(xs.size() );
//     std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );

// #ifdef DEBUG
// // counting how many xs to ps
// // int ind_debug = 0;
// // std::vector<int> size_p;
// // std::vector<double> unique_ps;
// // size_p.reserve(np);
// // unique_ps.reserve(np);
// // while (ind_debug < sortxs.size() ) {
// //     int cntr = 0;
// //     double debug_p = sortps[ind_debug];
// //     unique_ps.push_back(debug_p);
// //     cntr++;
// //     int ind_debug2 = ind_debug + 1;
// //     while (fabs(sortps[ind_debug2] - debug_p) < sort_threshold) {
// //         ind_debug2++;
// //         cntr++;
// //     }
// //     ind_debug += cntr;
// //     size_p.push_back(cntr);
// // }
// // cout << "np = " << np << ", and there are " << size_p.size() << " distinct p values" << endl;
// // for (int ii = 0; ii < size_p.size(); ++ii) {
// //     cout << "at v= " << unique_ps[ii] << " there are " << size_p[ii] << " points" << endl;
// // }
// // verbose = false;
// cout << "finding panel of first point" << endl;
// // verbose=true;
// #endif
//     bool beyond_boundary = false;
//     int leaf_ind = find_leaf_containing_xp_recursively(sortxs[0], sortps[0], beyond_boundary, 0, verbose);
// #ifdef DEBUG
// cout << "found first panel" << endl;
// #endif
//     std::vector<int> first_column_leaf_inds(np);
//     // std::vector<int> first_column_ind_in_old_mesh(np);

// #ifdef DEBUG 
// cout << "searching first column" << endl;
// // if (iter_num >= 4) { verbose = true;} 
//     // if (verbose) {
//         std::cout << "nx x np= " << nx << " x " << np << endl;
//         cout << "xs size: " << xs.size() << endl;
//     // }
// // cout << "First column points" << endl;
// // for (int ii = 0; ii < np; ++ii) {
// //     int point_ind = ii * np;
// //     cout << "sort point " << point_ind << " (x,p)=(" << sortxs[point_ind] << ", " << sortps[point_ind] << ")" << endl;
// // }
// #endif

//     // if (bcs == periodic_bcs) {
//         for (int ii =0; ii < np; ++ii) {
//             beyond_boundary = false;
//             int point_ind = ii * nx;
//             std::set<int> history;
//             history.emplace(leaf_ind);
//             #ifdef DEBUG_L2
//             cout << "testing point " << point_ind << ", x= " << sortxs[point_ind] << ", p= " << sortps[point_ind] << endl;
//             #endif
//             leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], beyond_boundary, 
//                                                                 leaf_ind, history, verbose);
//             first_column_leaf_inds[ii] = leaf_ind;
//             // point_in_leaf_panels_by_inds[leaf_ind].push_back(point_ind);
//             if (beyond_boundary) {
//                 leaf_panel_of_points[point_ind] = 0;
//             } else {
//                 leaf_panel_of_points[point_ind] = leaf_ind;
//             }
//             // cout << "point " << sort_indices[point_ind] << " in panel " << leaf_ind << " (sorted ind " << point_ind << ")" << endl;
//         }
// //     } else { // open bcs
// //         for (int ii =0; ii < np; ++ii) {
// //             int jj = 0;
// //             int point_ind = ii * np;
// //             std::set<int> history;
// //             history.emplace(leaf_ind);
// //             int temp_leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], leaf_ind, history, verbose);
// //             leaf_panel_of_points[point_ind] = temp_leaf_ind;

// //             while (temp_leaf_ind == 0 && jj < nx) {
// // #ifdef DEBUG
// // cout << "temp_leaf_ind was assigned 0!" << endl;
// // #endif
// //                 jj++;
// //                 point_ind++;
// //                 temp_leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], leaf_ind, history, verbose);
// //                 leaf_panel_of_points[point_ind] = temp_leaf_ind;
// //             }
// //             leaf_ind = temp_leaf_ind;
// //             first_column_ind_in_old_mesh[ii] = jj;
// //         }    
// //     }
// #ifdef DEBUG
// verbose = false;
// cout << "after first column" << endl;
// // cout << "first column ind in old mesh" << endl;
// // std::copy(first_column_ind_in_old_mesh.begin(), first_column_ind_in_old_mesh.end(),
// //         std::ostream_iterator<int>(cout, " "));
// // cout << endl;

// //     cout << "leaf_panel_of_points " << endl;
// //     for (int ii = 0; ii < xs.size(); ++ii) {
// //         cout << "point (sorted ind) " << ii <<", unsorted ind " << sort_indices[ii] << ": (x,p)=(" << sortxs[ii] << ", " << sortps[ii] << ") is in panel " << leaf_panel_of_points[ii] << endl;
// //     }
//     // if (verbose) {
//         cout << "panel search" << endl;
//     // }
// #endif
// // #pragma omp parallel
// // {
//     // int num_threads = omp_get_num_threads();
//     // if (omp_get_thread_num() == 0) {
//     //     cout << "Number of threads from calc_E: " << num_threads << endl;
//     // }
//     // #pragma omp for
//     for (int ii = 0; ii < np; ++ii) {
//         // int jj0 = first_column_ind_in_old_mesh[ii];
//         int jj0 = 0;
//         int point_ind = ii * nx + jj0;
//         int leaf_ind_c = first_column_leaf_inds[ii];//leaf_panel_of_points[point_ind];
//         for (int jj = jj0+1; jj < nx; ++jj) {
//             beyond_boundary = false;
//             point_ind++;
//             #ifdef DEBUG_L2
//             // if (verbose) {
//             cout << "Testing point " << point_ind << ", (x,p)= (" << sortxs[point_ind] << ", " << sortps[point_ind] << ")" << endl;
//             // }
//             #endif
//             std::set<int> history;
//             history.emplace(leaf_ind_c);
//             leaf_ind_c = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], beyond_boundary, leaf_ind_c, history, verbose);
//             // point_in_leaf_panels_by_inds[leaf_ind_c].push_back(point_ind);
//             if (beyond_boundary) {
//                 leaf_panel_of_points[point_ind] = 0;
//             } else {
//                 leaf_panel_of_points[point_ind] = leaf_ind_c;
//             }
//             #ifdef DEBUG_L3
//             cout << "sorted ind " << point_ind << ", ind " << sort_indices[point_ind] << " is in panel " << leaf_ind_c << endl;
//             #endif
//         }
//     }    
// // } // end omp parallel
// #ifdef DEBUG
// cout << "done with panel search" << endl;
// #endif
// #ifdef DEBUG_L3
//     cout << "leaf_panel_of_points " << endl;
//     for (int ii = 0; ii < xs.size(); ++ii) {
//         cout << "point (sorted ind) " << ii <<", unsorted ind " << sort_indices[ii] << ": (x,p)=(" << sortxs[ii] << ", " << sortps[ii] << ") is in panel " << leaf_panel_of_points[ii] << endl;
//     }
// #endif

//     for (int ii = 0; ii < leaf_panel_of_points.size(); ++ii) {
//         point_in_leaf_panels_by_inds[leaf_panel_of_points[ii]].push_back(ii);
//     }

// #ifdef DEBUG_L3

//     cout << "point in leaf panels by inds " << endl;
//     for (auto ii = 0; ii < point_in_leaf_panels_by_inds.size(); ++ii) {
//         std::vector<int> p_leaf = point_in_leaf_panels_by_inds[ii];
//         if (p_leaf.size() > 0 ) {
//             cout << "points in panel " << ii << ": ";
//             for (auto jj = 0; jj < p_leaf.size(); ++jj) {
//                 cout << p_leaf[jj] << ", ";
//             }
//             cout << endl;
//         }
//     }
//     // std::copy(leaf_panel_of_points.begin(), leaf_panel_of_points.end(), std::ostream_iterator<int>(cout, " "));
// // debug
//     // cout << "sort_indices[498] = " << sort_indices[498] << endl;
//     // cout << "leaf_panel_of_points[498] = " << leaf_panel_of_points[498] << endl;
//     // cout << "leaf_panel_of_points[sort_indicies[498]] = " << leaf_panel_of_points[sort_indices[498]] << endl;
// //
// #endif
//     auto search_stop = high_resolution_clock::now();
//     add_time(search_time, duration_cast<duration<double>>(search_stop - search_start) );
    
// #ifdef DEBUG
// cout << "eval interpolant" << endl;
// #endif
//     start = high_resolution_clock::now();
//     for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
//         if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
//             interpolate_from_panel_to_points(sortfs,sortxs,sortps,point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
//         }
//     }
//     for (int ii = 0; ii < fs.size(); ii++) {
//         // debugging
//         // if (sort_indices[ii] == 498) {
//             // cout << "(x,p)_498 = (" << xs[498] << ", " << ps[498] << ")" << endl;
//             // cout << "fs(498) = " << fs[498] << endl;
//         // }
//         // end debugging
//         fs[sort_indices[ii]] = sortfs[ii];
//     }
//     // debugging
//     // cout << "done unpacking sort" << endl;
//     // cout << "something's up with sort_indices: " << endl;
//     // std::copy(sort_indices.begin(), sort_indices.end(), std::ostream_iterator<int>(cout, " "));
//     // cout << endl;
//     // cout << "(x,p)_498 = (" << xs[498] << ", " << ps[498] << ")" << endl;
//     // cout << "fs(498) = " << fs[498] << endl;
//     // end debug element
//     stop = high_resolution_clock::now();
//     add_time(eval_time, duration_cast<duration<double>>(stop - start) );

//     #ifdef DEBUG
//     cout << "Done evaluating interpolant and done interpolating onto grid" << endl;
//     #endif
// }