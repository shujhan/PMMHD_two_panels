#include "AMRStructure.hpp"


void AMRStructure::generate_mesh(std::function<double (double,double)> f, 
                                 bool do_adaptive_refine, bool is_initial_step) 
{
    // bool verbose=false;

    // auto start = high_resolution_clock::now();
    create_prerefined_mesh();
    // auto stop = high_resolution_clock::now();
    // add_time(tree_build_time,  duration_cast<duration<double>>(stop - start) );
}


int AMRStructure::create_prerefined_mesh() {
    // create a level 1 panel, 2^2 + 1 = 5 points in x and y
    // in total 5 * 5 = 25 points 

    if (initial_height < 1) {
        throw std::invalid_argument("height must be greater than 1");
    }

    double dx = Lx / 4;
    double dy = Ly / 4;
    std::vector<double> xs_init, ys_init; 
    // get the x, y value for level 0 panel
    for(int ii = 0; ii < 5; ++ii) {
        xs_init.push_back(x_min + ii * dx);
        ys_init.push_back(y_min + ii * dy);
    }

    panels.clear();
    xs.clear();
    ys.clear();
    fs.clear();
    xs.reserve(25);
    ys.reserve(25);
    // put 9 points to xs and ys: level 0 panel 
    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 0; jj < 5; jj+=2) { 
            xs.push_back(xs_init[ii]); 
            ys.push_back(ys_init[jj]);
        }
    }
    // put points for level 1 panel 
    for (int ii = 0; ii < 2; ++ii) {
        xs.push_back(xs_init[2*ii]); xs.push_back(xs_init[2*ii]);
        ys.push_back(ys_init[1]); ys.push_back(ys_init[3]);
        for (int jj = 0; jj < 5; ++jj) {
            xs.push_back(xs_init[1 + 2*ii]);
            ys.push_back(ys_init[jj]);
        }
    }
    // finish level 1 panel 
    for (int jj = 1; jj < 5; jj+=2) {
        xs.push_back(xs_init[4]);
        ys.push_back(ys_init[jj]);
    }

    fs = std::vector<double>(25, 1.0);

    if (bcs == 0) { // periodic_bcs
        panels.push_back(Panel{});
        panels[0].is_left_bdry = true;
        panels[0].is_right_bdry = true;
        panels.push_back(Panel{1, 1, 0, 0, 3, 2, 3, -2});
        panels[1].is_left_bdry = true;
        panels.push_back(Panel{2, 1, 0, 1, 4, -2, 4, 1});
        panels[2].is_left_bdry = true;
        panels.push_back(Panel{3, 1, 0, 2, 1, 4, 1, -2});
        panels[3].is_right_bdry = true;
        panels.push_back(Panel{4, 1, 0, 3, 2,-2,2,3});
        panels[4].is_right_bdry = true;
    } 

    panels[1].set_point_inds(0,9,1,11,12,13,3,16,4);
    panels[1].needs_refinement = true;
    panels[2].set_point_inds(1,10,2,13,14,15,4,17,5);
    panels[2].needs_refinement = true;
    panels[3].set_point_inds(3,16,4,18,19,20,6,23,7);
    panels[3].needs_refinement = true;
    panels[4].set_point_inds(4,17,5,20,21,22,7,24,8);
    panels[4].needs_refinement = true;

    panels[0].set_child_inds_start(1);
    minimum_unrefined_index = 1;

    // finish level 1, 25 points and 5 panels in total 


    // call refine
    // first do uniform mesh 
    for (int level = 1; level < initial_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
        }
        refine_panels( [] (double x, double y) {return 1.0;} , false);
        minimum_unrefined_index = num_panels_pre_refine;
    }

    is_initial_mesh_set = true;
    return 0;

}



void AMRStructure::refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine) {
    std::vector <double> new_xs;
    std::vector <double> new_ys;
    std::vector <double> new_fs;
    std::vector <int> prospective_leaf_inds;
    // int new_vert_ind = particles.size();
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();
    

    // for each panel in this height 
    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        Panel* panel= &(panels[jj]); 
        if (panel->needs_refinement ) {
            // printf("refining panel %i\n", panel->get_panel_ind() );
            std::vector<double> panel_xs;
            std::vector<double> panel_ys;
            double dx, dy;

            const int* panel_points = panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int point_ind = panel_points[ii];
                panel_xs.push_back(xs[point_ind]);
                panel_ys.push_back(ys[point_ind]);
            }
            dx = panel_xs[3] - panel_xs[0];
            dy = panel_ys[1] - panel_ys[0];
            double sub_dx = 0.5 * dx;
            double sub_dy = 0.5 * dy;

            int num_new_panels = panels.size();
            // double x_left = panel_xs[0];
            // double x_mid = x_left + .5 * dx;
            // double x_right = panel_xs[2];
            // double p_bottom = panel_ps[0];
            // double p_mid = p_bottom + .5 * dp;
            // double p_top = panel_ps[1];
            double subpanel_xs[5], subpanel_ys[5];
            for (int ii = 0; ii < 5; ii ++) {
                subpanel_xs[ii] = panel_xs[0] + sub_dx * ii;
                subpanel_ys[ii] = panel_ys[0] + sub_dy * ii;
            }

            // int left_vert_ind, bottom_vert_ind, mid_vert_ind, top_vert_ind, right_vert_ind;
            int point_9_ind, point_10_ind, point_11_ind, point_15_ind, 
                point_18_ind, point_22_ind, point_23_ind, point_24_ind;
            int child_0_bottom_nbr_ind = -1;
            int child_0_left_nbr_ind = -1;
            int child_1_left_nbr_ind = -1;
            int child_1_top_nbr_ind = -1;
            int child_2_bottom_nbr_ind = -1;
            int child_2_right_nbr_ind = -1;
            int child_3_right_nbr_ind = -1;
            int child_3_top_nbr_ind = -1;

            // generate new vertices

            Panel* panel_parent;
            // check left neighbor
            if (panel->left_nbr_ind == -2) { // no left neighbor 
                child_0_left_nbr_ind = -2;
                child_1_left_nbr_ind = -2;
                point_9_ind = new_vert_ind++;
                point_10_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } 
            else if (panel->left_nbr_ind == -1) { // Not refined in parent panel, need go back to parent 
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_left = &(panels[panel_parent->left_nbr_ind]);
                if (!parent_left->is_refined_xy) {
                    parent_left->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_9_ind = new_vert_ind;
                // point_10_ind = new_vert_ind++;
                point_10_ind = point_9_ind + 1;
                new_vert_ind += 2;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } 
            else {
                Panel* panel_left = &(panels[panel->left_nbr_ind]);
                if (! panel_left->is_refined_xy) {
                    point_9_ind = new_vert_ind++;
                    point_10_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                }
                else {
                    child_0_left_nbr_ind = panel_left->child_inds_start +2;
                    Panel* child_0_left_nbr = &(panels[child_0_left_nbr_ind]);
                    child_0_left_nbr->right_nbr_ind = num_new_panels;
                    child_1_left_nbr_ind = panel_left->child_inds_start + 3;
                    Panel* child_1_left_nbr = &(panels[child_1_left_nbr_ind]);
                    child_1_left_nbr->right_nbr_ind = num_new_panels + 1;
                    // if (panel->is_left_bdry && bcs==periodic_bcs) {
                    if (panel->is_left_bdry && bcs==0) {
                        point_9_ind = new_vert_ind++;
                        point_10_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_9_ind = child_0_left_nbr->point_inds[7];
                        point_10_ind = child_1_left_nbr->point_inds[7];
                    }
                }
            }
            // check bottom neighbor
            if (panel->bottom_nbr_ind == -2) {
                child_0_bottom_nbr_ind = -2;
                child_2_bottom_nbr_ind = -2;
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
            } else if (panel->bottom_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_bottom = &(panels[panel_parent->bottom_nbr_ind]);
                if (!parent_bottom->is_refined_xy ) {
                    parent_bottom->needs_refinement = true;
                    need_further_refinement = true;
                    #ifdef DEBUG
                    cout << "refine: setting refinement flag in panel " << jj << endl;
                    #endif
                }
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
            } else {
                Panel* panel_bottom = &(panels[panel->bottom_nbr_ind]);
                if (! panel_bottom->is_refined_xy ) {
                    point_11_ind = new_vert_ind++;
                    point_18_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
                }
                else {
                    child_0_bottom_nbr_ind = panel_bottom->child_inds_start + 1;
                    Panel* child_0_bottom_nbr = &(panels[child_0_bottom_nbr_ind]);
                    child_0_bottom_nbr->top_nbr_ind = num_new_panels;
                    child_2_bottom_nbr_ind = panel_bottom->child_inds_start +3;
                    Panel* child_2_bottom_nbr = &(panels[child_2_bottom_nbr_ind]);
                    child_2_bottom_nbr->top_nbr_ind = num_new_panels + 2;

                    point_11_ind = child_0_bottom_nbr->point_inds[5];
                    point_18_ind = child_2_bottom_nbr->point_inds[5];
                }
            }

            // check top neighbor
            if (panel->top_nbr_ind == -2) {
                child_1_top_nbr_ind = -2;
                child_3_top_nbr_ind = -2;
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
            } else if (panel->top_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_top = &(panels[panel_parent->top_nbr_ind]);
                if (!parent_top->is_refined_xy ) {
                    parent_top->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
            } else {
                Panel* panel_top = &(panels[panel->top_nbr_ind]);
                if (! panel_top->is_refined_xy ) {
                    point_15_ind = new_vert_ind++;
                    point_22_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
                }
                else {
                    child_1_top_nbr_ind = panel_top->child_inds_start;
                    Panel* child_1_top_nbr = &(panels[child_1_top_nbr_ind]);
                    child_1_top_nbr->bottom_nbr_ind = num_new_panels + 1;
                    child_3_top_nbr_ind = panel_top->child_inds_start + 2;
                    Panel* child_3_top_nbr = &(panels[child_3_top_nbr_ind]);
                    child_3_top_nbr->bottom_nbr_ind = num_new_panels + 3;

                    point_15_ind = child_1_top_nbr->point_inds[3];
                    point_22_ind = child_3_top_nbr->point_inds[3];
                }
            }

            // check right neighbor
            if (panel->right_nbr_ind == -2) {
                child_2_right_nbr_ind = -2;
                child_3_right_nbr_ind = -2;
                point_23_ind = new_vert_ind++;
                point_24_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else if (panel->right_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_right = &(panels[panel_parent->right_nbr_ind]);
                if (!parent_right->is_refined_xy ) {
                    parent_right->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_23_ind = new_vert_ind++;
                point_24_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else {
                Panel* panel_right = &(panels[panel->right_nbr_ind]);
                if (! panel_right->is_refined_xy) {
                    point_23_ind = new_vert_ind++;
                    point_24_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    if (panel->right_nbr_ind == jj /* panel_ind = jj */) {
                        child_2_right_nbr_ind = num_new_panels;
                        child_0_left_nbr_ind = num_new_panels + 2;
                        child_3_right_nbr_ind = num_new_panels + 1;
                        child_1_left_nbr_ind = num_new_panels + 3;
                    }
                }
                else {
                    child_2_right_nbr_ind = panel_right->child_inds_start;
                    Panel* child_2_right_nbr = &(panels[child_2_right_nbr_ind]);
                    child_2_right_nbr->left_nbr_ind = num_new_panels+2;
                    child_3_right_nbr_ind = panel_right->child_inds_start + 1;
                    Panel* child_3_right_nbr = &(panels[child_3_right_nbr_ind]);
                    child_3_right_nbr->left_nbr_ind = num_new_panels + 3;

                    // if (panel->is_right_bdry and bcs==periodic_bcs) {
                    if (panel->is_right_bdry && bcs==0) {
                        point_23_ind = new_vert_ind++;
                        point_24_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_23_ind = child_2_right_nbr->point_inds[1];
                        point_24_ind = child_3_right_nbr->point_inds[1];
                    }
                }
            } // end check right neighbor

            // add interior points
            int point_12_ind = new_vert_ind;
            int point_16_ind = point_12_ind + 3;
            int point_19_ind = point_16_ind + 2;
            for (int ii = 0; ii < 3; ii++) { // point 12, 13, 14
                new_xs.push_back(subpanel_xs[1]);
                new_ys.push_back(subpanel_ys[1+ii]);
            }
            for (int ii = 0; ii < 2; ii++) { // point 16, 17
                new_xs.push_back(subpanel_xs[2]);
                new_ys.push_back(subpanel_ys[1+2*ii]);
            }
            for (int ii = 0; ii < 3; ii++) {// point 19, 20, 21
                new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[1+ii]);
            }
            new_vert_ind += 8;

            // generate new panels
            // add these to list of prospective_panel_indices
            if (do_adaptive_refine) {
                for (int ii = num_new_panels; ii < num_new_panels + 4; ++ii) {
                    prospective_leaf_inds.push_back(ii);
                }
            }
            // panel->child_inds_start = num_new_panels;
            panel->set_child_inds_start(num_new_panels);
            // printf("post refinement, panel looks like:\n");
            // panel->print_panel();
            int child_level = panel->level + 1;
            int panel_ind = panel->panel_ind;
            int* point_inds = panel->point_inds;
            // for (int ii = 0; ii < ; ii++) {
            //     panel_vertex_inds[ii] = panel->point_inds[ii];
            // }
            panels.push_back(Panel {num_new_panels, child_level, panel_ind, 0, 
                    point_inds[0], point_9_ind, point_inds[1],
                    point_11_ind, point_12_ind, point_12_ind + 1,
                    point_inds[3], point_16_ind, point_inds[4],
                    child_0_left_nbr_ind, num_new_panels + 1, 
                    num_new_panels + 2, child_0_bottom_nbr_ind,
                    panel->is_left_bdry, false});
            panels.push_back(Panel {num_new_panels+1, child_level, panel_ind, 1, 
                    point_inds[1], point_10_ind, point_inds[2],
                    point_12_ind+1, point_12_ind+2, point_15_ind,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    child_1_left_nbr_ind, child_1_top_nbr_ind,
                    num_new_panels + 3, num_new_panels,
                    panel->is_left_bdry, false});
            panels.push_back(Panel {num_new_panels+2, child_level, panel_ind, 2, 
                    point_inds[3], point_16_ind, point_inds[4],
                    point_18_ind, point_19_ind, point_19_ind+1,
                    point_inds[6], point_23_ind, point_inds[7],
                    num_new_panels, num_new_panels+3, 
                    child_2_right_nbr_ind, child_2_bottom_nbr_ind,
                    false, panel->is_right_bdry});
            panels.push_back(Panel {num_new_panels+3, child_level, panel_ind, 3,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    point_19_ind+1, point_19_ind+2, point_22_ind,
                    point_inds[7], point_24_ind, point_inds[8],
                    num_new_panels+1, child_3_top_nbr_ind,
                    child_3_right_nbr_ind, num_new_panels+2,
                    false, panel->is_right_bdry});

        } // end if panel is flagged
    } //end for loop through panels

    // set fs
    new_fs.reserve(new_xs.size() );
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_fs.push_back( f(new_xs.at(ii), new_ys.at(ii)) );
    }

    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); ys.push_back(new_ys[ii]); fs.push_back(new_fs[ii]);
        // particles.push_back(Particle(new_xs.at(ii), new_ps.at(ii), new_fs.at(ii), 0.0));
    }


    // test 
    // if (do_adaptive_refine) { 
    //     for (int ii = 0; ii < prospective_leaf_inds.size(); ++ii) {
    //         test_panel(prospective_leaf_inds.at(ii), false);
    //     }
    // }
}




