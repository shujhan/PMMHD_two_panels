#include "Panel.hpp"
// Constructors
Panel::Panel() 
    : panel_ind(0)
{
    for(int i = 0; i < 9; i++) {
        point_inds[i] = i;
    }
}

Panel::Panel(int panel_ind)
    
{   this->panel_ind = panel_ind;
    // for(int i = 0; i < 9; i++) {
    //     this->point_inds[i] = point_inds[i];
    // }
}

void Panel::set_point_inds(int p0, int p1, int p2, int p3, int p4, 
                    int p5, int p6, int p7, int p8) {
    point_inds[0] = p0;
    point_inds[1] = p1;
    point_inds[2] = p2;
    point_inds[3] = p3;
    point_inds[4] = p4;
    point_inds[5] = p5;
    point_inds[6] = p6;
    point_inds[7] = p7;
    point_inds[8] = p8;
}