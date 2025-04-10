#include <iostream>

#include "tinympc/admm.hpp"

#define DEBUG_MODULE "TINYALG"

extern "C"
{


double a_max = 17;
double theta = 25 * M_PI / 180;
double scix = 1;
double sciy = 1;
double sciz = 0.1;
/**
 * Update linear terms from Riccati backward pass
*/
void backward_pass_grad(TinySolver *solver)
{
    for (int i = solver->work->N - 2; i >= 0; i--)
    {
        (solver->work->d.col(i)).noalias() = solver->cache->Quu_inv[i] * (solver->work->Bdyn.transpose() * solver->work->p.col(i + 1) + solver->work->r.col(i) + solver->cache->BPf);
        (solver->work->p.col(i)).noalias() = solver->work->q.col(i) + solver->cache->AmBKt[i].lazyProduct(solver->work->p.col(i + 1)) - (solver->cache->Kinf[i].transpose()).lazyProduct(solver->work->r.col(i)) + solver->cache->APf; 
    }
}

/**
 * Use LQR feedback policy to roll out trajectory
*/
void forward_pass(TinySolver *solver)
{
    for (int i = 0; i < solver->work->N - 1; i++)
    {
        (solver->work->u.col(i)).noalias() = -solver->cache->Kinf[i].lazyProduct(solver->work->x.col(i)) - solver->work->d.col(i);
        (solver->work->x.col(i + 1)).noalias() = solver->work->Adyn.lazyProduct(solver->work->x.col(i)) + solver->work->Bdyn.lazyProduct(solver->work->u.col(i)) + solver->work->fdyn;
        
        // double sci = 1;

        // primal update for var d
        // general format for d update d1 = 0.5 * (sd1x1 + sdx2 - ld1x1/rho - ld1x2/rho)
            
        /////// Should be correct and sufficient for ball + cone
        double coeff = 1.0 / (solver->cache->rhod1 * scix + solver->cache->rhod1 / scix + solver->cache->rhod1 / scix);
        solver->work->d1.col(i) = coeff * (solver->work->ld11.col(i) - solver->work->ld13.col(i) / scix - solver->work->ld16.col(i) / scix - solver->cache->rhod1 * solver->work->sd11.col(i)
        + solver->cache->rhod1 * solver->work->sd13.col(i) / scix + solver->cache->rhod1 * solver->work->sd16.col(i) / scix);
        solver->work->d1.col(i)(0) += solver->cache->rhod1 * a_max;

        coeff = 1.0 / (solver->cache->rhod2 * sciz + solver->cache->rhod2 / sciz + solver->cache->rhod2 / sciz);
        solver->work->d2.col(i) = coeff * (solver->work->ld22.col(i) - solver->work->ld24.col(i) / sciz - solver->work->ld29.col(i) / sciz - solver->cache->rhod2 * solver->work->sd22.col(i)
        + solver->cache->rhod2 * solver->work->sd24.col(i) / sciz + solver->cache->rhod2 * solver->work->sd29.col(i) / sciz);
        solver->work->d2.col(i)(0) += solver->cache->rhod2 * a_max;

        coeff = 1.0 / (solver->cache->rhod7 * sciy + solver->cache->rhod7 / sciy + solver->cache->rhod7 / sciy);
        solver->work->d7.col(i) = coeff * (solver->work->ld710.col(i) - solver->work->ld711.col(i) / sciy - solver->work->ld713.col(i) / sciy - solver->cache->rhod7 * solver->work->sd710.col(i)
        + solver->cache->rhod7 * solver->work->sd711.col(i) / sciy + solver->cache->rhod7 * solver->work->sd713.col(i) / sciy);
        solver->work->d7.col(i)(0) += solver->cache->rhod7 * a_max;

        solver->work->d3.col(i) = (1.0 / (2)) * (solver->work->sd35.col(i) + solver->work->sd36.col(i) - solver->work->ld35.col(i) / solver->cache->rhod3 - solver->work->ld36.col(i) / solver->cache->rhod3);        
        solver->work->d4.col(i) = (1.0 / (3)) * (solver->work->sd46.col(i) + solver->work->sd49.col(i) + solver->work->sd413.col(i) - solver->work->ld46.col(i) / solver->cache->rhod4 - solver->work->ld49.col(i) / solver->cache->rhod4 - solver->work->ld413.col(i) / solver->cache->rhod4);
        solver->work->d8.col(i) = (1.0 / (2)) * (solver->work->sd812.col(i) + solver->work->sd813.col(i) - solver->work->ld812.col(i) / solver->cache->rhod8 - solver->work->ld813.col(i) / solver->cache->rhod8);
    }
}


/**
 * Project slack (auxiliary) variables into their feasible domain, defined by
 * projection functions related to each constraint
 * TODO: pass in meta information with each constraint assigning it to a
 * projection function
*/
void update_slack(TinySolver *solver)
{

    double g = 9.81;
    // double sci = 1;
    
    // Update bound constraint slack variables for state
    solver->work->vnew = solver->work->x + solver->work->g;
    // Update bound constraint slack variables for input
    solver->work->znew = solver->work->u + solver->work->y;

    // slack update for gp
    // general format for slack update sd1x1 = d1 + ld1x1/rho
    solver->work->sd13 = solver->work->d1 + solver->work->ld13 / solver->cache->rhod1;
    solver->work->sd16 = solver->work->d1 + solver->work->ld16 / solver->cache->rhod1;
    solver->work->sd24 = solver->work->d2 + solver->work->ld24 / solver->cache->rhod2;
    solver->work->sd29 = solver->work->d2 + solver->work->ld29 / solver->cache->rhod2;
    solver->work->sd35 = solver->work->d3 + solver->work->ld35 / solver->cache->rhod3;
    solver->work->sd36 = solver->work->d3 + solver->work->ld36 / solver->cache->rhod3;
    solver->work->sd46 = solver->work->d4 + solver->work->ld46 / solver->cache->rhod4;
    solver->work->sd49 = solver->work->d4 + solver->work->ld49 / solver->cache->rhod4;
    solver->work->sd413 = solver->work->d4 + solver->work->ld413 / solver->cache->rhod4;
    solver->work->sd711 = solver->work->d7 + solver->work->ld711 / solver->cache->rhod7;
    solver->work->sd713 = solver->work->d7 + solver->work->ld713 / solver->cache->rhod7;
    solver->work->sd812 = solver->work->d8 + solver->work->ld812 / solver->cache->rhod8;
    solver->work->sd813 = solver->work->d8 + solver->work->ld813 / solver->cache->rhod8;
    solver->work->sz9 = solver->work->x + solver->work->lz9 / solver->cache->rhox;


    for (int k = 0; k < solver->work->N - 1; ++k) {

        solver->work->sd11.col(k) = -scix * solver->work->d1.col(k) + solver->work->ld11.col(k) / solver->cache->rhod1;
        solver->work->sd11.col(k)(0) += a_max;

        solver->work->sd22.col(k) = -sciz * solver->work->d2.col(k) + solver->work->ld22.col(k) / solver->cache->rhod2;
        solver->work->sd22.col(k)(0) += a_max; 

        solver->work->sd710.col(k) = -sciy * solver->work->d7.col(k) + solver->work->ld710.col(k) / solver->cache->rhod7;
        solver->work->sd710.col(k)(0) += a_max; 

        solver->work->sz1.col(k) = solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b 
        + solver->work->lz1.col(k) / solver->cache->rhox;

        solver->work->sz2.col(k) = solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b 
        + solver->work->lz2.col(k) / solver->cache->rhox;

        solver->work->sz10.col(k) = solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b 
        + solver->work->lz10.col(k) / solver->cache->rhox;

        solver->work->sz3.col(k) = solver->work->SC[k][1].A * solver->work->x.col(k) + solver->work->SC[k][1].b 
        + solver->work->lz3.col(k) / solver->cache->rhox;

        solver->work->sz4.col(k) = solver->work->SC[k][2].A * solver->work->x.col(k) + solver->work->SC[k][2].b 
        + solver->work->lz4.col(k) / solver->cache->rhox;

        solver->work->sz11.col(k) = solver->work->SC[k][4].A * solver->work->x.col(k) + solver->work->SC[k][4].b 
        + solver->work->lz11.col(k) / solver->cache->rhox;

        solver->work->sz5.col(k) = solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b 
        + solver->work->lz5.col(k) / solver->cache->rhox;

        solver->work->sz12.col(k) = solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b 
        + solver->work->lz12.col(k) / solver->cache->rhox;
    }

    if(solver->settings->en_state_bound == 5) {
        for (int k = 0; k < solver->work->N - 1; ++k) {
            // soc1
            solver->work->szt1.col(k) = solver->work->SC[k][5].A * solver->work->x.col(k) + solver->work->SC[k][5].b 
            + solver->work->lzt1.col(k) / solver->cache->rhox;
            
            solver->work->szt2.col(k) = (solver->work->SC[k][5].c.transpose() * solver->work->x.col(k)) 
            + solver->work->SC[k][5].d + solver->work->lzt2.col(k) / solver->cache->rhox;
        }
    }

    if (solver->settings->en_state_bound == 5) {
        for (int k = 0; k < solver->work->N - 1; ++k) {
            
            if(solver->work->szt1.col(k).norm() <= -solver->work->szt2.col(k)(0)) {
                std::cout << "below " << std::endl;
                solver->work->szt1.col(k).setZero();
                solver->work->szt2.col(k).setZero();
            }
            else if(solver->work->szt1.col(k).norm() <= solver->work->szt2.col(k)(0)) {
                solver->work->szt1.col(k) = solver->work->szt1.col(k);
                solver->work->szt2.col(k) = solver->work->szt2.col(k);
            }
            else if(solver->work->szt1.col(k).norm() > solver->work->szt2.col(k)(0)){
                // std::cout << "projecting " << std::endl;
                double norm = solver->work->szt1.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->szt2.col(k)(0) / norm);
                solver->work->szt1.col(k) = coeff * solver->work->szt1.col(k);
                solver->work->szt2.col(k)(0) = coeff * norm;
            }
        }
    }
    if (solver->settings->en_state_bound == 1) {

        for (int k = 0; k < solver->work->N - 1; ++k) {
            
            // C1: ||H_xyz Z + m_xyz|| <= a_max - sci * d1
            if(solver->work->sz1.col(k).norm() <= -solver->work->sd11.col(k)(0)) {
                // std::cout << "below 1" << std::endl;
                solver->work->sz1.col(k).setZero();
                solver->work->sd11.col(k).setZero();
            }
            else if(solver->work->sz1.col(k).norm() <= solver->work->sd11.col(k)(0)) {
                solver->work->sz1.col(k) = solver->work->sz1.col(k);
                solver->work->sd11.col(k) = solver->work->sd11.col(k);
            }
            else if(solver->work->sz1.col(k).norm() > solver->work->sd11.col(k)(0)){
                // std::cout << "projecting 1" << std::endl;
                double norm = solver->work->sz1.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd11.col(k)(0) / norm);
                solver->work->sz1.col(k) = coeff * solver->work->sz1.col(k);
                solver->work->sd11.col(k)(0) = coeff * norm;
            }

            // C2: ||H_xyz Z + m_xyz|| <= a_max - sci * d2
            if(solver->work->sz2.col(k).norm() <= -solver->work->sd22.col(k)(0)) {
                // std::cout << "below 2" << std::endl;
                solver->work->sz2.col(k).setZero();
                solver->work->sd22.col(k).setZero();
            }
            else if(solver->work->sz2.col(k).norm() <= solver->work->sd22.col(k)(0)) {
                solver->work->sz2.col(k) = solver->work->sz2.col(k);
                solver->work->sd22.col(k) = solver->work->sd22.col(k);
            }
            else if(solver->work->sz2.col(k).norm() > solver->work->sd22.col(k)(0)){
                // std::cout << "projecting 2" << std::endl;
                double norm = solver->work->sz2.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd22.col(k)(0) / norm);
                solver->work->sz2.col(k) = coeff * solver->work->sz2.col(k);
                solver->work->sd22.col(k)(0) = coeff * norm;
            }

             // C10: ||H_xyz Z + m_xyz|| <= a_max - sci * d7
            if(solver->work->sz10.col(k).norm() <= -solver->work->sd710.col(k)(0)) {
                // std::cout << "below 1" << std::endl;
                solver->work->sz10.col(k).setZero();
                solver->work->sd710.col(k).setZero();
            }
            else if(solver->work->sz10.col(k).norm() <= solver->work->sd710.col(k)(0)) {
                solver->work->sz10.col(k) = solver->work->sz10.col(k);
                solver->work->sd710.col(k) = solver->work->sd710.col(k);
            }
            else if(solver->work->sz10.col(k).norm() > solver->work->sd710.col(k)(0)){
                // std::cout << "projecting 1" << std::endl;
                double norm = solver->work->sz10.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd710.col(k)(0) / norm);
                solver->work->sz10.col(k) = coeff * solver->work->sz10.col(k);
                solver->work->sd710.col(k)(0) = coeff * norm;
            }

            // C3: ||L_x Z_bar|| <= d1
            if(solver->work->sz3.col(k).norm() <= -solver->work->sd13.col(k)(0)) {
                // std::cout << "below 3" << std::endl;
                solver->work->sz3.col(k).setZero();
                solver->work->sd13.col(k).setZero();
            }
            else if(solver->work->sz3.col(k).norm() <= solver->work->sd13.col(k)(0)) {
                solver->work->sz3.col(k) = solver->work->sz3.col(k);
                solver->work->sd13.col(k) = solver->work->sd13.col(k);
            }
            else if(solver->work->sz3.col(k).norm() > solver->work->sd13.col(k)(0)){
                // std::cout << "projecting 3" << std::endl;
                double norm = solver->work->sz3.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd13.col(k)(0) / norm);
                solver->work->sz3.col(k) = coeff * solver->work->sz3.col(k);
                solver->work->sd13.col(k)(0) = coeff * norm;
            }

            // C4: ||L_z Z_bar|| <= d2
            if(solver->work->sz4.col(k).norm() <= -solver->work->sd24.col(k)(0)) {
                // std::cout << "below 4" << std::endl;
                solver->work->sz4.col(k).setZero();
                solver->work->sd24.col(k).setZero();
            }
            else if(solver->work->sz4.col(k).norm() <= solver->work->sd24.col(k)(0)) {
                solver->work->sz4.col(k) = solver->work->sz4.col(k);
                solver->work->sd24.col(k) = solver->work->sd24.col(k);
            }
            else if(solver->work->sz4.col(k).norm() > solver->work->sd24.col(k)(0)){
                // std::cout << "projecting 4" << std::endl;
                double norm = solver->work->sz4.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd24.col(k)(0) / norm);
                solver->work->sz4.col(k) = coeff * solver->work->sz4.col(k);
                solver->work->sd24.col(k)(0) = coeff * norm;
            }

            // C11: ||L_y Z_bar|| <= d7
            if(solver->work->sz11.col(k).norm() <= -solver->work->sd711.col(k)(0)) {
                // std::cout << "below 4" << std::endl;
                solver->work->sz11.col(k).setZero();
                solver->work->sd711.col(k).setZero();
            }
            else if(solver->work->sz11.col(k).norm() <= solver->work->sd711.col(k)(0)) {
                solver->work->sz11.col(k) = solver->work->sz11.col(k);
                solver->work->sd711.col(k) = solver->work->sd711.col(k);
            }
            else if(solver->work->sz11.col(k).norm() > solver->work->sd711.col(k)(0)){
                // std::cout << "projecting 4" << std::endl;
                double norm = solver->work->sz11.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd711.col(k)(0) / norm);
                solver->work->sz11.col(k) = coeff * solver->work->sz11.col(k);
                solver->work->sd711.col(k)(0) = coeff * norm;
            }

            // cone constraint
            // C5: ||H_xy Z + m_xy|| <= d3
            if(solver->work->sz5.col(k).norm() <= -solver->work->sd35.col(k)(0)) {
                solver->work->sz5.col(k).setZero();
                solver->work->sd35.col(k).setZero();
            }
            else if(solver->work->sz5.col(k).norm() <= solver->work->sd35.col(k)(0)) {
                solver->work->sz5.col(k) = solver->work->sz5.col(k);
                solver->work->sd35.col(k) = solver->work->sd35.col(k);
            }
            else if(solver->work->sz5.col(k).norm() > solver->work->sd35.col(k)(0)){
                double norm = solver->work->sz5.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd35.col(k)(0) / norm);
                solver->work->sz5.col(k) = coeff * solver->work->sz5.col(k);
                solver->work->sd35.col(k)(0) = coeff * norm;
            }

            // C12: ||H_xy Z + m_xy|| <= d8
            if(solver->work->sz12.col(k).norm() <= -solver->work->sd812.col(k)(0)) {
                solver->work->sz12.col(k).setZero();
                solver->work->sd812.col(k).setZero();
            }
            else if(solver->work->sz12.col(k).norm() <= solver->work->sd812.col(k)(0)) {
                solver->work->sz12.col(k) = solver->work->sz12.col(k);
                solver->work->sd812.col(k) = solver->work->sd812.col(k);
            }
            else if(solver->work->sz12.col(k).norm() > solver->work->sd812.col(k)(0)){
                double norm = solver->work->sz12.col(k).norm();
                double coeff = 0.5 * (1 + solver->work->sd812.col(k)(0) / norm);
                solver->work->sz12.col(k) = coeff * solver->work->sz12.col(k);
                solver->work->sd812.col(k)(0) = coeff * norm;
            }

            // C6: d3 <= d4 - sci d1
            Eigen::VectorXd a(3);
            a << 1, -1, scix;
            if (solver->work->sd36.col(k)(0) - solver->work->sd46.col(k)(0) + scix * solver->work->sd16.col(k)(0) <= 0) {
                solver->work->sd36.col(k) = solver->work->sd36.col(k);
                solver->work->sd46.col(k) = solver->work->sd46.col(k);
                solver->work->sd16.col(k) = solver->work->sd16.col(k); 
            }
            else {
                double coeff = solver->work->sd36.col(k)(0) - solver->work->sd46.col(k)(0) + scix * solver->work->sd16.col(k)(0);
                solver->work->sd36.col(k)(0) = solver->work->sd36.col(k)(0) - coeff * a(0) / pow(a.norm(), 2);
                solver->work->sd46.col(k)(0) = solver->work->sd46.col(k)(0) - coeff * a(1) / pow(a.norm(), 2);
                solver->work->sd16.col(k)(0) = solver->work->sd16.col(k)(0) - coeff * a(2) / pow(a.norm(), 2);
            }

            // C13: d8 <= d4 - sci d7
            a << 1, -1, sciy;
            if (solver->work->sd813.col(k)(0) - solver->work->sd413.col(k)(0) + sciy * solver->work->sd713.col(k)(0) <= 0) {
                solver->work->sd813.col(k) = solver->work->sd813.col(k);
                solver->work->sd413.col(k) = solver->work->sd413.col(k);
                solver->work->sd713.col(k) = solver->work->sd713.col(k); 
            }
            else {
                double coeff = solver->work->sd813.col(k)(0) - solver->work->sd413.col(k)(0) + sciy * solver->work->sd713.col(k)(0);
                solver->work->sd813.col(k)(0) = solver->work->sd813.col(k)(0) - coeff * a(0) / pow(a.norm(), 2);
                solver->work->sd413.col(k)(0) = solver->work->sd413.col(k)(0) - coeff * a(1) / pow(a.norm(), 2);
                solver->work->sd713.col(k)(0) = solver->work->sd713.col(k)(0) - coeff * a(2) / pow(a.norm(), 2);
            }

            // C9: d2 <= c8 z + c9 d4
            Eigen::VectorXd a_9(12);
            a_9(0) = sciz;
            a_9.segment(1, 10) = -solver->work->H_z.col(k);
            a_9(11) = 1.0 / std::tan(theta);
            Eigen::VectorXd za1(12);
            za1(0) = solver->work->sd29.col(k)(0);
            za1.segment(1, 10) = solver->work->sz9.col(k);
            za1(11) = solver->work->sd49.col(k)(0);
            double b = solver->work->m_hat_z.col(k)(0);
            if (a_9.dot(za1) - b <= 0) {
                solver->work->sd29.col(k) = solver->work->sd29.col(k);
                solver->work->sz9.col(k) = solver->work->sz9.col(k);
                solver->work->sd49.col(k) = solver->work->sd49.col(k); 
            }
            else {
                // std::cout << "projecting 9" << std::endl;
                double coeff = a_9.dot(za1);
                solver->work->sd29.col(k)(0) = solver->work->sd29.col(k)(0) + (b - coeff) * a_9(0) / pow(a_9.norm(), 2);
                solver->work->sz9.col(k) = solver->work->sz9.col(k) + (b - coeff) * a_9.segment(1, 10) / pow(a_9.norm(), 2);
                solver->work->sd49.col(k)(0) = solver->work->sd49.col(k)(0) + (b - coeff) * a_9(11) / pow(a_9.norm(), 2);
            }
        }
    }

    // // Box constraints on input
    // if (solver->settings->en_input_bound) {
    //     // solver->work->znew = solver->work->u_max.cwiseMin(solver->work->u_min.cwiseMax(solver->work->znew));
    // }
}

/**
 * Update next iteration of dual variables by performing the augmented
 * lagrangian multiplier update
*/
void update_dual(TinySolver *solver)
{
    // double sci = 1;
    
    
    // dual update for const
    // general format for dual update ld1x1 = lx1x1 + rho(d1 - sd1x1)
    
    solver->work->ld13 = solver->work->ld13 + solver->cache->rhod1 * (solver->work->d1 - solver->work->sd13);
    solver->work->ld16 = solver->work->ld16 + solver->cache->rhod1 * (solver->work->d1 - solver->work->sd16);
    solver->work->ld24 = solver->work->ld24 + solver->cache->rhod2 * (solver->work->d2 - solver->work->sd24);
    solver->work->ld35 = solver->work->ld35 + solver->cache->rhod3 * (solver->work->d3 - solver->work->sd35);
    solver->work->ld36 = solver->work->ld36 + solver->cache->rhod3 * (solver->work->d3 - solver->work->sd36);
    solver->work->ld46 = solver->work->ld46 + solver->cache->rhod4 * (solver->work->d4 - solver->work->sd46);
    solver->work->ld49 = solver->work->ld49 + solver->cache->rhod4 * (solver->work->d4 - solver->work->sd49);
    solver->work->ld413 = solver->work->ld413 + solver->cache->rhod4 * (solver->work->d4 - solver->work->sd413);
    solver->work->ld711 = solver->work->ld711 + solver->cache->rhod7 * (solver->work->d7 - solver->work->sd711);
    solver->work->ld812 = solver->work->ld812 + solver->cache->rhod8 * (solver->work->d8 - solver->work->sd812);
    solver->work->ld813 = solver->work->ld813 + solver->cache->rhod8 * (solver->work->d8 - solver->work->sd813);
    solver->work->lz9 = solver->work->lz9 + solver->cache->rhox * (solver->work->x - solver->work->sz9);

     if (solver->settings->en_state_bound == 5) {
        for (int k = 0; k < solver->work->N - 1; ++k) {
            // soc1
            solver->work->lzt1.col(k) = solver->work->lzt1.col(k) + 
            solver->cache->rhox * (solver->work->SC[k][5].A * solver->work->x.col(k) + solver->work->SC[k][5].b - solver->work->szt1.col(k));
            
            solver->work->lzt2.col(k) = solver->work->lzt2.col(k) + 
            solver->cache->rhox * (solver->work->SC[k][5].c.transpose() * solver->work->x.col(k) + solver->work->SC[k][5].d - solver->work->szt2.col(k));
        }
     }
    
    for (int k = 0; k < solver->work->N - 1; ++k) {

        solver->work->ld11.col(k)(0) += solver->cache->rhod1 * 
        (a_max - scix * solver->work->d1.col(k)(0) - solver->work->sd11.col(k)(0));

        solver->work->ld22.col(k)(0) += solver->cache->rhod2 * 
        (a_max - sciz * solver->work->d2.col(k)(0) - solver->work->sd22.col(k)(0)); 

        solver->work->ld710.col(k)(0) += solver->cache->rhod7 * 
        (a_max - sciy * solver->work->d7.col(k)(0) - solver->work->sd710.col(k)(0));
        
        solver->work->lz1.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz1.col(k));
        
        solver->work->lz2.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz2.col(k));

        solver->work->lz10.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz10.col(k));

        solver->work->lz3.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][1].A * solver->work->x.col(k) + solver->work->SC[k][1].b - solver->work->sz3.col(k));

        solver->work->lz4.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][2].A * solver->work->x.col(k) + solver->work->SC[k][2].b - solver->work->sz4.col(k));

        solver->work->lz11.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][4].A * solver->work->x.col(k) + solver->work->SC[k][4].b - solver->work->sz11.col(k));

        solver->work->lz5.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b - solver->work->sz5.col(k));

        solver->work->lz12.col(k) +=
        solver->cache->rhox * (solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b - solver->work->sz12.col(k));
    }

    // Update bound constraint slack variables for state
    solver->work->g = solver->work->g + solver->work->x - solver->work->vnew;
    // Update bound constraint slack variables for input
    solver->work->y = solver->work->y + solver->work->u - solver->work->znew;
}

/**
 * Update linear control cost terms in the Riccati feedback using the changing
 * slack and dual variables from ADMM
*/
void update_linear_cost(TinySolver *solver)
{
    for (int i =0 ; i < solver->work->N ;i++){
        solver->work->q.col(i) = -(solver->work->Q[i] * solver->work->Xref.col(i));
        // ::Serial.printf("xrefr: %.6f\n", solver->work->Xref.col(i)(0)); 
    }
    if (solver->settings->en_state_bound == 0) {
        // (solver->work->q).noalias() -= solver->cache->rhox0 * (solver->work->vnew - solver->work->g);
    }
    if(solver->settings->en_state_bound == 5) {
        for (int k = 0; k < solver->work->N - 1; ++k) {
            (solver->work->q.col(k)).noalias() += solver->work->lzt1.col(k).transpose() * solver->work->SC[k][5].A +
            solver->cache->rhox * solver->work->SC[k][5].b.transpose() * solver->work->SC[k][5].A -
            solver->cache->rhox * solver->work->szt1.col(k).transpose() * solver->work->SC[k][5].A;
            
            (solver->work->q.col(k)).noalias() += solver->work->lzt2.col(k).transpose() * solver->work->SC[k][5].c.transpose() +
            solver->cache->rhox * solver->work->SC[k][5].d.transpose() * solver->work->SC[k][5].c.transpose() -
            solver->cache->rhox * solver->work->szt2.col(k).transpose() * solver->work->SC[k][5].c.transpose();
        }
    }
    if (solver->settings->en_state_bound == 1) {
        (solver->work->q).noalias() -= solver->cache->rhox * (solver->work->sz9 - solver->work->lz9);
        for (int k = 0; k < solver->work->N - 1; ++k) {
            (solver->work->q.col(k)).noalias() += solver->work->lz1.col(k).transpose() * solver->work->SC[k][0].A +
            solver->cache->rhox * solver->work->SC[k][0].b.transpose() * solver->work->SC[k][0].A -
            solver->cache->rhox * solver->work->sz1.col(k).transpose() * solver->work->SC[k][0].A;
            
            (solver->work->q.col(k)).noalias() += solver->work->lz2.col(k).transpose() * solver->work->SC[k][0].A +
            solver->cache->rhox * solver->work->SC[k][0].b.transpose() * solver->work->SC[k][0].A -
            solver->cache->rhox * solver->work->sz2.col(k).transpose() * solver->work->SC[k][0].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz10.col(k).transpose() * solver->work->SC[k][0].A +
            solver->cache->rhox * solver->work->SC[k][0].b.transpose() * solver->work->SC[k][0].A -
            solver->cache->rhox * solver->work->sz10.col(k).transpose() * solver->work->SC[k][0].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz3.col(k).transpose() * solver->work->SC[k][1].A +
            solver->cache->rhox * solver->work->SC[k][1].b.transpose() * solver->work->SC[k][1].A -
            solver->cache->rhox * solver->work->sz3.col(k).transpose() * solver->work->SC[k][1].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz4.col(k).transpose() * solver->work->SC[k][2].A +
            solver->cache->rhox * solver->work->SC[k][2].b.transpose() * solver->work->SC[k][2].A -
            solver->cache->rhox * solver->work->sz4.col(k).transpose() * solver->work->SC[k][2].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz11.col(k).transpose() * solver->work->SC[k][4].A +
            solver->cache->rhox * solver->work->SC[k][4].b.transpose() * solver->work->SC[k][4].A -
            solver->cache->rhox * solver->work->sz11.col(k).transpose() * solver->work->SC[k][4].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz5.col(k).transpose() * solver->work->SC[k][3].A +
            solver->cache->rhox * solver->work->SC[k][3].b.transpose() * solver->work->SC[k][3].A -
            solver->cache->rhox * solver->work->sz5.col(k).transpose() * solver->work->SC[k][3].A;

            (solver->work->q.col(k)).noalias() += solver->work->lz12.col(k).transpose() * solver->work->SC[k][3].A +
            solver->cache->rhox * solver->work->SC[k][3].b.transpose() * solver->work->SC[k][3].A -
            solver->cache->rhox * solver->work->sz12.col(k).transpose() * solver->work->SC[k][3].A;
        }
    }

    solver->work->r = -(solver->work->Uref.array().colwise() * solver->work->R.array());
    // (solver->work->r).noalias() -= solver->cache->rhou * (solver->work->znew - solver->work->y);

    solver->work->p.col(solver->work->N - 1) = -(solver->work->Xref.col(solver->work->N - 1).transpose().lazyProduct(solver->cache->Pinf[solver->work->N-1]));
    if (solver->settings->en_state_bound == 0) {
        // (solver->work->p.col(solver->work->N - 1)).noalias() -= solver->cache->rhox0 * (solver->work->vnew.col(solver->work->N - 1) - solver->work->g.col(solver->work->N - 1));
    }
    if (solver->settings->en_state_bound == 5) {
        (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lzt1.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][5].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][5].b.transpose() * solver->work->SC[solver->work->N - 1][5].A -
            solver->cache->rhox * solver->work->szt1.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][5].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lzt2.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][5].c.transpose() +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][5].d.transpose() * solver->work->SC[solver->work->N - 1][5].c.transpose() -
            solver->cache->rhox * solver->work->szt2.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][5].c.transpose();
    }
    if (solver->settings->en_state_bound == 1) {
            (solver->work->p.col(solver->work->N - 1)).noalias() -= solver->cache->rhox * (solver->work->sz9.col(solver->work->N - 1) - solver->work->lz9.col(solver->work->N - 1));            
            
            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz1.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][0].b.transpose() * solver->work->SC[solver->work->N - 1][0].A -
            solver->cache->rhox * solver->work->sz1.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz2.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][0].b.transpose() * solver->work->SC[solver->work->N - 1][0].A -
            solver->cache->rhox * solver->work->sz2.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz10.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][0].b.transpose() * solver->work->SC[solver->work->N - 1][0].A -
            solver->cache->rhox * solver->work->sz10.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][0].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz3.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][1].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][1].b.transpose() * solver->work->SC[solver->work->N - 1][1].A -
            solver->cache->rhox * solver->work->sz3.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][1].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz4.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][2].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][2].b.transpose() * solver->work->SC[solver->work->N - 1][2].A -
            solver->cache->rhox * solver->work->sz4.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][2].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz11.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][4].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][4].b.transpose() * solver->work->SC[solver->work->N - 1][4].A -
            solver->cache->rhox * solver->work->sz11.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][4].A;

            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz5.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][3].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][3].b.transpose() * solver->work->SC[solver->work->N - 1][3].A -
            solver->cache->rhox * solver->work->sz5.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][3].A;
            
            (solver->work->p.col(solver->work->N - 1)).noalias() += solver->work->lz12.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][3].A +
            solver->cache->rhox * solver->work->SC[solver->work->N - 1][3].b.transpose() * solver->work->SC[solver->work->N - 1][3].A -
            solver->cache->rhox * solver->work->sz12.col(solver->work->N - 1).transpose() * solver->work->SC[solver->work->N - 1][3].A;
    }
}

/**
 * Check for termination condition by evaluating whether the largest absolute
 * primal and dual residuals for states and inputs are below threhold.
*/
bool termination_condition(TinySolver *solver)
{   
    // double sci = 1;
    if (solver->work->iter % solver->settings->check_termination == 0)
    {
        if (solver->settings->en_state_bound == 0) {
            solver->work->primal_residual_state = (solver->work->x - solver->work->vnew).cwiseAbs().maxCoeff();
        }
        
        solver->work->primal_residual_input = (solver->work->u - solver->work->znew).cwiseAbs().maxCoeff();
        if (solver->settings->en_state_bound == 5) {
            for (int k = 0; k < solver->work->N - 1; ++k) {
                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][5].A * solver->work->x.col(k) + solver->work->SC[k][5].b - solver->work->szt1.col(k)).cwiseAbs().maxCoeff());
                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, 
                (solver->work->SC[k][5].c.transpose() * solver->work->x.col(k) + solver->work->SC[k][5].d - solver->work->szt2.col(k)).cwiseAbs().maxCoeff());
            }
        }
        if (solver->settings->en_state_bound == 1) {
            for (int k = 0; k < solver->work->N - 1; ++k) {
                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz1.col(k)).cwiseAbs().maxCoeff());
                
                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz2.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][0].A * solver->work->x.col(k) + solver->work->SC[k][0].b - solver->work->sz10.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][1].A * solver->work->x.col(k) + solver->work->SC[k][1].b - solver->work->sz3.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][2].A * solver->work->x.col(k) + solver->work->SC[k][2].b - solver->work->sz4.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][4].A * solver->work->x.col(k) + solver->work->SC[k][4].b - solver->work->sz11.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b - solver->work->sz5.col(k)).cwiseAbs().maxCoeff());

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state,
                (solver->work->SC[k][3].A * solver->work->x.col(k) + solver->work->SC[k][3].b - solver->work->sz12.col(k)).cwiseAbs().maxCoeff());
                
                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, 
                std::abs((a_max - scix * solver->work->d1.col(k)(0) - solver->work->sd11.col(k)(0))));

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, 
                std::abs((a_max - sciz * solver->work->d2.col(k)(0) - solver->work->sd22.col(k)(0))));

                solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, 
                std::abs((a_max - sciy * solver->work->d7.col(k)(0) - solver->work->sd710.col(k)(0))));
            }
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->x - solver->work->sz9).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d1 - solver->work->sd13).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d1 - solver->work->sd16).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d2 - solver->work->sd24).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d2 - solver->work->sd29).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d3 - solver->work->sd35).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d3 - solver->work->sd36).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d4 - solver->work->sd46).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d4 - solver->work->sd49).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d4 - solver->work->sd413).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d7 - solver->work->sd711).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d7 - solver->work->sd713).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d8 - solver->work->sd812).cwiseAbs().maxCoeff());
            solver->work->primal_residual_state = std::max(solver->work->primal_residual_state, (solver->work->d8 - solver->work->sd813).cwiseAbs().maxCoeff());
        }
        if (solver->settings->en_state_bound == 0) {
            solver->work->dual_residual_state = ((solver->work->v - solver->work->vnew).cwiseAbs().maxCoeff()) * solver->cache->rhox0;
        }
        if (solver->settings->en_state_bound == 5) {
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->szt1_p - solver->work->szt1).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->szt2_p - solver->work->szt2).cwiseAbs().maxCoeff()) * solver->cache->rhox);
        }
        if (solver->settings->en_state_bound == 1) {
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz1_p - solver->work->sz1).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz2_p - solver->work->sz2).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz3_p - solver->work->sz3).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz4_p - solver->work->sz4).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz5_p - solver->work->sz5).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz9_p - solver->work->sz9).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz10_p - solver->work->sz10).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz11_p - solver->work->sz11).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sz12_p - solver->work->sz12).cwiseAbs().maxCoeff()) * solver->cache->rhox);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd11_p - solver->work->sd11).cwiseAbs().maxCoeff()) * solver->cache->rhod1);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd13_p - solver->work->sd13).cwiseAbs().maxCoeff()) * solver->cache->rhod1);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd16_p - solver->work->sd16).cwiseAbs().maxCoeff()) * solver->cache->rhod1);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd22_p - solver->work->sd22).cwiseAbs().maxCoeff()) * solver->cache->rhod2);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd24_p - solver->work->sd24).cwiseAbs().maxCoeff()) * solver->cache->rhod2);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd29_p - solver->work->sd29).cwiseAbs().maxCoeff()) * solver->cache->rhod2);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd35_p - solver->work->sd35).cwiseAbs().maxCoeff()) * solver->cache->rhod3);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd36_p - solver->work->sd36).cwiseAbs().maxCoeff()) * solver->cache->rhod3);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd46_p - solver->work->sd46).cwiseAbs().maxCoeff()) * solver->cache->rhod4);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd49_p - solver->work->sd49).cwiseAbs().maxCoeff()) * solver->cache->rhod4);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd413_p - solver->work->sd413).cwiseAbs().maxCoeff()) * solver->cache->rhod4);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd710_p - solver->work->sd710).cwiseAbs().maxCoeff()) * solver->cache->rhod7);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd711_p - solver->work->sd711).cwiseAbs().maxCoeff()) * solver->cache->rhod7);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd713_p - solver->work->sd713).cwiseAbs().maxCoeff()) * solver->cache->rhod7);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd812_p - solver->work->sd812).cwiseAbs().maxCoeff()) * solver->cache->rhod8);
            solver->work->dual_residual_state = std::max(solver->work->dual_residual_state, ((solver->work->sd813_p - solver->work->sd813).cwiseAbs().maxCoeff()) * solver->cache->rhod8);
        }
        solver->work->dual_residual_input = ((solver->work->z - solver->work->znew).cwiseAbs().maxCoeff()) * solver->cache->rhou;
        
        // std::cout << "prim res: " << solver->work->primal_residual_state << std::endl;
        // std::cout << "dual res: " << solver->work->dual_residual_state << std::endl; 
        if (solver->work->primal_residual_state < solver->settings->abs_pri_tol &&
            solver->work->primal_residual_input < solver->settings->abs_pri_tol &&
            solver->work->dual_residual_state < solver->settings->abs_dua_tol &&
            solver->work->dual_residual_input < solver->settings->abs_dua_tol)
        {
            return true;                 
        }
    }
    return false;
}

int solve(TinySolver *solver)
{
    // Initialize variables
    solver->solution->solved = 0;
    solver->solution->iter = 0;
    solver->work->status = 11; // TINY_UNSOLVED
    solver->work->iter = 0;

    for (int i = 0; i < solver->settings->max_iter; i++)
    {
        // Solve linear system with Riccati and roll out to get new trajectory
        forward_pass(solver);

        // Project slack variables into feasible domain
        update_slack(solver);

        // Compute next iteration of dual variables
        update_dual(solver);

        // Update linear control cost terms using reference trajectory, duals, and slack variables
        update_linear_cost(solver);

        solver->work->iter += 1;

        // // Check for whether cost is minimized by calculating residuals
        // if (termination_condition(solver)) {
        //     solver->work->status = 1; // TINY_SOLVED

        //     // Save solution
        //     solver->solution->iter = solver->work->iter;
        //     solver->solution->solved = 1;
        //     solver->solution->x = solver->work->vnew;
        //     solver->solution->u = solver->work->znew;
        //     return 0;
        // }

        // Save previous slack variables
        solver->work->v = solver->work->vnew;
        solver->work->z = solver->work->znew;


        solver->work->sz1_p = solver->work->sz1;
        solver->work->sz2_p = solver->work->sz2;
        solver->work->sz3_p = solver->work->sz3;
        solver->work->sz4_p = solver->work->sz4;
        solver->work->sz5_p = solver->work->sz5;
        solver->work->sz9_p = solver->work->sz9;
        solver->work->sz10_p = solver->work->sz10;
        solver->work->sz11_p = solver->work->sz11;
        solver->work->sz12_p = solver->work->sz12;
        solver->work->szt1_p = solver->work->szt1;
        solver->work->szt2_p = solver->work->szt2;
        solver->work->sd11_p = solver->work->sd11;
        solver->work->sd13_p = solver->work->sd13;
        solver->work->sd16_p = solver->work->sd16;
        solver->work->sd22_p = solver->work->sd22;
        solver->work->sd24_p = solver->work->sd24;
        solver->work->sd29_p = solver->work->sd29;
        solver->work->sd35_p = solver->work->sd35;
        solver->work->sd36_p = solver->work->sd36;
        solver->work->sd46_p = solver->work->sd46;
        solver->work->sd49_p = solver->work->sd49;
        solver->work->sd413_p = solver->work->sd413;
        solver->work->sd710_p = solver->work->sd710;
        solver->work->sd711_p = solver->work->sd711;
        solver->work->sd713_p = solver->work->sd713;
        solver->work->sd812_p = solver->work->sd812;
        solver->work->sd813_p = solver->work->sd813;

        backward_pass_grad(solver);

    }
    solver->solution->iter = solver->work->iter;
    solver->solution->solved = 0;
    solver->solution->x = solver->work->vnew;
    solver->solution->u = solver->work->znew;
    return 1;
}

} /* extern "C" */
