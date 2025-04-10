#pragma once


#include <Eigen/Eigen.h>
#include <vector>
// #include <Eigen/Core>
// #include <Eigen/LU>

using namespace Eigen;


#ifdef __cplusplus
extern "C"
{
#endif

    typedef double tinytype;  // should be double if you want to generate code
    typedef Matrix<tinytype, Dynamic, Dynamic> tinyMatrix;
    typedef Matrix<tinytype, Dynamic, 1> tinyVector;

    typedef struct {
      tinyMatrix A;
      tinyVector b;
      tinyVector c;
      tinyVector d;
    } SOC;
    // typedef Matrix<tinytype, NSTATES, 1> tiny_VectorNx;
    // typedef Matrix<tinytype, NINPUTS, 1> tiny_VectorNu;
    // typedef Matrix<tinytype, NSTATES, NSTATES> tiny_MatrixNxNx;
    // typedef Matrix<tinytype, NSTATES, NINPUTS> tiny_MatrixNxNu;
    // typedef Matrix<tinytype, NINPUTS, NSTATES> tiny_MatrixNuNx;
    // typedef Matrix<tinytype, NINPUTS, NINPUTS> tiny_MatrixNuNu;

    // typedef Matrix<tinytype, NSTATES, NHORIZON> tiny_MatrixNxNh;       // Nu x Nh
    // typedef Matrix<tinytype, NINPUTS, NHORIZON - 1> tiny_MatrixNuNhm1; // Nu x Nh-1

    /**
     * Solution
     */
    typedef struct {
        int iter;
        int solved;
        tinyMatrix x; // nx x N
        tinyMatrix u; // nu x N-1
    } TinySolution;

    /**
     * Matrices that must be recomputed with changes in time step, rho
     */
    typedef struct {
        tinytype rhou;
        tinytype rhox;
        tinytype rhox0;

        tinytype rhod1;
        tinytype rhod2;
        tinytype rhod3;
        tinytype rhod4;
        tinytype rhod7;
        tinytype rhod8;

        tinyMatrix Kinf[10];       // nu x nx
        tinyMatrix Pinf[10];       // nx x nx
        tinyMatrix Quu_inv[10];    // nu x nu
        tinyMatrix AmBKt[10];      // nx x nx
        tinyVector APf;        // nx x 1
        tinyVector BPf;        // nu x 1
    } TinyCache;

    /**
     * User settings
     */
    typedef struct {
        tinytype abs_pri_tol;
        tinytype abs_dua_tol;
        int max_iter;
        int check_termination;
        int en_state_bound;
        int en_input_bound;
    } TinySettings;

    /**
     * Problem variables
     */
    typedef struct {
        int nx; // Number of states
        int nu; // Number of control inputs
        int N;  // Number of knotpoints in the horizon

        tinyMatrix z_star;

        tinyMatrix m_hat_z;
        tinyMatrix H_z;

        // State and input
        tinyMatrix x;    // nx x N
        tinyMatrix u;    // nu x N-1
        tinyMatrix d1;    // nu x N-1
        tinyMatrix d2;    // nu x N-1
        tinyMatrix d3;    // nu x N-1
        tinyMatrix d4;    // nu x N-1
        tinyMatrix d7;    // nu x N-1
        tinyMatrix d8;    // nu x N-1

        // Linear control cost terms
        tinyMatrix q;    // nx x N
        tinyMatrix r;    // nu x N-1

        // Linear Riccati backward pass terms
        tinyMatrix p;    // nx x N
        tinyMatrix d;    // nu x N-1

        // Auxiliary variables
        tinyMatrix v;    // nx x N

        // prev slack vars for const 3 4
        tinyMatrix sd11_p;
        tinyMatrix sd13_p;
        tinyMatrix sd16_p;
        tinyMatrix sd22_p;
        tinyMatrix sd24_p;
        tinyMatrix sd29_p;
        tinyMatrix sd35_p;
        tinyMatrix sd36_p;
        tinyMatrix sd46_p;
        tinyMatrix sd49_p;
        tinyMatrix sd413_p;
        tinyMatrix sd710_p;
        tinyMatrix sd711_p;
        tinyMatrix sd713_p;
        tinyMatrix sd812_p;
        tinyMatrix sd813_p;
        tinyMatrix sz1_p;
        tinyMatrix sz2_p;
        tinyMatrix sz3_p;
        tinyMatrix sz4_p;
        tinyMatrix sz5_p;
        tinyMatrix sz9_p;
        tinyMatrix sz10_p;
        tinyMatrix sz11_p;
        tinyMatrix sz12_p;
        tinyMatrix szt1_p;
        tinyMatrix szt2_p;

        // slack vars for const 3 4
        tinyMatrix sd11;
        tinyMatrix sd13;
        tinyMatrix sd16;
        tinyMatrix sd22;
        tinyMatrix sd24;
        tinyMatrix sd29;
        tinyMatrix sd35;
        tinyMatrix sd36;
        tinyMatrix sd46;
        tinyMatrix sd49;
        tinyMatrix sd413;
        tinyMatrix sd710;
        tinyMatrix sd711;
        tinyMatrix sd713;
        tinyMatrix sd812;
        tinyMatrix sd813;

        tinyMatrix sz1;
        tinyMatrix sz2;
        tinyMatrix sz3;
        tinyMatrix sz4;
        tinyMatrix sz5;
        tinyMatrix sz9;
        tinyMatrix sz10;
        tinyMatrix sz11;
        tinyMatrix sz12;
        tinyMatrix szt1;
        tinyMatrix szt2;

        
        tinyMatrix vnew; // nx x N
        tinyMatrix z;    // nu x N-1
        tinyMatrix znew; // nu x N-1
        tinyMatrix vcnew; // nx x N
        tinyMatrix zcnew; // nu x N-1
        // Dual variables

        // dual vars for const 3 4
        tinyMatrix ld11;
        tinyMatrix ld13;
        tinyMatrix ld16;
        tinyMatrix ld22;
        tinyMatrix ld24;
        tinyMatrix ld29;
        tinyMatrix ld35;
        tinyMatrix ld36;
        tinyMatrix ld46;
        tinyMatrix ld49;
        tinyMatrix ld413;
        tinyMatrix ld710;
        tinyMatrix ld711;
        tinyMatrix ld713;
        tinyMatrix ld812;
        tinyMatrix ld813;
        tinyMatrix lz1;
        tinyMatrix lz2;
        tinyMatrix lz3;
        tinyMatrix lz4;
        tinyMatrix lz5;
        tinyMatrix lz9;
        tinyMatrix lz10;
        tinyMatrix lz11;
        tinyMatrix lz12;
        tinyMatrix lzt1;
        tinyMatrix lzt2;
        
        tinyMatrix g;    // nx x N
        tinyMatrix y;    // nu x N-1

        // Q, R, A, B given by user
        tinyMatrix Q[10];       // nx x 1
        tinyVector R;       // nu x 1
        tinyMatrix Adyn;    // nx x nx
        tinyMatrix Bdyn;    // nx x nu
        tinyVector fdyn;    // nx x 1 affine vector

        // State and input bounds
        tinyMatrix x_min;   // nx x N
        tinyMatrix x_max;   // nx x N
        tinyMatrix u_min;   // nu x N-1
        tinyMatrix u_max;   // nu x N-1

        // Reference trajectory to track for one horizon
        tinyMatrix Xref;    // nx x N
        tinyMatrix Uref;    // nu x N-1

        // Temporaries
        tinyVector Qu;      // nu x 1

        // State SOC
        int numStateCones;
        std::vector <std::vector<SOC>> SC;
        std::vector<SOC> socx;
        // Input SOC
        int numInputCones;
        std::vector<SOC> socu;

        // Variables for keeping track of solve status
        tinytype primal_residual_state;
        tinytype primal_residual_input;
        tinytype dual_residual_state;
        tinytype dual_residual_input;
        int status;
        int iter;
    } TinyWorkspace;

    /**
     * Main TinyMPC solver structure that holds all information.
     */
    typedef struct {
        TinySolution *solution; // Solution
        TinySettings *settings; // Problem settings
        TinyCache *cache;       // Problem cache
        TinyWorkspace *work;    // Solver workspace
    } TinySolver;

#ifdef __cplusplus
}
#endif
