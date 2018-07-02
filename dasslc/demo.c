#include "dasslc.h"

// The residual function declaration
DASSLC_RES residuals;

void main(){
    // The root structure
    PTR_ROOT root;

    // The integration interval (time span)
    REAL t0 = 0.0, dt = 0.1, tf = 1, tout = dt;

    // The configuration file
    char *inputfile = "demo.dat";

    // The system's rank (number of equations)
    int rank = 2;

    // The initial conditions
    REAL y0[2];
    y0[0] = y0[1] = 1;

    // The setup call
    BOOL error = daSetup(inputfile, &root, residuals, rank, t0, y0, NULL, NULL, NULL, NULL);
    if (error){
        printf ("Setup error = %d\n", error);
        return;
    }

    // Find the initial conditions (for the derivatives)
    error = dasslc(INITIAL_COND, &root, residuals, &t0, tout, NULL, NULL);
    if (error < 0){
        printf ("error = %d\n", error);
        return;
    }

    // Integrate over the time span
    do{
        error = dasslc(TRANSIENT, &root, residuals, &t0, tout, NULL, NULL);
        if (error < 0){
            printf ("error = %d\n", error);
            return;
        }
        tout += dt;
    }while(tout < tf);

    // Solve for the steady state case
    if (error >= 0)
        error = dasslc(STEADY_STATE, &root, residuals, &t0, tf, NULL, NULL);
    if (error < 0)
        printf ("error = %d\n", error);

    // Save the steady state in file
    daStat(root.savefile, &root);
    daFree (&root);

    return;

} // end main

BOOL residuals (PTR_ROOT *root, REAL t, REAL *y, REAL *yp, REAL *res, BOOL *jac){
    BOOL error = FALSE;
   
    res[0] = yp[0] + 2*y[0];
    res[1] = yp[1] + y[1];

    return (error);
} // end residuals

