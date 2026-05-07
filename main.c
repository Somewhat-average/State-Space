#include <stdio.h>
#include "state_space.h"

int main(void)
{
    stateSpace * sys;

    double k = 1;
    double a = 0;
    double m = 10;
    double x0 = 5;
    double v0 = 0;
    double t = 0;
    double t_stop = 100;

    const int dimx = 2;
    const int dimu = 1;
    const int dimy = 1;

    sys = stateSpaceCreate(dimx, dimy, dimu);
    if (!sys) { return -1; }

    /* set up state space eqns */
    setElement(sys->A, 1, 2, 1);
    setElement(sys->A, 2, 1, -k/m);
    setElement(sys->A, 2, 2, -a/m);

    setElement(sys->B, 2, 1, 1/m);
    setElement(sys->C, 1, 1, 1);

    setElement(sys->x, 1, 1, x0);
    setElement(sys->x, 2, 1, v0);

    /* apply the inital_state */
    calculateRungeKutta(sys);
    stateSpaceUpdate(sys);

    while (t <= t_stop) {
        printf("%.4f,%.4f\n", t, ELEM(sys->y, 1, 1));
        stateSpaceUpdate(sys);
        t += sys->h;
    }

    return stateSpaceDestroy(sys);
}
