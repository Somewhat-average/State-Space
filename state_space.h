#include "matrix.h"

typedef struct {
    int size_x;
    int size_y;
    int size_u;

    /* state transistion */
    matrix *A;
    matrix *B;
    matrix *C;
    matrix *D;

    /* state */
    matrix *x;
    matrix *x_dot;
    matrix *u;
    matrix *y;

    /* helpers */
    matrix *vec_x1;
    matrix *vec_x2;
    matrix *vec_x3;
    matrix *vec_y1;
    matrix *vec_y2;
    matrix *k1;
    matrix *k2;
    matrix *k3;
    matrix *k4;
    matrix *A1;
    matrix *A2;
    matrix *A3;
    matrix *A4;
    matrix *RK4_A;
    matrix *RK4_B;
    matrix *I;

    /* time step */
    double h;
} stateSpace;

stateSpace * stateSpaceCreate(const int size_x, const int size_u, const int size_y);
int stateSpaceDestroy(stateSpace * ss);
void stateSpaceUpdate(stateSpace * ss);
void calculateRungeKutta(stateSpace * ss);


stateSpace * stateSpaceCreate(const int size_x, const int size_u, const int size_y) {
    stateSpace * ss;

    /* must have some size */
    if (size_x < 1) { return NULL; }
    if (size_u < 1) { return NULL; }
    if (size_y < 1) { return NULL; }

    ss = (stateSpace *)malloc(sizeof(stateSpace));
    if (ss == NULL) { return NULL; }

    ss->size_x = size_x;
    ss->size_y = size_y;
    ss->size_u = size_u;

    ss->A = NULL;
    ss->B = NULL;
    ss->C = NULL;
    ss->D = NULL;
    ss->x = NULL;
    ss->x_dot = NULL;
    ss->u = NULL;
    ss->y = NULL;
    ss->vec_x1 = NULL;
    ss->vec_x2 = NULL;
    ss->vec_x3 = NULL;
    ss->vec_y1 = NULL;
    ss->vec_y2 = NULL;
    ss->k1 = NULL;
    ss->k2 = NULL;
    ss->k3 = NULL;
    ss->k4 = NULL;
    ss->A1 = NULL;
    ss->A2 = NULL;
    ss->A3 = NULL;
    ss->A4 = NULL;
    ss->RK4_A = NULL;
    ss->RK4_B = NULL;
    ss->I = NULL;

    /* rows = output len; cols = vector len */
    ss->A = matrixCreate(size_x, size_x);
    ss->B = matrixCreate(size_x, size_u);
    ss->C = matrixCreate(size_y, size_x);
    ss->D = matrixCreate(size_y, size_u);

    /* vectors */
    ss->x = matrixCreate(size_x, 1);
    ss->x_dot = matrixCreate(size_x, 1);
    ss->y = matrixCreate(size_y, 1);
    ss->u = matrixCreate(size_u, 1);

    /* helpers */
    ss->vec_x1 = matrixCreate(size_x, 1);
    ss->vec_x2 = matrixCreate(size_x, 1);
    ss->vec_x3 = matrixCreate(size_x, 1);
    ss->vec_y1 = matrixCreate(size_y, 1);
    ss->vec_y2 = matrixCreate(size_y, 1);
    ss->k1 = matrixCreate(size_x, 1);
    ss->k2 = matrixCreate(size_x, 1);
    ss->k3 = matrixCreate(size_x, 1);
    ss->k4 = matrixCreate(size_x, 1);
    ss->A1 = matrixCreate(size_x, size_x);
    ss->A2 = matrixCreate(size_x, size_x);
    ss->A3 = matrixCreate(size_x, size_x);
    ss->A4 = matrixCreate(size_x, size_x);
    ss->RK4_A = matrixCreate(size_x, size_x);
    ss->RK4_B = matrixCreate(size_x, size_x);

    ss->I = matrixCreate(size_x, size_x);
    identity(ss->I);

    if (!ss->A ||
        !ss->B ||
        !ss->C ||
        !ss->D ||
        !ss->x ||
        !ss->y ||
        !ss->u ||
        !ss->x_dot ||
        !ss->vec_x1 ||
        !ss->vec_x2 ||
        !ss->vec_x3 ||
        !ss->vec_y1 ||
        !ss->vec_y2 ||
        !ss->k1 ||
        !ss->k2 ||
        !ss->k3 ||
        !ss->k4 ||
        !ss->A1 ||
        !ss->A2 ||
        !ss->A3 ||
        !ss->A4 ||
        !ss->RK4_A ||
        !ss->RK4_B ||
        !ss->I) {
        stateSpaceDestroy(ss);
        return NULL;
    }

    /* default step size */
    ss->h = 1e-2;
    
    return ss;
}


int stateSpaceDestroy(stateSpace *ss) {
    int failed = 0;
    if (ss == NULL) { return 0; }
    if (ss->A)
        if (matrixDestroy(ss->A)) { failed = -1; }
    if (ss->B)
        if (matrixDestroy(ss->B)) { failed = -2; }
    if (ss->C)
        if (matrixDestroy(ss->C)) { failed = -3; }
    if (ss->D)
        if (matrixDestroy(ss->D)) { failed = -4; }
    if (ss->x)
        if (matrixDestroy(ss->x)) { failed = -5; }
    if (ss->x_dot)
        if (matrixDestroy(ss->x_dot)) { failed = -6; }
    if (ss->u)
        if (matrixDestroy(ss->u)) { failed = -7; }
    if (ss->y)
        if (matrixDestroy(ss->y)) { failed = -8; }
    if (ss->vec_x1)
        if (matrixDestroy(ss->vec_x1)) { failed = -9; }
    if (ss->vec_x2)
        if (matrixDestroy(ss->vec_x2)) { failed = -10; }
    if (ss->vec_x3)
        if (matrixDestroy(ss->vec_x3)) { failed = -11; }
    if (ss->vec_y1)
        if (matrixDestroy(ss->vec_y1)) { failed = -12; }
    if (ss->vec_y2)
        if (matrixDestroy(ss->vec_y2)) { failed = -13; }
    if (ss->k1)
        if (matrixDestroy(ss->k1)) { failed = -14; }
    if (ss->k2)
        if (matrixDestroy(ss->k2)) { failed = -15; }
    if (ss->k3)
        if (matrixDestroy(ss->k3)) { failed = -16; }
    if (ss->k4)
        if (matrixDestroy(ss->k4)) { failed = -17; }
    if (ss->A1)
        if (matrixDestroy(ss->A1)) { failed = -18; }
    if (ss->A2)
        if (matrixDestroy(ss->A2)) { failed = -19; }
    if (ss->A3)
        if (matrixDestroy(ss->A3)) { failed = -20; }
    if (ss->A4)
        if (matrixDestroy(ss->A4)) { failed = -21; }
    if (ss->RK4_A)
        if (matrixDestroy(ss->RK4_A)) { failed = -22; }
    if (ss->RK4_B)
        if (matrixDestroy(ss->RK4_B)) { failed = -23; }
    if (ss->I)
        if (matrixDestroy(ss->I)) { failed = -24; }

    free(ss);
    return failed; 
}

/* advance the system by h seconds */
void stateSpaceUpdate(stateSpace * ss) {
    /* update y */
    product(ss->C, ss->x, ss->vec_y1);
    product(ss->D, ss->u, ss->vec_y2);
    sum(ss->vec_y1, ss->vec_y2, ss->y);

    /* prepare to update x */
    product(ss->RK4_A, ss->x, ss->vec_x1);
    product(ss->B, ss->u, ss->vec_x2);
    product(ss->RK4_B, ss->vec_x2, ss->vec_x3);

    /* store old data in unused matrix */
    matrixCopyData(ss->x, ss->vec_x2);
    /* update x_dot (unneccesary) */
    difference(ss->x, ss->vec_x2, ss->x_dot);

    /* update x */
    sum(ss->vec_x1, ss->vec_x3, ss->x);

    return;
}


void calculateRungeKutta(stateSpace * ss) {
    int row, col;

    double c1 = ss->h;
    double c2 = c1*c1/2.0;
    double c3 = c1*c1*c1/6.0;
    double c4 = c1*c1*c1*c1/24.0;

    /* generate powers of A */
    ss->A1 = matrixCopy(ss->A);
    product(ss->A1, ss->A, ss->A2);
    product(ss->A2, ss->A, ss->A3);
    product(ss->A3, ss->A, ss->A4);

    /* scale powers of A */
    multiply(c1, ss->A1, ss->A1);
    multiply(c2, ss->A2, ss->A2);
    multiply(c3, ss->A3, ss->A3);
    multiply(c4, ss->A4, ss->A4);

    for (col = 1; col <= ss->I->cols; col++)
        for (row = 1; row <= ss->I->rows; row++)
            ELEM(ss->RK4_A, row, col) = ELEM(ss->I, row, col)
                                      + ELEM(ss->A1, row, col)
                                      + ELEM(ss->A2, row, col)
                                      + ELEM(ss->A3, row, col)
                                      + ELEM(ss->A4, row, col);
    /* scale I */
    multiply(c1, ss->I, ss->I);

    /* rescale powers of A */
    multiply(c2/c1, ss->A1, ss->A1);
    multiply(c3/c2, ss->A2, ss->A2);
    multiply(c4/c3, ss->A3, ss->A3);


    for (col = 1; col <= ss->I->cols; col++)
        for (row = 1; row <= ss->I->rows; row++)
            ELEM(ss->RK4_B, row, col) = ELEM(ss->I, row, col)
                                      + ELEM(ss->A1, row, col)
                                      + ELEM(ss->A2, row, col)
                                      + ELEM(ss->A3, row, col); 

    return;
}

