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
    matrix *Ad;
    matrix *Bd;

    /* time step */
    double h;
} stateSpace;

stateSpace * stateSpaceCreate(const int size_x, const int size_u, const int size_y);
int stateSpaceDestroy(stateSpace * ss);
void stateSpaceUpdate(stateSpace * ss);
int calculateRungeKutta(stateSpace * ss);


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
    ss->Ad = NULL;
    ss->Bd = NULL;

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
    ss->Ad = matrixCreate(size_x, size_x);
    ss->Bd = matrixCreate(size_x, size_x);

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
        !ss->Ad ||
        !ss->Bd) {
        stateSpaceDestroy(ss);
        return NULL;
    }

    /* default step size */
    ss->h = 1e-3;
    
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
    if (ss->Ad)
        if (matrixDestroy(ss->Ad)) { failed = -22; }
    if (ss->Bd)
        if (matrixDestroy(ss->Bd)) { failed = -23; }

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
    product(ss->Ad, ss->x, ss->vec_x1);
    product(ss->B, ss->u, ss->vec_x2);
    product(ss->Bd, ss->vec_x2, ss->vec_x3);

    /* store old data in unused matrix */
    matrixCopyData(ss->x, ss->vec_x2);
    /* update x_dot (unneccesary) */
    difference(ss->x, ss->vec_x2, ss->x_dot);

    /* update x */
    sum(ss->vec_x1, ss->vec_x3, ss->x);

    return;
}


int calculateRungeKutta(stateSpace * ss) {
    int row, col;
    matrix *A1, *A2, *A3, *A4, *I;

    double c1 = ss->h;
    double c2 = c1*c1/2.0;
    double c3 = c1*c1*c1/6.0;
    double c4 = c1*c1*c1*c1/24.0;

    A1 = matrixCreate(ss->size_x, ss->size_x);
    A2 = matrixCreate(ss->size_x, ss->size_x);
    A3 = matrixCreate(ss->size_x, ss->size_x);
    A4 = matrixCreate(ss->size_x, ss->size_x);
    I = matrixCreate(ss->size_x, ss->size_x);
    identity(I);

    if (!A1 ||
        !A2 ||
        !A3 ||
        !A4 ||
        !I) {
        return 1;
    }

    /* generate powers of A */
    A1 = matrixCopy(ss->A);
    product(A1, ss->A, A2);
    product(A2, ss->A, A3);
    product(A3, ss->A, A4);

    /* scale powers of A */
    multiply(c1, A1, A1);
    multiply(c2, A2, A2);
    multiply(c3, A3, A3);
    multiply(c4, A4, A4);

    for (col = 1; col <= I->cols; col++)
        for (row = 1; row <= I->rows; row++)
            ELEM(ss->Ad, row, col) = ELEM(I, row, col)
                                      + ELEM(A1, row, col)
                                      + ELEM(A2, row, col)
                                      + ELEM(A3, row, col)
                                      + ELEM(A4, row, col);
    /* scale I */
    multiply(c1, I, I);

    /* rescale powers of A */
    multiply(c2/c1, A1, A1);
    multiply(c3/c2, A2, A2);
    multiply(c4/c3, A3, A3);

    for (col = 1; col <= I->cols; col++)
        for (row = 1; row <= I->rows; row++)
            ELEM(ss->Bd, row, col) = ELEM(I, row, col)
                                      + ELEM(A1, row, col)
                                      + ELEM(A2, row, col)
                                      + ELEM(A3, row, col); 

    matrixDestroy(A1);
    matrixDestroy(A2);
    matrixDestroy(A3);
    matrixDestroy(A4);
    matrixDestroy(I);

    return 0;
}

