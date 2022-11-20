#include <stdio.h>
#include <malloc.h>
#include <math.h>

#define syst 10
#define EPS 1e-2

// Gauss func
int chooseMain(double *a, int n, int num);
void swapLines(double *a, double *b, int n, int m_n, int num);
void solveLineralSystem(int n, double *a, double *b, double *x, int *err);
void printer(int n, double *a, double *b);

// Runge-Kutta func
double g(double x);
double sol(double x);
void RungeKutta(int n, double alph, double bett, double **u);
void findEdge(int n, double *edge);