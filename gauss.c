#include "mylib.h"

int chooseMain(double *a, int n, int num)
{
    double v_max = a[num * n + num];
    int n_max = num;
    for (int i = num + 1; i < n; i++) {
        if (fabs(a[i * n + num]) - fabs(v_max) > EPS) {
            n_max = i;
            v_max = a[i * n + num];
        }
    }
    return n_max;
}

void swapLines(double *a, double *b, int n, int m_n, int num)
{
    double buf_a;
    double buf_b;
    // put line into right place
    for (int i = 0; i < n; i++) {
        buf_a = a[num * n + i];
        a[num * n + i] = a[m_n * n + i];
        a[m_n * n + i] = buf_a;
    }
    buf_b = b[num];
    b[num] = b[m_n];
    b[m_n] = buf_b;
    // subsraction line from another
    for (int i = num + 1; i < n; i++) {
        buf_a = a[i * n + num];
        a[i * n + num] = 0;
        for (int j = num + 1; j < n; j++) {
            a[i * n + j] =
                a[i * n + j] - a[num * n + j] / a[num * n + num] * buf_a;
        }
        b[i] = b[i] - b[num] / a[num * n + num] * buf_a;
    }
}

void solveLineralSystem(int n, double *a, double *b, double *x, int *err)
{
    printer(n, a, b);
    for (int i = 0; i < n; i++) {
        int col_max;
        col_max = chooseMain(a, n, i);
        if ((a[i * n + col_max] > EPS) || (a[i * n + col_max] < -EPS)) {
            swapLines(a, b, n, col_max, i);
        }
        printer(n, a, b);
    }
    for (int i = 0; i < n; i++) {
        if ((a[i * n + i] > -EPS) && (a[i * n + i] < EPS)) {
            *err = -1;
            return;
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= a[i * n + j] * x[j];
        }
        x[i] /= a[i * n + i];
    }
}

void printer(int n, double *a, double *b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.3e", a[i * n + j]);
        }
        printf("%10.3e", b[i]);
        printf("\n");
    }
    printf("\n");
}