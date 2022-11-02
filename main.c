#include <stdio.h>
#include <math.h>

#define left 0
#define right 1

int main(int argc, char *argv[]) {

}

void RungeKutta(int n, double alph, double bett, double *u[]) { // this func count value of u in spec points
	double C[4][4];
	double h = 1. / n;
	
	u[0][0] = 0;
	u[1][0] = 0;
	u[2][0] = alph;
	u[3][0] = bett;

	for (int i = 0; i < n - 1; i ++) {
		C[0][0] = u[i][1];
		C[0][1] = u[i][2];
		C[0][2] = u[i][3];
		C[0][3] = /*g(x)*/ - u[i][3] - u[i - 1][0] * sin(i*h);

		C[1][0] = u[i][1] + C[0][1] * h / 2.;
		C[1][1] = u[i][2] + C[0][2] * h / 2.;
		C[1][2] = u[i][3] + C[0][3] * h / 2.;
		C[1][3] = /*g(x)*/ - (u[i][3] + C[0][3] * h / 2.) - (u[i][0] + C[0][0] * h / 2.) * sin((i + 1./2.)*h);

		C[2][0] = u[i][1] + C[1][1] * h / 2.;
		C[2][1] = u[i][2] + C[1][2] * h / 2.;
		C[2][2] = u[i][3] + C[1][3] * h / 2.;
		C[2][3] = /*g(x)*/ - (u[i][3] + C[1][3] * h / 2.) - (u[i][0] + C[1][0] * h / 2.) * sin((i + 1./2.)*h);

		C[3][0] = u[i][1] + C[0][1] * h / 2.;
		C[3][1] = u[i][2] + C[0][2] * h / 2.;
		C[3][2] = u[i][3] + C[0][3] * h / 2.;
		C[3][3] = /*g(x)*/ - (u[i][3] + C[0][3] * h) - (u[i][0] + C[0][0] * h) * sin((i + 1)*h);

		for (int j = 0; j < 4; j ++) {
			u[i + 1][j] = u[i][j] + (h / 2.)*(C[0][j] + 2*C[1][j] + 2*C[2][j] + C[3][j]);
		}
	}
}

double RKRes(int n, int k, double x, double *u[]) { // find value of k-th func in point x
    int rb = 1;
	double h  = 1. / n;
    for (int i = 1; i < n; i++) {
        if (x < i*h) {
            rb = i;
            break;
        }
    }
    
	return (x - h*(rb - 1)) / h * (u[k][rb] - u[k][rb - 1]);
}

void findEdge(int n, double *edge) { // func ready, it count edge parametrs
	
	int err = 0;
	double A_0[3][3], b_0[3], A_1[3][3], b_1[3];
	double abc_0[3], abc_1[3];
	double A[2][2], b[2];
	double u[4][n];

	for (int i = 0; i < 3; i ++) {
		
		A_0[i][0] = i / 2.;
		A_0[i][1] = 1 - i / 2.;
		A_0[i][2] = 1;
        RungeKutta(n, A_0[i][0], A_0[i][1], u);
		b_0[i] = u[0][n - 1];

		A_1[i][0] = i / 2.;
		A_1[i][1] = 1 - i / 2.;
		A_1[i][2] = 1;
        RungeKutta(n, A_1[i][0], A_1[i][1], u);
		b_1[i] = u[1][n - 1];
	}

	solveLineralSystem(3, A_0, b_0, abc_0, &err);
	solveLineralSystem(3, A_1, b_1, abc_1, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}

	A[0][0] = abc_0[0];
	A[0][1] = abc_0[1];
	b[0] = (-1)*abc_0[2];

	A[1][0] = abc_1[0];
	A[1][1] = abc_1[1];
	b[1] = (-1)*abc_1[2];

	solveLineralSystem(2, A, b, edge, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}
}