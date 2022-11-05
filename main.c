#include "mylib.h"

// #define left 0
// #define right 1
#define syst 10

int main() {
	FILE *out;
	double e[2], **u;
	int n = 20;
	char *endptr;
	double h;

	printf("nana\n");

	// scanf("%d", n);
	// n = (int)strtol(argv[1], &endptr, syst);
    // if (n <= 0) {
    //     printf("Incorrect value for argument n\n");
    //     return -1;
    // }

	h = 1. / n;

	u = (double **)malloc(sizeof(double *)*4);
	for (int i = 0; i < n; i ++) {
		u[i] = (double *)malloc(sizeof(double)*n);
	}
	
	// findEdge(n, e);
	RungeKutta(n, -pow(M_PI, 2) / 2., 0, u);

	out = fopen("res.txt", "w");

	for (int i = 0; i < n; i ++) {
		fprintf(out, "%lf %lf\n", i*h, u[0][i]);
	}

	printf("finish\n");
	return 0;
}

double g(double x) {
	return  - pow(M_PI, 4) / 2. * cos(M_PI*x) - pow(M_PI, 3) * sin(M_PI*x) / 2. + sin(x)*(1 - cos(M_PI*x)) / 2.;
}

void RungeKutta(int n, double alph, double bett, double **u) { // this func count value of u in spec points
	double C[4][4];
	double h = 1. / n;
	
	u[0][0] = 0;
	u[1][0] = 0;
	u[2][0] = alph;
	u[3][0] = bett;

	for (int i = 0; i < n - 1; i ++) {
		C[0][0] = u[1][i];
		C[0][1] = u[2][i];
		C[0][2] = u[3][i];
		C[0][3] = g(i*h) - u[3][i] - u[0][i] * sin(i*h);

		C[1][0] = u[1][i] + C[0][1] * h / 2.;
		C[1][1] = u[2][i] + C[0][2] * h / 2.;
		C[1][2] = u[3][i] + C[0][3] * h / 2.;
		C[1][3] = g((i + 1./2.)*h) - (u[3][i] + C[0][3] * h / 2.) - (u[0][i] + C[0][0] * h / 2.) * sin((i + 1./2.)*h);

		C[2][0] = u[1][i] + C[1][1] * h / 2.;
		C[2][1] = u[2][i] + C[1][2] * h / 2.;
		C[2][2] = u[3][i] + C[1][3] * h / 2.;
		C[2][3] = g((i + 1./2.)*h) - (u[3][i] + C[1][3] * h / 2.) - (u[0][i] + C[1][0] * h / 2.) * sin((i + 1./2.)*h);

		C[3][0] = u[1][i] + C[2][1] * h;
		C[3][1] = u[2][i] + C[2][2] * h;
		C[3][2] = u[3][i] + C[2][3] * h;
		C[3][3] = g((i + 1)*h) - (u[3][i] + C[2][3] * h) - (u[0][i] + C[2][0] * h) * sin((i + 1)*h);

		printf("step %d: %lf %lf %lf %lf\n", i, C[0][0], C[1][0], C[2][0], C[3][0]);

		for (int j = 0; j < 4; j ++) {
			u[j][i + 1] = u[j][i] + (h / 6.)*(C[0][j] + 2*C[1][j] + 2*C[2][j] + C[3][j]);
		}
	}
}

/*double RKRes(int n, int k, double x, double *u[]) { // find value of k-th func in point x
    int rb = 1;
	double h  = 1. / n;
    for (int i = 1; i < n; i++) {
        if (x < i*h) {
            rb = i;
            break;
        }
    }
    
	return (x - h*(rb - 1)) / h * (u[k][rb] - u[k][rb - 1]);
}*/

void findEdge(int n, double *edge) { // func ready, it count edge parametrs
	
	int err = 0;
	double A_0[9], b_0[3], A_1[9], b_1[3];
	double abc_0[3], abc_1[3];
	double A[4], b[2];
	double **u;

	u = (double **)malloc(sizeof(double *)*4);
	for (int i = 0; i < 4; i ++) {
		u[i] = (double *)malloc(sizeof(double)*n);
	}

	for (int i = 0; i < 3; i ++) {
		
		A_0[i*3] = (i + 1)*(i + 1) / 2.;
		A_0[i*3 + 1] = 1 - i / 3.;
		A_0[i*3 + 2] = 1;
        RungeKutta(n, A_0[i*3], A_0[i*3 + 1], u);
		b_0[i] = u[0][n - 1];
		printf("b_0: %lf, ", b_0[i]);

		A_1[i*3] = (i + 1)*(i + 1) / 2.;
		A_1[i*3 + 1] = 1 - i / 3.;
		A_1[i*3 + 2] = 1;
        RungeKutta(n, A_1[i*3], A_1[i*3 + 1], u);
		b_1[i] = u[1][n - 1];
		printf("b_1: %lf\n", b_1[i]);
	}

	solveLineralSystem(3, A_0, b_0, abc_0, &err);
	printf("%d\n", err);
	solveLineralSystem(3, A_1, b_1, abc_1, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}

	A[0] = abc_0[0];
	A[1] = abc_0[1];
	b[0] = (-1)*abc_0[2];

	A[2] = abc_1[0];
	A[3] = abc_1[1];
	b[1] = (-1)*abc_1[2];

	solveLineralSystem(2, A, b, edge, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}
}