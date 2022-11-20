#include "mylib.h"

int main() {
	FILE *out, *func;
	double e[2], **u;
	int n = 10, max_rate = 0;
	char *endptr;
	double h, err_rate = 0, solve_rate = 0;

	printf("Enter the number of points: ");
	scanf("%d", &n);

	h = 1. / (n - 1);

	u = (double **)malloc(sizeof(double *)*4);
	for (int i = 0; i < n; i ++) {
		u[i] = (double *)malloc(sizeof(double)*n);
	}
	
	findEdge(n, e);
	RungeKutta(n, e[0], e[1], u);

	out = fopen("res.txt", "w");
	func = fopen("targ.txt", "w");

	max_rate = 0;

	for (int i = 0; i < n; i ++) {
		fprintf(out, "%lf %lf\n", i*h, u[0][i]);
	}

	for (int i = 0; i < 100; i ++) {
		fprintf(func, "%lf %lf\n", i*(1 / 99.), sol(i*(1 / 99.)));
	}

	fclose(out);
	fclose(func);

	for (int i = 0; i < n; i ++) {
		err_rate += h * (sol(i*h) - u[0][i])*(sol(i*h) - u[0][i]);
	}

	printf("Error rate: %le\n", sqrt(err_rate));
	printf("Check results in the file\n");
	return 0;
}

double g(double x) {
	return  - pow(M_PI, 4) / 2. * cos(M_PI*x) - pow(M_PI, 3) * sin(M_PI*x) / 2. + sin(x)*(1 - cos(M_PI*x)) / 2.;
}

double sol(double x) {
	return (1 - cos(M_PI*x)) / 2.;
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

		for (int j = 0; j < 4; j ++) {
			u[j][i + 1] = u[j][i] + (h / 6.)*(C[0][j] + 2*C[1][j] + 2*C[2][j] + C[3][j]);
		}
	}
}

void findEdge(int n, double *edge) { // it count edge parametrs
	
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

		A_1[i*3] = (i + 1)*(i + 1) / 2.;
		A_1[i*3 + 1] = 1 - i / 3.;
		A_1[i*3 + 2] = 1;
        RungeKutta(n, A_1[i*3], A_1[i*3 + 1], u);
		b_1[i] = u[1][n - 1];
	}

	solveLineralSystem(3, A_0, b_0, abc_0, &err);
	solveLineralSystem(3, A_1, b_1, abc_1, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}

	A[0] = abc_0[0];
	A[1] = abc_0[1];
	b[0] = 1 - abc_0[2];

	A[2] = abc_1[0];
	A[3] = abc_1[1];
	b[1] = (-1)*abc_1[2];

	solveLineralSystem(2, A, b, edge, &err);

	if (err == -1) {
		printf("Change parametrs pls\n");
		return;
	}
}