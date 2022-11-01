#include <stdio.h>

int main(int argc, char *argv[]) {

}

void RungeKutta(double alph, double bett, double *C) {
	// here will be formulas, which count value R-K in points
}

double RKRes(int n, int k, double x, double a, double b, double *C) {
    
}

void findEdge(double *edge) { // func ready, it count edge parametrs
	
	int err = 0;
	double A_0[3][3], b_0[3], A_1[3][3], b_1[3];
	double abc_0[3], abc_1[3];
	double A[2][2], b[2];

	for (int i = 0; i < 3; i ++) {
		
		A_0[i][0] = i / 2.;
		A_0[i][1] = 1 - i / 2.;
		A_0[i][2] = 1;
        RungeKutta(A_0[i][0], A_0[i][1]);
		b_0[i] = RKRes(/*arguments*/);

		A_1[i][0] = i / 2.;
		A_1[i][1] = 1 - i / 2.;
		A_1[i][2] = 1;
        RungeKutta(A_0[i][0], A_0[i][1]);
		b_1[i] = RKRes(/*arguments*/);
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