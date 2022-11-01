int main(int argc, char *argv[]) {

}

void RungeKutta(double alph, double bett, double *C) {
	// here will be formulas
}

double RKRes(int n, int k, double x, double a, double b, double *C) {
    
}

void findEdge(double *edge) {
	
	double C_0[3], C_1[3];
	double A_0[3][3], b_0[3], A_1[3][3], b_1[3];

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
}