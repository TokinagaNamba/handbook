#include <stdio.h>
#include <math.h>

#define N 100

double f(double t) {
	return 1;
}
int main() {
	int i, j;
	// \partial_t^\alpha u = 1 in (0, 2), f(0) = 0;
	// u(t) = t^\alpha/\Gamma(1 + \alpha)
	// Lin, Y., & Xu, C. (2007). Finite difference/spectral approximations for the time-fractional diffusion equation. J. Comput. Phys., 225(2), 1533–1552. article. http://doi.org/10.1016/j.jcp.2007.02.001
	// double a = tgamma(alpha);
	double alpha = 0.5;
	double T = 2.0;
	double h = T / N;
	double u[N + 1];
	double b[N + 1];
	for (i = 0; i <= N; i++) {
		double a = pow(i + 1, 1 - alpha) - pow(i, 1 - alpha);
		b[i] = a;
	}
	double alpha0 = tgamma(2 - alpha) * pow(h, alpha);
	// initial condition
	u[0] = 0;
	for (i = 1; i <= N; i++) {
		// k + 1 = i
		double a = 0;
		a = (1 - b[1]) * u[i - 1];
		for (j = 1; j <= i - 2; j++) {
			a += (b[j] - b[j + 1]) * u[i - 1 - j];
		}
		a += b[i - 1] * u[0];
		a += alpha0 * f(i * h);
		u[i] = a;
	}
	for (i = 0; i <= N; i++) {
		// pow(u[i] * tgamma(1 + alpha), 1 / alpha)
		printf("%f %f\n", i * h, u[i]);
	}
	return 0;
}
