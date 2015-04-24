#include "model_SVA.h"
#include <time.h>
#include <iostream>

using namespace std;
using namespace Eigen;

SVA_model::SVA_model(int K, int N, int dim, double C, double l, double nu, double alpha, double lambda): K(K), N(N), dim(dim), C(C), l(l), nu(nu), alpha(alpha), lambda(lambda) {
	rng = gsl_rng_alloc(gsl_rng_rand48);
	long seed = clock();
	gsl_rng_set(rng, seed);
}

SVA_model::~SVA_model() {
	gsl_rng_free(rng);
}

void SVA_model::initialize() {

}

void SVA_model::find_bacteria_solution() {
	//DP Means++ & SVM
	/*
	   DP Means++
	 */
	bac_sol.clear();
	double *P = new double[N];
	for (int i = 0; i < N; i++)
		P[i] = 1.0 / N;
	gsl_ran_discrete_t *grd;
	grd = gsl_ran_discrete_preproc(N, P);
	int num_first = gsl_ran_discrete(rng, grd);
	bac_sol.push_back(num_first);
	gsl_ran_discrete_free(grd);
	delete[] P;
}

void SVA_model::compute_coreset() {

}

void SVA_model::M2DPM() {

}

void SVA_model::map_back() {

}

void SVA_model::cross_validiction() {

}


