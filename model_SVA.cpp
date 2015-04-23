#include "model_SVA.h"
#include "time.h"

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
	
}

void SVA_model::compute_coreset() {

}

void SVA_model::M2DPM() {

}

void SVA_model::map_back() {

}

void SVA_model::cross_validiction() {

}


