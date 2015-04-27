

#include "model_SVA.h"
#include "utilities.h"
#include <time.h>
#include <iostream>

using namespace std;
using namespace Eigen;

SVA_model::SVA_model(int K, int N, int dim, double C, double l, double nu, double alpha, double lambda, double SVM_learning_rate, int SVM_iterations): K(K), N(N), dim(dim), C(C), l(l), nu(nu), alpha(alpha), lambda(lambda), SVM_learning_rate(SVM_learning_rate), SVM_iterations(SVM_iterations) {
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
	double *distance = new double[N];
	for (int i = 0; i < N; i++)
		distance[i] = 1.0;
	gsl_ran_discrete_t *grd;
	grd = gsl_ran_discrete_preproc(N, distance);
	int num_first = gsl_ran_discrete(rng, grd);
	bool *mark = new bool[N];
	for (int i = 0; i < N; i++)
		mark[i] = false;
	mark[num_first] = true;
	bac_sol.push_back(num_first);
	vector<VectorXd> bac_sol_point;
	bac_sol_point.push_back(x[num_first]);
	gsl_ran_discrete_free(grd);
	double dis = 0;
	int *assign = new int[N];
	cal_dp_means_dis(x, bac_sol_point, dis, assign, distance);
	while (bac_sol.size() < (unsigned int)N && dis > 16 * lambda * bac_sol.size() * (log(bac_sol.size()) / log(2) + 2)) {
		for (int i = 0; i < N; i++)
			distance[i] += epsilon;
		grd = gsl_ran_discrete_preproc(N, distance);
		do {
			num_first = gsl_ran_discrete(rng, grd);
		}while (!mark[num_first]);
		mark[num_first] = true;
		bac_sol.push_back(num_first);
		bac_sol_point.push_back(x[num_first]);
		cal_dp_means_dis(x, bac_sol_point, dis, assign, distance);
	}
	for (int i = 0; i < N; i++)
		bac_sol_assign.push_back(assign[i]);
	delete[] assign;
	delete[] distance;
	delete[] mark;
	bac_sol_classifier.clear();
	for (int i = 0; (unsigned int)i < bac_sol.size(); i++) {
		VectorXd temp = MatrixXd::Zero(dim, 1);
		bac_sol_classifier.push_back(temp);
	}
	for (int iter = 0; iter < SVM_iterations; iter++) {
		vector<VectorXd> temp_classifier;
		for (int i = 0; (unsigned int) i < bac_sol.size(); i++)
			temp_classifier.push_back(bac_sol_classifier[i]);
		for (int i = 0; (unsigned int) i < bac_sol.size(); i++)
			bac_sol_classifier[i] -= SVM_learning_rate * temp_classifier[i];
		for (int i = 0; i < N; i++)
			if (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
				bac_sol_classifier[i] += SVM_learning_rate * y[i] * x[i];					
	}
}

void SVA_model::compute_coreset() {
	
}

void SVA_model::M2DPM() {

}

void SVA_model::map_back() {

}

void SVA_model::cross_validiction() {

}


