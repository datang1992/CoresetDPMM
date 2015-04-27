

#include "model_SVA.h"
#include "utilities.h"
#include <time.h>
#include <map>
#include <iostream>

using namespace std;
using namespace Eigen;

SVA_model::SVA_model(int K, int N, int dim, double C, double S, double nu, double alpha, double lambda, double SVM_learning_rate, int SVM_iterations, double coreset_epsilon): K(K), N(N), dim(dim), C(C), S(S), nu(nu), alpha(alpha), lambda(lambda), SVM_learning_rate(SVM_learning_rate), SVM_iterations(SVM_iterations), coreset_epsilon(coreset_epsilon) {
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
	// SVM, using subgradient optimization
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
			bac_sol_classifier[i] -= SVM_learning_rate * temp_classifier[i]; // Subgradient Update
		for (int i = 0; i < N; i++)
			if (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
				bac_sol_classifier[i] += 2 * nu * nu * C * SVM_learning_rate * y[i] * x[i];	// Subgradient Update			
	}
}

void SVA_model::compute_coreset() {
	double alpha = 16 * (log(bac_sol.size()) / log(2) + 2) + 10;
	double *s = new double[N];
	double *dis_square_sum = new double[bac_sol.size()];
	int *num_point = new int[bac_sol.size()];
	memset(dis_square_sum, 0, sizeof(double) * bac_sol.size());
	memset(num_point, 0, sizeof(int) * bac_sol.size());
	for (int i = 0; (unsigned int) i < bac_sol.size(); i++) {
		num_point[bac_sol_assign[i]]++;
		dis_square_sum[bac_sol_assign[i]] += 2 * S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[i] - x[bac_sol[bac_sol_assign[i]]]);
		if (y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 1)
			dis_square_sum[bac_sol_assign[i]] += C * (y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) - 1);
	}
	double Loss = lambda * K;
	for (int i = 0; (unsigned int) i < bac_sol.size(); i++)
		Loss += bac_sol_classifier[i].dot(bac_sol_classifier[i]) / (2 * nu * nu);
	for (int i = 0; i < N; i++)
		if (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
			Loss += 2 * C * (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]));
	for (int i = 0; i < N; i++)
		Loss += S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[bac_sol[bac_sol_assign[i]]]);
	for (int i = 0; i < N; i++) {
		s[i] = 1 + 2 * alpha * N * S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[i] - x[bac_sol[bac_sol_assign[i]]]) / Loss;
		if (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
			s[i] += 2 * alpha * N * C * (1 - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i])) / Loss;
		s[i] += 2 * alpha * N * dis_square_sum[bac_sol_assign[i]] / (num_point[bac_sol_assign[i]] * Loss);
		s[i] += 4 * N / num_point[bac_sol_assign[i]];
	}
	int coreset_point_num = 10 * dim * bac_sol.size() * bac_sol.size() * bac_sol.size() * log(bac_sol.size()) / (coreset_epsilon * coreset_epsilon);
	
	double sum_s = 0;
	for (int i = 0; i < N; i++)
		sum_s += s[i];
	for (int i = 0; i < N; i++)
		s[i] /= sum_s;

	map<int, double> point_weight;
	gsl_ran_discrete_t *grd;
	grd = gsl_ran_discrete_preproc(N, s);

	for (int i = 0; i < coreset_point_num; i++) {
		int sample = gsl_ran_discrete(rng, grd);
		if (point_weight.find(sample) != point_weight.end())
			point_weight[sample] += 1 / (coreset_point_num * s[sample]);
		else
			point_weight[sample] = 1 / (coreset_point_num * s[sample]);	
	}

	for (map<int, double>::iterator iter = point_weight.begin(); iter != point_weight.end(); iter++) {
		coreset.push_back(iter -> first);
		coreset_weight.push_back(iter -> second);
	}

	gsl_ran_discrete_free(grd);
	delete[] dis_square_sum;
	delete[] num_point;
	delete[] s;
}

void SVA_model::M2DPM() {

}

void SVA_model::map_back() {

}

void SVA_model::cross_validiction() {

}


