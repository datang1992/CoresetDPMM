

#include "model_SVA.h"
#include "utilities.h"
#include <time.h>
#include <map>
#include <iostream>

using namespace std;
using namespace Eigen;

SVA_model::SVA_model(double K, int N, int dim, double C, double l, double S, double nu, double nu2, double lambda, double SVM_learning_rate, int SVM_iterations, double coreset_epsilon, int initial_cluster_number, int number_of_iterations): K(K), N(N), dim(dim), C(C), l(l), S(S), nu(nu), nu2(nu2),  lambda(lambda), SVM_learning_rate(SVM_learning_rate), SVM_iterations(SVM_iterations), coreset_epsilon(coreset_epsilon), initial_cluster_number(initial_cluster_number), number_of_iterations(number_of_iterations) {
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
	while (bac_sol.size() < (unsigned int) N && dis > 16 * lambda * bac_sol.size() * (log(bac_sol.size() + epsilon) / log(2) + 2)) {
		for (int i = 0; i < N; i++)
			distance[i] = epsilon;
		grd = gsl_ran_discrete_preproc(N, distance);
		do {
			num_first = gsl_ran_discrete(rng, grd);
		}while (mark[num_first]);
		mark[num_first] = true;
		bac_sol.push_back(num_first);
		bac_sol_point.push_back(x[num_first]);
		if (bac_sol.size() % 10 == 0)
			cout << "Having adding " << bac_sol.size() << " points into the bacteria solution set!" << endl;
		cal_dp_means_dis(x, bac_sol_point, dis, assign, distance);
	}
	for (int i = 0; i < N; i++)
		bac_sol_assign.push_back(assign[i]);
	delete[] assign;
	delete[] distance;
	delete[] mark;
	cout << "Finish DP Means++! We got " << bac_sol.size() << " centers in the bacteria solution!" << endl;
	// SVM, using subgradient optimization
	bac_sol_classifier.clear();
	for (int i = 0; (unsigned int) i < bac_sol.size(); i++) {
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
			if (l - y[i] * temp_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
				bac_sol_classifier[bac_sol_assign[i]] += 2 * nu * nu * C * SVM_learning_rate * y[i] * x[i];	// Subgradient Update			
	}
	cout << "Finish SVM!" << endl;
}

void SVA_model::compute_coreset() {
	double alpha = 16 * (log(bac_sol.size()) / log(2) + 2) + 10;
	double *s = new double[N];
	double *dis_square_sum = new double[bac_sol.size()];
	int *num_point = new int[bac_sol.size()];
	memset(dis_square_sum, 0, sizeof(double) * bac_sol.size());
	memset(num_point, 0, sizeof(int) * bac_sol.size());
	for (int i = 0; i < N; i++) {
		num_point[bac_sol_assign[i]]++;
		dis_square_sum[bac_sol_assign[i]] += 2 * S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[i] - x[bac_sol[bac_sol_assign[i]]]);
		if (y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > l)
			dis_square_sum[bac_sol_assign[i]] += C * (y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) - l);
	}
	double Loss = lambda * K;
	for (int i = 0; (unsigned int) i < bac_sol.size(); i++)
		Loss += bac_sol_classifier[i].dot(bac_sol_classifier[i]) / (2 * nu * nu);
	for (int i = 0; i < N; i++)
		if (l - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
			Loss += 2 * C * (l - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]));
	for (int i = 0; i < N; i++)
		Loss += S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[i] - x[bac_sol[bac_sol_assign[i]]]);
	for (int i = 0; i < N; i++) {
		s[i] = 1 + 2 * alpha * N * S * (x[i] - x[bac_sol[bac_sol_assign[i]]]).dot(x[i] - x[bac_sol[bac_sol_assign[i]]]) / Loss;
		if (l - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i]) > 0)
			s[i] += 2 * alpha * N * C * (l - y[i] * bac_sol_classifier[bac_sol_assign[i]].dot(x[i])) / Loss;
		s[i] += 2 * alpha * N * dis_square_sum[bac_sol_assign[i]] / (num_point[bac_sol_assign[i]] * Loss);
		s[i] += 4 * N / num_point[bac_sol_assign[i]];
	}
	int coreset_point_num = min(2 * dim * bac_sol.size() * bac_sol.size() * bac_sol.size() * log(bac_sol.size()) / (coreset_epsilon * coreset_epsilon), N / 10.0);
	
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

void SVA_model::collect_coreset() {

}

void SVA_model::M2DPM() {
	// Initilization
	int n = coreset.size();
	double *prob = new double[initial_cluster_number];
	for (int i = 0; i < initial_cluster_number; i++)
		prob[i] = 1;
	gsl_ran_discrete_t *grd = gsl_ran_discrete_preproc(initial_cluster_number, prob);
	int *assign = new int[n];
	double *omega = new double[n];
	for	(int i = 0; i < n; i++) {
		assign[i] = gsl_ran_discrete(rng, grd);
		omega[i] = gsl_ran_gaussian(rng, 1 / nu);
	}
	eta.clear();
	mu.clear();
	for (int i = 0; i < initial_cluster_number; i++) {
		VectorXd temp(dim);
		for (int j = 0; j < dim; j++)
			temp(j) = gsl_ran_gaussian(rng, 1 / nu);
		eta.push_back(temp);
		for (int j = 0; j < dim; j++)
			temp(j) = gsl_ran_gaussian(rng, 1 / nu2);
		mu.push_back(temp);
	}
	
	vector<VectorXd> weighted_sum_x;
	vector<double> weight_sum;
	for (int i = 0; (unsigned int) i < mu.size(); i++) {
		weighted_sum_x.push_back(MatrixXd::Zero(dim, 1));
		weight_sum.push_back(0);
	}

	// Sampling Process
	for (int iter = 0; iter < number_of_iterations; iter++) {
		
		// Update z
		for (int i = 0; i < n; i++) {
			double *Q = new double[eta.size() + 1];
			for (int j = 0; (unsigned int) j < mu.size(); j++) {
				Q[j] = S * 0.5 * (x[coreset[i]] - mu[j]).dot(x[coreset[i]] - mu[j]) / (nu2 * nu2);
				if (l - y[coreset[i]] * eta[j].dot(x[coreset[i]]) > 0)
					Q[j] += 2 * C * (l - y[coreset[i]] * eta[j].dot(x[coreset[i]]));
			}
			VectorXd eta_star;
			if (2 * C * nu * nu >= 1 / (x[coreset[i]].dot(x[coreset[i]])))
				eta_star = 1 / (x[coreset[i]].dot(x[coreset[i]])) * y[coreset[i]] * x[coreset[i]];
			else
				eta_star = 2 * C * nu * nu * y[coreset[i]] * x[coreset[i]];
			Q[eta.size()] = lambda + eta_star.dot(eta_star) / (2 * nu * nu);
			if (l - y[coreset[i]] * eta_star.dot(x[coreset[i]]) > 0)
				Q[eta.size()] += 2 * C * (l - y[coreset[i]] * eta_star.dot(x[coreset[i]]));
			int best_n;
			double best_v = 1e100;
			for (int j = 0; (unsigned int) j <= mu.size(); j++)
				if (Q[j] < best_v) {
					best_v = Q[j];
					best_n = j;
				}
			if ((unsigned int) best_n < eta.size())
				assign[i] = best_n;
			else {
				assign[i] = mu.size();
				eta.push_back(eta_star);
				mu.push_back(x[coreset[i]]);
				weighted_sum_x.push_back(MatrixXd::Zero(dim, 1));
				weight_sum.push_back(0);
			}
		}
		
		// Update mu
		for (int i = 0; (unsigned int) i < mu.size(); i++) {
			weighted_sum_x[i].setZero(dim, 1);
			weight_sum[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			weighted_sum_x[assign[i]] += coreset_weight[i] * x[coreset[i]];
			weight_sum[assign[i]] += coreset_weight[i];
		}
		for (int i = 0; (unsigned int) i < eta.size(); i++) {
			if (fabs(weight_sum[i]) > epsilon)
				mu[i] = weighted_sum_x[i] / weight_sum[i];
			else {
				for (int j = 0; j < dim; j++)
					mu[i](j) = gsl_ran_gaussian(rng, 1 / nu2);
			}
		}	
		
		// Update omega
		for (int i = 0; i < n; i++)
			omega[i] = C * hinge_plus(l - y[coreset[i]] * eta[assign[i]].dot(x[coreset[i]]));
		
		// Update eta
		MatrixXd *Lambda = new MatrixXd[eta.size()];
		for (int i = 0; (unsigned int) i < mu.size(); i++)
			Lambda[i] = MatrixXd::Identity(dim, dim)  / (nu * nu);
		for (int i = 0; i < n; i++)
			Lambda[assign[i]] += C * C * x[coreset[i]] * x[coreset[i]].transpose() / omega[i];
		MatrixXd *Lambda2 = new MatrixXd[mu.size()];
		for (int i = 0; (unsigned int) i < mu.size(); i++)
			Lambda2[i].setZero(dim, 1);
		for (int i = 0; i < n; i++)
			Lambda2[assign[i]] += C * y[coreset[i]] * (omega[i] + C * l) / omega[i] * x[coreset[i]];
		for (int i = 0; (unsigned int) i < eta.size(); i++)
			eta[i] = Lambda[i].inverse() * Lambda2[i];
		delete[] Lambda;
		delete[] Lambda2;
	}
	
	gsl_ran_discrete_free(grd);
	delete[] prob;
	delete[] assign;
	delete[] omega;
}

void SVA_model::map_back() {

}

void SVA_model::compute_assignment() {
	z.clear();
	predicted_y.clear();
	for (int i = 0; i < N; i++) {
		double min_loss = 1e100;
		z.push_back(0);
		for (int j = 0; (unsigned int) j < mu.size(); j++)
			if (2 * C * hinge_plus(l - y[i] * eta[j].dot(x[i])) + S * (x[i] - mu[j]).dot(x[i] - mu[j]) < min_loss) {
				min_loss = 2 * C * hinge_plus(l - y[i] * eta[j].dot(x[i])) + S * (x[i] - mu[j]).dot(x[i] - mu[j]);
				z[i] = j;
			}
		if (eta[z[i]].dot(x[i]) > 0)
			predicted_y.push_back(1);
		else
			predicted_y.push_back(-1);
	}
}

double SVA_model::cross_validation(int num_fold) {
	return 0;
}

double SVA_model::validation() {
	double acc = 0;
	for (int i = 0; i < N; i++)
		if (y[i] * predicted_y[i] > 0)
			acc += 1.0 / N;
	return acc;
}
