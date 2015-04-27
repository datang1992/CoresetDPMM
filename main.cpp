#include "model_SVA.h"
#include <iostream>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <eigen3/Eigen/Dense>

using namespace std;

int main () {
	int N = 100;
	double K = 1;
	int dim = 2;
	double C = 1, l = 1, S = 1, nu = 1, nu2 = 1;
	double lambda = 1;
	double SVM_learning_rate = 0.2;
	double SVM_iterations = 100;
	double coreset_epsilon = 1e-3;
	int initial_cluster_number = 10;
	int number_of_iterations = 100;
	SVA_model model(N, K, dim, C, l, S, nu, nu2, lambda, SVM_learning_rate, SVM_iterations, coreset_epsilon, initial_cluster_number, number_of_iterations);
	for (int i = 0; i < N; i++) {
		Vector2d v;
		v(0) = gsl_ran_gaussian(model.rng, 1);
		v(1) = gsl_ran_gaussian(model.rng, 1);
		model.x.push_back(v);
		if (v(0) + v(1) > 0)
			model.y.push_back(1);
		else
			model.y.push_back(-1);
	}
	model.find_bacteria_solution();
	model.compute_coreset();
	model.M2DPM();
	cout << "Accuracy: " << setiosflags(ios::fixed) << setprecision(2) << 100 * model.validation() << "%." << endl;
	return 0;
}
