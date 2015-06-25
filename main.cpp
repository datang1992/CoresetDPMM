#include "model_SVA.h"
#include <iostream>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <eigen3/Eigen/Dense>

using namespace std;

int main () {
	int N = 10000;
	double K = 1;
	int dim = 2;
	double C = 1, l = 10000, S = 0.01, nu = 0.1, nu2 = 0.1;
	double lambda = 1;
	double SVM_learning_rate = 0.2;
	double SVM_iterations = 1000;
	double coreset_epsilon = 1e-2;
	double omega_min = 1e-3;
	int initial_cluster_number = 10;
	int number_of_iterations = 2000;
	int weight_type = 3;
	SVA_model model(K, N, dim, C, l, S, nu, nu2, lambda, SVM_learning_rate, SVM_iterations, coreset_epsilon, initial_cluster_number, number_of_iterations, omega_min, weight_type);
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
	cout << "Finish initialization!" << endl;
	model.find_bacteria_solution();
	cout << "Finish finding bacteria solution!" << endl;
	model.compute_coreset();
	
	cout << "Finish computing the coreset! Coreset size: " << model.coreset.size() << endl;
	model.M2DPM();
	cout << "Finish the M2DPM algorithm!" << endl;
	model.compute_assignment();
	cout << "Finish computing the assignments!" << endl;
	cout << "Coreset accuracy: " << setiosflags(ios::fixed) << setprecision(2) << 100 * model.coreset_acc2 << "%." << endl;
	cout << "Accuracy: " << setiosflags(ios::fixed) << setprecision(2) << 100 * model.validation() << "%." << endl;
	int *coreset_size = new int[model.mu.size()];
	memset(coreset_size, 0, sizeof(int) * model.mu.size());
	for (int i = 0; i < N; i++)
		coreset_size[model.z[i]]++;
	int used_clusters = 0;
	for (int i = 0; (unsigned int) i < model.mu.size(); i++)
		if (coreset_size[i] > 0) {
			cout << i << '\t' << setiosflags(ios::fixed) << setprecision(2) << 100.0 * coreset_size[i] / N << "%\t" << setiosflags(ios::fixed) << setprecision(6) << model.eta[i](0) << '\t' << model.eta[i](1) << '\t' << model.mu[i](0) << '\t' << model.mu[i](1) << endl;
			used_clusters++;
		}
	cout << "Used clusters: " << used_clusters << endl;
	//cout << "Coreset used clusters: " << model.coreset_used_clusters << endl;
	delete[] coreset_size;
	return 0;
}
