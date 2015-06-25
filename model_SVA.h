#ifndef __MODEL_SVA_H__
#define __MODEL_SVA_H__


#include <vector>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class SVA_model {
	public:
		const double K;
		const int N;
		const int dim; 
		vector<double> y;
		vector<VectorXd> x;
		vector<int> z;
		vector<VectorXd> mu;
		vector<VectorXd> eta;

		vector<int> thread_assignment;
		vector<int> origin_number;
		vector<SVA_model*> sub_models;

		vector<double> omega;

		const double C, l, S, nu, nu2;

		const double lambda; //DP Means parameter
		
		const double SVM_learning_rate; //For SVM;
		const int SVM_iterations; //For SVM;
		
		const double coreset_epsilon; // Coreset Error Bound

		int initial_cluster_number;

		int number_of_iterations;

		int Nthreads;

		long rand_seed;
		gsl_rng *rng;

		vector<double> predicted_y;
		double coreset_acc2;
		int coreset_used_clusters;

		const double omega_min;
		const int weight_type;

		static const double PI = 3.14159265358979323846264338327950288;
		static const double epsilon = 1e-6;

		vector<int> bac_sol;
		vector<int> bac_sol_assign;
		vector<VectorXd> bac_sol_classifier;
		vector<int> coreset;
		vector<double> coreset_weight;

		void initialize();
		void find_bacteria_solution();
		void compute_coreset();
		void collect_coreset();
		void M2DPM();
		void map_back();
		void compute_assignment();
		double cross_validation(int);
		double validation();
		double coreset_acc();

		SVA_model(double, int, int, double, double, double, double, double, double, double, int, double, int, int, double, int);
		~SVA_model();
};


#endif
