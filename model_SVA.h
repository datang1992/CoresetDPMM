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
		const int K,  N;
		const int dim; 
		vector<double> y;
		vector<VectorXd> x;
		vector<int> z;
		vector<double> zeta;
		vector<int> Num;
		vector<VectorXd> eta;

		vector<int> thread_assignment;
		vector<int> origin_number;
		vector<SVA_model*> sub_models;

		int initial_cluster_number;
		
		int Nthreads;

		int number_of_iterations;
		
		vector<double> omega;

		const double C, S, nu;
		const double alpha;

		const double lambda; //DP Means parameter
		
		const double SVM_learning_rate; //For SVM;
		const int SVM_iterations; //For SVM;
		
		const double coreset_epsilon; // Coreset Error Bound

		long rand_seed;
		gsl_rng *rng;

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
		void cross_validiction();
		
		SVA_model(int, int, int, double, double, double, double, double, double, int, double);
		~SVA_model();
};


#endif
