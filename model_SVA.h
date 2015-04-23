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
		vector<MatrixXd> x;
		vector<int> z;
		vector<double> zeta;
		vector<int> Num;
		vector<MatrixXd> eta;

		vector<int> thread_assignment;
		vector<int> origin_number;
		vector<SVA_model*> sub_models;

		int initial_cluster_number;
		
		int Nthreads;

		int number_of_iterations;
		
		vector<double> omega;

		const double C, l, nu;
		const double alpha;

		const double lambda; //DP Means parameter

		long rand_seed;
		gsl_rng *rng;

		vector<MatrixXd> bac_sol;
		vector<MatrixXd> coreset;


		void initialize();
		void find_bacteria_solution();
		void compute_coreset();
		void collect_coreset();
		void M2DPM();
		void map_back();
		void cross_validiction();
		
		SVA_model(int, int, int, double, double, double, double, double);
		~SVA_model();
};


#endif
