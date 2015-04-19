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
		vectot<int> origin_number;
		vector<SVA_model*> sub_models;

		int initial_cluster_number;
		
		int Nthreads;

		int number_of_iterations;
		
		vector<double> omega;

		const double C, l, nu;
		const double alpha;

		long rand_seed;
		gsl_rng *rng;


		void initialize();
		
		SVA_model();
		~SVA_model();
};


#endif
