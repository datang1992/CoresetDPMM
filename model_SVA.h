#ifndef __MODEL_SVA_H__
#define __MODEL_SVA_H__


#include <vector>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

using namespace std;


class SVA_model {
	public:
		int K,  N;
		int dim;
		vector<double> y;
		vector<gsl_matirx*> x;
		vector<int> z;
		vectore<double> zeta;
		vector<int> Num;
		vector<gsl_matrix*> eta;
		
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
