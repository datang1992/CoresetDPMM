#include "utilities.h"
#include <iostream>

using namespace std;
using namespace Eigen;


void cal_dp_means_dis(vector<VectorXd> &p, vector<VectorXd> &q, double &dis, int *assign, double *distance) {
	dis = 0;
	//#pragma omp parallel for
	for (int i = 0; (unsigned int) i < p.size(); i++) {
		distance[i] = 1e100;
		for (int j = 0; (unsigned int) j < q.size(); j++) {
			double dis_square = (p[i] - q[j]).dot(p[i] - q[j]);
			if (dis_square < distance[i]) {
				distance[i] = dis_square;
				assign[i] = j;
			}
		}
		dis += distance[i];
	}
}

double hinge_plus(double x) {
	if (x > 0)
		return x;
	return 0;
}
