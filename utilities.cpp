#include "utilities.h"

using namespace std;
using namespace Eigen;


void cal_dp_means_dis(vector<VectorXd> &p, vector<VectorXd> &q, double &dis, int *assign, double *distance) {
	dis = 0;
	//#pragma omp parallel for
	for (unsigned int i = 0; i < p.size(); i++) {
		distance[i] = 1e100;
		for (unsigned int j = 0; j < q.size(); j++) {
			double dis_square = (p[i] - q[i]).dot(p[i] - q[i]);
			if (dis_square < distance[i]) {
				distance[i] = dis_square;
				assign[i] = j;
			}
		}
		dis += distance[i];
	}
}
