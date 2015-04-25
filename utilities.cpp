#include "utilities.h"

using namespace std;
using namespace Eigen;


double dp_means_dis(vector<VectorXd> &p, vector<VectorXd> &q) {
	double dis = 0;
	//#pragma omp parallel for
	for (unsigned int i = 0; i < p.size(); i++) {
		double min_dis = 1e100;
		for (unsigned int j = 0; j < q.size(); j++) {
			double dis_square = (p[i] - q[i]).dot(p[i] - q[i]);
			if (dis_square < min_dis)
				min_dis = dis_square;
		}
		dis += min_dis;
	}
	return 0;
}
