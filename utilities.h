#ifndef __UTILITIES_H__
#define __UTLIITIES_H__


#include <eigen3/Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;


void cal_dp_means_dis(vector<VectorXd>&, vector<VectorXd>&, double&, int*, double*);

double hinge_plus(double);

#endif
