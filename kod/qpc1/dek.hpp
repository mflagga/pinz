#ifndef DEK_HPP_
#define DEK_HPP_

#include <iostream>
#include <cmath>
#include <complex>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

typedef std::complex<double> cmp;
const cmp I(0.0,1.0);

void initHdlaKx(cmp *H, double kx, int Ny, double alpha, double dx);
void updateHdlaKx(cmp *H, double kx, int Ny, double alpha, double dx);
void diagRealSymmEvals(cmp *H, int N, cmp *eval);
void upw_inv(double E, double alpha, int Ny, cmp *evals, cmp *evecs);
// void MrazyM(cmp *A, cmp *B, cmp *C, int n, int m, int p);
void znajdzPsi(cmp *psi, double alpha, double c_in, cmp *evecs_pop, int l, int Nx, int Ny, cmp lambdanminus, int liczba, double E, double *V);

#endif
