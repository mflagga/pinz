#include "dek.hpp"

void initHdlaKx(cmp *H, double kx, int Ny, double alpha, double dx){ // zakladajac ze jest calloc
    for (int j=0;j<Ny;j++){
        H[j*Ny+j]=4.0*alpha-alpha*2.0*std::cos(kx*dx);
        if (j!=Ny-1){
            H[(j+1)*Ny+j]=-alpha;
            H[j*Ny+j+1]=-alpha;
        }
    }
}

void updateHdlaKx(cmp *H, double kx, int Ny, double alpha, double dx){ // tylko przekątna
    for (int j=0;j<Ny;j++){
        H[j*Ny+j]=4.0*alpha-alpha*2.0*std::cos(kx*dx);
    }
}

void diagRealSymmEvals(cmp *H, int N, cmp *eval){
    // inicjalizacja macierzy gsl
    gsl_matrix *Hgsl = gsl_matrix_alloc(N,N);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            gsl_matrix_set(Hgsl,i,j,std::real(H[i*N+j]));
        }
    }
    // wektor własny
    gsl_vector *evalgsl = gsl_vector_alloc(N);
    // rozwiazanie
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(N);
    gsl_eigen_symm(Hgsl,evalgsl,w);
    gsl_eigen_symm_free(w);
    gsl_sort_vector(evalgsl);
    // przekaz
    for (int i=0;i<N;i++) eval[i]=gsl_vector_get(evalgsl,i);
    // czystki
    gsl_matrix_free(Hgsl);
    gsl_vector_free(evalgsl);
}

void upw_inv(double E, double alpha, int Ny, cmp *evals, cmp *evecs){
    // macierz B
    double inv_alpha = 1.0/alpha;
    Eigen::VectorXcd Binv_diag(2*Ny);
    Binv_diag.head(Ny).setOnes();
    Binv_diag.tail(Ny).setConstant(-inv_alpha);
    // maceirz H
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(Ny,Ny);
    H.diagonal().setConstant(4.0*alpha);
    H.diagonal(-1).setConstant(-alpha);
    H.diagonal(+1).setConstant(-alpha);
    // maceirz A
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(2*Ny,2*Ny);
    A.block(0,Ny,Ny,Ny)=Eigen::MatrixXcd::Identity(Ny,Ny);
    A.block(Ny,0,Ny,Ny)=alpha*Eigen::MatrixXcd::Identity(Ny,Ny);
    A.block(Ny,Ny,Ny,Ny)=E*Eigen::MatrixXcd::Identity(Ny,Ny) - H;
    // iloczyn
    Eigen::MatrixXcd BinvA = Binv_diag.asDiagonal() * A;
    // rozwiązanie
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(BinvA);
    // wektory i wartości wł
    Eigen::VectorXcd evals_eig = solver.eigenvalues();
    Eigen::MatrixXcd evecs_eig = solver.eigenvectors();
    for (int i=0;i<Ny;i++){
        evals[i] = evals_eig(i);
        for (int j=0;j<Ny;j++){
            evecs[i*Ny+j] = evecs_eig(j,i);
        }
    }
}

// void MrazyM(cmp *A, cmp *B, cmp *C, int n, int m, int p){
//     for (int i=0;i<n;i++){
//         for (int j=0;j<p;j++){
//             C[i*p+j]=0.0;
//             for (int k=0;k<m;k++){
//                 C[i*p+j] += A[i*m+k]*B[k*p+j];
//             }
//         }
//     }
// }
