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
    // zapisanie 2Ny wartosci i wektorow (o polowe krotszych)
    for (int i=0;i<2*Ny;i++){
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

void znajdzPsi(cmp *psi, double alpha, double c_in, cmp *evecs_pop, int l, int Nx, int Ny, cmp *evals_pop, int liczba, double E, double *V){
    cmp lambdanminus = evals_pop[l];
    cmp lambdanplus = 1.0/lambdanminus;
    cmp deltaplus = 1.0-1.0/lambdanplus;
    cmp deltaminus = 1.0-1.0/lambdanminus;
    // wektor wyrazow wolnych
    Eigen::VectorXcd www_eig = Eigen::VectorXcd::Zero(Nx*Ny);
    for (int j=0;j<Ny;j++){
        www_eig(j) = -alpha*c_in*evecs_pop[l*Ny+j]*(deltaplus-deltaminus);
    }
    // macierz
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(Nx*Ny,Nx*Ny);
    // hamiltonian
    Eigen::MatrixXcd H_eig = Eigen::MatrixXcd::Zero(Ny,Ny);
    H_eig.diagonal(-1).setConstant(-alpha);
    H_eig.diagonal(+1).setConstant(-alpha);
    // macierz tau
    Eigen::MatrixXcd tau = Eigen::MatrixXcd::Zero(Ny,Ny);
    tau.diagonal().setConstant(-alpha);
    // macierze alfa i beta
    Eigen::MatrixXcd alphaM(Ny,Ny);
    Eigen::MatrixXcd betaM(Ny,Ny);
    for (int mu=0;mu<Ny;mu++){
        for (int ni=0;ni<Ny;ni++){
            alphaM(mu,ni)=0.0;
            betaM(mu,ni)=0.0;
            for (int n=0;n<liczba;n++){
                alphaM(mu,ni) += std::conj(evecs_pop[n*Ny+ni])*evecs_pop[n*Ny+mu]*(1.0-1.0/evals_pop[n]);
                betaM(mu,ni) += std::conj(evecs_pop[n*Ny+ni])*evecs_pop[n*Ny+mu]*(1.0-1.0/evals_pop[n]);
            }
        }
    }
    // wypelnianie macierzy
    for (int i=0;i<Nx;i++){
        // wypelnienie diagonali hamiltonianu
        for (int j=0;j<Ny;j++) H_eig(j,j) = 4.0*alpha + V[i*Ny+j];
        // wypelnienie macierzy
        if (i==0){
            M.block(i*Ny,i*Ny,Ny,Ny) = H_eig - E*Eigen::MatrixXcd::Identity(Ny,Ny) + tau - tau*alphaM;
        }
        else if (i==Nx-1){
            M.block(i*Ny,i*Ny,Ny,Ny) = H_eig - E*Eigen::MatrixXcd::Identity(Ny,Ny) + tau.adjoint() - tau.adjoint()*betaM;
        }
        else{
            M.block(i*Ny,i*Ny,Ny,Ny) = H_eig - E*Eigen::MatrixXcd::Identity(Ny,Ny);
        }
        if (i!=Nx-1){
            M.block((i+1)*Ny,i*Ny,Ny,Ny) = tau;
            M.block(i*Ny,(i+1)*Ny,Ny,Ny) = tau.adjoint();
        }
    }
    // rozwiazanie
    Eigen::VectorXcd psi_eig = M.lu().solve(www_eig);
    // zapisanie
    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            psi[i*Ny+j] = psi_eig(i*Ny+j);
        }
    }
}

int deltaK(int n, int m){
    return n==m?1:0;
}

void initV_QPC(double *V, int Nx, int Ny, double *x, double *y, double ymin, double ymax, double nm_to_bohr){
    const double sigmax=300.0*nm_to_bohr;
    const double sigmay=300.0*nm_to_bohr;
    const double Vg=-0.05;
    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            V[i*Ny+j]=-0.035*Vg*(std::exp(-std::pow(std::pow(x[i]/sigmax,2)+std::pow((y[j]-ymin)/sigmay,2),2))+std::exp(-std::pow(std::pow(x[i]/sigmax,2)+std::pow((y[j]-ymax)/sigmay,2),2)));
        }
    }
}
