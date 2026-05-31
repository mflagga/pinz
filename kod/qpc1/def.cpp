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
