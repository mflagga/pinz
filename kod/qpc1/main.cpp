#include "dek.hpp"

int main(){

    // przeliczniki
    double nm_to_bohr=18.8973;
    double bohr_to_nm=1.0/nm_to_bohr;
    double eV_to_hartree=0.036749;
    double hartree_to_eV=1.0/eV_to_hartree;

    // parametry układu
    const double m=0.017;
    const double xmin=-5.0*nm_to_bohr;
    // const double xmax=-xmin;
    const double ymin=-10.0*nm_to_bohr;
    const double ymax=-ymin;
    // const double sigmax=300.0;
    // const double sigmay=300.0;

    // parametry symulacji
    const int Nx=11;
    const int Ny=21;
    
    // parametry wtórne
    const double dx=(ymax-ymin)/(Ny-1);
    // const double dy=(ymax-ymin)/(Ny-1);
    const double alpha=1.0/(2.0*m*dx*dx);

    // alokacja / inicjalizacja
    double *x=new double[Nx];
    for (int i=0;i<Nx;i++) x[i]=xmin+i*dx;
    double *y=new double[Ny];
    for (int j=0;j<Ny;j++) y[j]=ymin+j*dx;
    double *V=new double[Nx*Ny]{};
    cmp **psi = new cmp*[Nx];
    for (int i=0;i<Nx;i++) psi[i] = new cmp[Ny]{};
    cmp *H=new cmp[Ny*Ny]{};

    // wyznacznanie realcji dyspersji

    // alokacja rozwiazania
    int Nk=50;
    cmp **dysp = new cmp*[Nk];
    for (int p=0;p<Nk;p++) dysp[p]=new cmp[Ny];

    // rozwiazanie + zapisanie kx
    FILE *kxfile=fopen("kxfile.csv","w");
    double dk=(2.0*M_PI/dx)/(Nk-1);
    double kx;
    for (int p=0;p<Nk;p++){
        // nowe k
        kx=-M_PI/dx + p*dk;
        // zapisanie k
        if (p!=0) fprintf(kxfile,",");
        fprintf(kxfile,"%lf",kx*1.0/bohr_to_nm);
        // zmiana hamiltonianu
        if (p==0){
            initHdlaKx(H,kx,Ny,alpha,dx);
        }
        else{
            updateHdlaKx(H,kx,Ny,alpha,dx);
        }
        // diagonalizacja
        diagRealSymmEvals(H,Ny,dysp[p]);
    }
    fclose(kxfile);

    // zapisanie energii
    FILE *dyspfile=fopen("dyspfile.csv","w");
    for (int p=0;p<Nk;p++){
        for (int j=0;j<Ny;j++){
            if (j!=0) fprintf(dyspfile,",");
            fprintf(dyspfile,"%lf",std::real(dysp[p][j])*hartree_to_eV);
        }
        fprintf(dyspfile,"\n");
    }
    fclose(dyspfile);

    // uogólniony problem własny

    cmp *evals = new cmp[2*Ny];
    cmp *evecs = new cmp[4*Ny*Ny];
    upw_inv(0.02,alpha,Ny,evals,evecs);

    FILE *lfile=fopen("lfile.csv","w");
    FILE *ufile=fopen("ufile.csv","w");
    for (int i=0;i<2*Ny;i++){
        if (i!=0) fprintf(lfile,",");
        fprintf(lfile,"%lf",std::abs(evals[i]));
        for (int j=0;j<2*Ny;j++){
            fprintf(ufile,"%lf,%lf\n",std::real(evecs[i*2*Ny+j]),std::imag(evecs[i*2*Ny+j]));
        }
    }
    fclose(lfile);
    fclose(ufile);

    // czystki 
    delete [] x;
    delete [] y;
    delete [] V;
    for (int i=0;i<Nx;i++) delete [] psi[i];
    delete [] psi;
    delete [] H;
    for (int i=0;i<Nk;i++) delete [] dysp[i];
    delete [] dysp;
    delete [] evals;
    delete [] evecs;

    // return zero
    return 0;
}
