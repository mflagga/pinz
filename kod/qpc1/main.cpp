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
    const double E=0.4*eV_to_hartree;

    // parametry symulacji
    const int Nx=1*10-1;
    const int Ny=1*20-1;
    
    // parametry wtórne
    const double dx=(ymax-ymin)/(Ny+1);
    // const double dy=(ymax-ymin)/(Ny-1);
    const double alpha=1.0/(2.0*m*dx*dx);

    // alokacja / inicjalizacja
    double *x=new double[Nx];
    for (int i=0;i<Nx;i++) x[i]=xmin+(i+1)*dx;
    double *y=new double[Ny];
    for (int j=0;j<Ny;j++) y[j]=ymin+(j+1)*dx;
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
        fprintf(kxfile,"%e",kx*1.0/bohr_to_nm);
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
            fprintf(dyspfile,"%e",std::real(dysp[p][j])*hartree_to_eV);
        }
        fprintf(dyspfile,"\n");
    }
    fclose(dyspfile);

    // uogólniony problem własny

    cmp *evals = new cmp[Ny];
    cmp *evecs = new cmp[Ny*Ny];
    upw_inv(E,alpha,Ny,evals,evecs);

    // normowanie funckji wlasnych
    double sum, invsqrtsum;
    for (int i=0;i<Ny;i++){
        sum=0.0;
        for (int j=0;j<Ny;j++) sum += std::norm(evecs[i*Ny+j])*dx;
        invsqrtsum=1.0/std::sqrt(sum);
        for (int j=0;j<Ny;j++) evecs[i*Ny+j]*=invsqrtsum;
    }

    // znalezienie modów poprzecznych
    int liczba=0;
    bool *poprzeczne = new bool[Ny]{};
    double eps=1e-8;
    for (int i=0;i<Ny;i++){
        if (std::abs(evals[i])>=1.0-eps && std::abs(evals[i])<=1.0+eps){
            poprzeczne[i]=true;
            liczba++;
        }
    }
    printf("Liczba modów poprzecznych: %d\n",liczba);

    // zapisanie
    FILE *lfile=fopen("lfile.csv","w");
    FILE *ufile=fopen("ufile.csv","w");
    FILE *poprzecznefile=fopen("poprzecznefile.csv","w");
    for (int i=0;i<Ny;i++){
        if (i!=0) fprintf(lfile,",");
        fprintf(lfile,"%e",std::abs(evals[i]));
        if (i!=0) fprintf(poprzecznefile,",");
        fprintf(poprzecznefile,"%d",poprzeczne[i]);
        for (int j=0;j<Ny;j++){
            fprintf(ufile,"%e,%e\n",std::real(evecs[i*Ny+j])*1.0/std::sqrt(bohr_to_nm),std::imag(evecs[i*Ny+j])*1.0/std::sqrt(bohr_to_nm));
        }
    }
    fclose(lfile);
    fclose(ufile);
    fclose(poprzecznefile);

    cmp *evals_pop = new cmp[liczba];
    cmp *evecs_pop = new cmp[liczba*Ny];

    int l=0;
    for (int i=0;i<Ny;i++){
        if (poprzeczne[i]){
            evals_pop[l]=evals[i];
            for (int j=0;j<Ny;j++){
                evecs_pop[l*Ny+j] = evecs[i*Ny+j];
            }
            l++;
        }
    }

    // prędkości i pędy
    cmp is;
    double *v = new double[liczba];
    double *kx2 = new double[liczba];
    for (int i=0;i<liczba;i++){
        // pęd
        kx2[i]=std::imag(std::log(evals_pop[i])/dx);
        // iloczyn skalarny
        is=0.0+0.0*I;
        for (int j=0;j<Ny;j++){
            is += evecs_pop[i*Ny+j]*std::conj(evecs_pop[i*Ny+j]);
        }
        // predkosc
        v[i]=-2.0*dx*(-alpha)*std::imag(evals_pop[i]*is);
        // printf("%lf\n",v[i]);
    }

    // zapisanie 
    FILE *kx2vfile=fopen("kx2vfile.csv","w");
    for (int i=0;i<liczba;i++) fprintf(kx2vfile,"%e,%e\n",kx2[i]*nm_to_bohr,v[i]);
    fclose(kx2vfile);

    // wspolczynnik transmisji

    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%lf,%d,%d,%e,%e",E*hartree_to_eV,Ny,Nx,xmin*bohr_to_nm,ymin*bohr_to_nm);
    fclose(misc);
    

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
    delete [] v;
    delete [] poprzeczne;
    delete [] evals_pop;
    delete [] evecs_pop;
    delete [] kx2;

    // return zero
    return 0;
}
