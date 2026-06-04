#include "dek.hpp"

int main(){

    // przeliczniki
    double nm_to_bohr=18.8973;
    double bohr_to_nm=1.0/nm_to_bohr;
    double eV_to_hartree=0.036749;
    double hartree_to_eV=1.0/eV_to_hartree;

    // parametry układu
    const double m=0.017;
    const double xmin=-500.0*nm_to_bohr;
    // const double xmax=-xmin;
    const double ymin=-350.0*nm_to_bohr;
    const double E=0.015*eV_to_hartree;
    const double ymax=-ymin;
    double Vparam[]={
        300.0*nm_to_bohr, // sigma_x
        300.0*nm_to_bohr, // sigma_y
        -1.3*eV_to_hartree // V_gates
    };
    Vparam[0]*=1; // żeby -Werror puściło bez potencjału

    // parametry symulacji
    const int Nx=int(50*0.6)-1;
    const int Ny=int(35*0.6)-1;

    // parametry wtórne
    const double dx=(ymax-ymin)/(Ny+1);
    // const double dy=(ymax-ymin)/(Ny-1);
    const double alpha=1.0/(2.0*m*dx*dx);

    // alokacja / inicjalizacja
    double *x=new double[Nx];
    for (int i=0;i<Nx;i++) x[i]=xmin+(i+1)*dx;
    double *y=new double[Ny];
    for (int j=0;j<Ny;j++) y[j]=ymin+(j+1)*dx;
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

    cmp *evals = new cmp[2*Ny];
    cmp *evecs = new cmp[2*Ny*Ny];
    upw_inv(E,alpha,Ny,evals,evecs);

    // normowanie funckji wlasnych
    double sum, invsqrtsum;
    for (int i=0;i<2*Ny;i++){
        sum=0.0;
        for (int j=0;j<Ny;j++) sum += std::norm(evecs[i*Ny+j]);
        invsqrtsum=1.0/std::sqrt(sum);
        for (int j=0;j<Ny;j++) evecs[i*Ny+j]*=invsqrtsum;
    }

    // znalezienie modów poprzecznych
    int liczba=0;
    bool *poprzeczne = new bool[2*Ny]{};
    double eps=1e-8;
    for (int i=0;i<2*Ny;i++){
        if (std::abs(evals[i])>=1.0-eps && std::abs(evals[i])<=1.0+eps){
            poprzeczne[i]=true;
            liczba++;
        }
    }
    liczba/=2;
    printf("\tLiczba modów poprzecznych: %d\n",liczba);

    // zapisanie
    FILE *lfile=fopen("lfile.csv","w");
    FILE *ufile=fopen("ufile.csv","w");
    FILE *poprzecznefile=fopen("poprzecznefile.csv","w");
    for (int i=0;i<2*Ny;i++){
        if (i!=0) fprintf(lfile,",");
        fprintf(lfile,"%e",std::abs(evals[i]));
        if (i!=0) fprintf(poprzecznefile,",");
        fprintf(poprzecznefile,"%d",poprzeczne[i]);
        for (int j=0;j<Ny;j++){
            fprintf(ufile,"%e,%e\n",std::real(evecs[i*Ny+j]),std::imag(evecs[i*Ny+j]));
        }
    }
    fclose(lfile);
    fclose(ufile);
    fclose(poprzecznefile);

    // zapisanie tylko ujemnych wartosci i wektorow
    cmp *evals_pop = new cmp[liczba];
    cmp *evecs_pop = new cmp[liczba*Ny];

    int l=0;
    for (int i=0;i<2*Ny;i++){
        if (poprzeczne[i] && std::imag(evals[i])<0.0){
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

    // petla po Vg 
    FILE *TRfile=fopen("TRfile.csv","w");
    for (double Vg=-1.3*eV_to_hartree;Vg<(-0.7+eps)*eV_to_hartree;Vg+=0.1*eV_to_hartree){
    auto t0=std::chrono::high_resolution_clock::now();

    // obliczenie psi
    double *V=new double[Nx*Ny]{};
    Vparam[2]=Vg;
    initV_QPC(V,Nx,Ny,x,y,ymin,ymax,Vparam);
    cmp **psin = new cmp*[liczba];
    for (int l=0;l<liczba;l++) psin[l] = new cmp[Nx*Ny];

    double c_in=1.0;

    for (int l=0;l<liczba;l++){
        znajdzPsi(psin[l],alpha,c_in,evecs_pop,l,Nx,Ny,evals_pop,liczba,E,V);
    }

    // zapisanie psi
    // FILE *psifile = fopen("psifile.csv","w");
    // for (int l=0;l<liczba;l++){
    //     for (int i=0;i<Nx;i++){
    //         for (int j=0;j<Ny;j++){
    //             if (i!=0 || j!=0) fprintf(psifile,",");
    //             fprintf(psifile,"%lf",std::norm(psin[l][i*Ny+j]));
    //         }
    //     }
    //     fprintf(psifile,"\n");
    // }
    // fclose(psifile);

    // wspolczynnik transmisji

    // wsp c_out d_out
    cmp *c_out = new cmp[liczba]{};
    cmp *d_out = new cmp[liczba]{};
    double *Tn = new double[liczba]{};
    double *Rn = new double[liczba]{};
    double T{}, R{};

    // wzory od 12 do 15
    for (int n=0;n<liczba;n++){
        for (int m=0;m<liczba;m++){
            c_out[m]=0.0; d_out[m]=0.0;
            for (int ni=0;ni<Ny;ni++){
                c_out[m] += std::conj(evecs_pop[m*Ny+ni])*psin[n][ni];
                d_out[m] += std::conj(evecs_pop[m*Ny+ni])*psin[n][(Nx-1)*Ny+ni];
            }
            c_out[m] -= cmp(deltaK(n,m));
            Tn[n] += std::norm(d_out[m]/c_in)*std::abs(v[m]/v[n]);
            Rn[n] += std::norm(c_out[m]/c_in)*std::abs(v[m]/v[n]);
        }
        T+=Tn[n];
        R+=Rn[n];
    }
    auto t1=std::chrono::high_resolution_clock::now();
    std::cout<<"\tczas: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<'\n';

    fprintf(TRfile,"%lf,%e,%e\n",Vg*hartree_to_eV,T,R);
    delete [] V;
    delete [] c_out;
    delete [] d_out;
    delete [] Tn;
    delete [] Rn;
    for (int l=0;l<liczba;l++) delete [] psin[l];
    delete [] psin;
    }
    fclose(TRfile);

    // zsumowane prawdopodobienstwo
    // double *rho = new double[Nx*Ny];
    // FILE *rhofile=fopen("rhofile.csv","w");
    // for (int i=0;i<Nx;i++){
    //     for (int j=0;j<Ny;j++){
    //         rho[i*Ny+j]=0.0;
    //         for (int n=0;n<liczba;n++){
    //             rho[i*Ny+j] += std::norm(psin[n][i*Ny+j]);
    //         }
    //         if (i!=0 || j!=0) fprintf(rhofile,",");
    //         fprintf(rhofile,"%lf",rho[i*Ny+j]);
    //     }
    // }
    // fclose(rhofile);

    // zapisanie potecjalu
    // FILE *potfile=fopen("potfile.csv","w");
    // for (int i=0;i<Nx;i++){
    //     for (int j=0;j<Ny;j++){
    //         if (i!=0 || j!=0) fprintf(potfile,",");
    //         fprintf(potfile,"%lf",V[i*Ny+j]*hartree_to_eV);
    //     }
    // }
    // fclose(potfile);

    // zapisanie parametrów
    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%lf,%d,%d,%e,%e,%d",E*hartree_to_eV,Ny,Nx,xmin*bohr_to_nm,ymin*bohr_to_nm,liczba);
    fclose(misc);

    // czystki 
    delete [] x;
    delete [] y;
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
