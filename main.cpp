#include "Film.h" 
#include "Film3D.h"

Film F_min(Film F, double A, double Nxx, double Nyy, double Nzz, double K, double phi, double H_a_norm, unsigned int N) {   //C'est la fonction qui minimise l'énergie d'un film d'aimantation 2D
    srand48(time(NULL));
    for (unsigned int n=0; n<N; n++) {
        unsigned int i=rand()%(F.nx*F.ny);
        double E_ini=F.E_tot(A, Nxx, Nyy, Nzz, K, phi, H_a_norm);
        double theta_ini=F(i);
        F(i)=drand48()*2*M_PI;
        if (F.E_tot(A, Nxx, Nyy, Nzz, K, phi, H_a_norm)>=E_ini) F(i)=theta_ini;
    }
    return F;
}

Film3D F_min3D(Film3D F, double A, double Nxx, double Nyy, double Nzz, double K, double phiz, double phix, double H_a_norm, unsigned int N) {
    srand48(time(NULL));
    for (unsigned int n=0; n<N; n++) {
        unsigned int i=rand()%(F.nx*F.ny);
        double E_ini=F.E_tot(A, Nxx, Nyy, Nzz, K, phiz, phix, H_a_norm);
        double thetaz_ini=real(F(i));
	double thetax_ini=imag(F(i));
        F(i)=drand48()*2*M_PI+I*thetax_ini;
        if (F.E_tot(A, Nxx, Nyy, Nzz, K, phiz, phix, H_a_norm)>=E_ini) {
		F(i)=real(F(i))+I*(-drand48()*M_PI/2+M_PI/2);
        }
        else {
            double E_ini=F.E_tot(A, Nxx, Nyy, Nzz, K, phiz, phix, H_a_norm);
	    F(i)=real(F(i))+I*drand48()*M_PI;
            if (F.E_tot(A, Nxx, Nyy, Nzz, K, phiz, phix, H_a_norm)>=E_ini) F(i)=real(F(i))+I*thetax_ini;
	}
    }
    return F;
}

int main() {
    //FePt
    double 
        Nzz=0.3,
        Nxx=(1-Nzz)/2,
        Nyy=Nxx,
        M_norm=1.14*1e6, 
        A=2.7*1e-11, 
        K=6.6*1e6, 
        H_a_norm=0,
        phi=0;
    unsigned int  
        N=1*1e5;   //Nombre d'itérations de la fonction F_min
    Film F("F", 20, 20, M_norm);
    Film Fm=F_min(F, A, Nxx, Nyy, Nzz, K, 0, 0, N);
    Fm.write("FePt2020.txt");   //Fichier de données contenant les directions d'aimantation pour lesquelles l'énergie est minimale

    return 0;
}
