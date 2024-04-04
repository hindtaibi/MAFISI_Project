#ifndef FILM3D_H
#define FILM3D_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <array>
#include <time.h> 
#include <stdlib.h>
#include <complex.h>
using namespace std;

#define mu_0 4*M_PI*1e-7   //Perméabilité du vide
//Dimensions de l'ellipsoïde de révolution
#define lx 0.5*1e-9
#define ly 0.5*1e-9
#define lz 1*1e-9

class Film3D {
private:           
    double M_norm;  //Norme de l'aimantation (la même dans chaque grain) 
    complex<double> theta_def;   //Valeur de theta par défaut (thetaz, thetax)
    complex<double>* theta_data;   //Tableau des valeurs de theta (thetaz, thetax)	  
public:
    unsigned int 
        nx,   //Nombre de colonnes 
        ny;   //Nombre de lignes
    string name;   //Nom du film 
    //Constructors
    Film3D(string name_, unsigned int nx_, unsigned int ny_, double M_norm_, complex<double> theta_def_);
    Film3D(string name_, unsigned int nx_, unsigned int ny_, double M_norm_=1);   //Génère un film avec des directions de M aléatoires (random constructor)
    Film3D(const Film3D&);   //Copy Constructor
    //Destructor
    ~Film3D();
    //Operators
    //Accéder à un élément du film
    complex<double> operator()(int i) const;   //En utilisant l'indice i
    complex<double> operator()(int x, int y) const;   //En utilisant les coordonnées x et y
    //Modifier un élément du film   
    complex<double>& operator()(int i);   //En utilisant l'indice i
    complex<double>& operator()(int x, int y);   //En utilisant les coordonnées x et y 
    //Methods
    unsigned int find_i(int x, int y);   //Trouve l'indice à partir des coordonnées x et y dans le film
    array<unsigned int, 2> find_xy(int i);   //Trouve (x,y) à partir de l'indice i
    void show();
    void write(string fichz_name="thetaz_data", string fichy_name="thetax_data");
    //Energies
    double E_ech(double A);   //Energie d'échange (version intégrale)
    double E_d(double Nxx, double Nyy, double Nzz);   //Energie d'interaction avec le champ démagnétisant
    double E_ani(double K);   //Energie d'anisotropie
    double E_z(double phiz, double phix, double H_a_norm);   //Energie Zeeman
    double E_tot(double A, double Nxx, double Nyy, double Nzz, double K, double phiz, double phix, double H_a_norm);   //Energie totale
};

#endif
