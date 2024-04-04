#ifndef FILM_H
#define FILM_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <time.h> 
#include <stdlib.h>
#include <array>
using namespace std;

#define mu_0 4*M_PI*1e-7   //Perméabilité du vide
//Dimensions de l'ellipsoïde de révolution
#define lx 0.5*1e-9
#define ly 0.5*1e-9
#define lz 1*1e-9

class Film {
private:           
    double 
        M_norm,   //Norme de l'aimantation (la même dans chaque grain) 
        theta_def;   //Valeur de theta par défaut 
    double* theta_data;   //Tableau des valeurs de theta   
public:
    unsigned int 
        nx,   //Nombre de colonnes 
        ny;   //Nombre de lignes
    string name;   //Nom du film
    //Constructors
    Film(string name_, unsigned int nx_, unsigned int ny_, double M_norm_, double theta_def_);
    Film(string name_, unsigned int nx_, unsigned int ny_, double M_norm_=1);   //Génère un film avec des directions de M aléatoires (random constructor)
    Film(const Film&);   //Copy Constructor
    //Destructor
    ~Film();
    //Operators
    //Accéder à un élément du film
    double operator()(int i) const;   //En utilisant l'indice i
    double operator()(int x, int y) const;   //En utilisant les coordonnées x et y
    //Modifier un élément du film   
    double& operator()(int i);   //En utilisant l'indice i
    double& operator()(int x, int y);   //En utilisant les coordonnées x et y 
    //Methods
    unsigned int find_i(int x, int y);   //Trouve l'indice à partir des coordonnées x et y dans le film
    array<unsigned int, 2> find_xy(int i);   //Trouve (x,y) à partir de l'indice i
    void show();
    void write(string fich_name="theta_data");
    //Energies
    double E_ech(double A);   //Energie d'échange (version intégrale)
    double E_d(double Nxx, double Nyy, double Nzz);   //Energie d'interaction avec le champ démagnétisant
    double E_ani(double K);   //Energie d'anisotropie
    double E_z(double phi, double H_a_norm);   //Energie Zeeman. Phi est l'angle entre l'axe Oz et H_a
    double E_tot(double A, double Nxx, double Nyy, double Nzz, double K, double phi, double H_a_norm);   //Energie totale
    double E_tot(double A, double Nxx, double Nyy, double Nzz, double K);   //Energie totale sans la présence d'un champ magnétique extérieur
    //Advanced writing
    void write(string fich_name, double A, double Nxx, double Nyy, double Nzz, double K, double phi, double H_a_norm); 
};

#endif
