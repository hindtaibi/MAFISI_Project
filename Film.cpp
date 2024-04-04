#include "Film.h"

//Constructors
Film::Film(string name_, unsigned int nx_, unsigned int ny_, double M_norm_, double theta_def_) {
    name=name_;
    nx=nx_;
    ny=ny_;
    M_norm=M_norm_;
    theta_def=theta_def_;
    theta_data=new double[nx*ny];
    for (unsigned int i=0; i<nx*ny; i++) theta_data[i]=theta_def;
    cout<<"Film "<<name<<" created @ "<<theta_data<<endl;
}

Film::Film(string name_, unsigned int nx_, unsigned int ny_, double M_norm_) {   //Random constructor
    name=name_;
    nx=nx_;
    ny=ny_;
    M_norm=M_norm_;
    theta_data=new double[nx*ny];
    srand48(time(NULL));
    for (unsigned int i=0; i<nx*ny; i++) theta_data[i]=drand48()*2*M_PI;
    cout<<"Film "<<name<<" created @ "<<theta_data<<endl;
}

Film::Film(const Film& other)   //Copy constructor
    : name("copy of "+other.name), ny(other.ny), nx(other.nx), M_norm(other.M_norm), theta_data(nullptr) {
        theta_data=new double[ny*nx];
        for (unsigned int k = 0; k < ny*nx; k++) theta_data[k]=other.theta_data[k];
        cout<<"Film "<<name<<" created @ "<<theta_data<<endl;
    }

//Destructor
Film::~Film() {
    if (theta_data!=nullptr) delete[] theta_data;
    std::cout <<"Film "<<name<<" @ "<<theta_data<<" destroyed"<<endl;
    theta_data = nullptr;
}

//Operators
double Film::operator()(int i) const{  
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);    
    return theta_data[i];
}

double Film::operator()(int x, int y) const {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    return theta_data[nx*(y%ny)+x%nx];
}

double& Film::operator()(int i) {
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);   
    return theta_data[i];
}

double& Film::operator()(int x, int y) {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    return theta_data[nx*(y%ny)+x%nx];
}

//Methods
unsigned int Film::find_i(int x, int y) {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    unsigned int i=nx*(y%ny)+x%nx;  
    return i;    
}

array<unsigned int,2> Film::find_xy(int i) {
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);
    unsigned int x=i*(1-nx)/(1-nx*ny);
    unsigned int y=i*(1-ny)/(1-nx*ny);
    return {x,y};
}

void Film::show() {
    for (unsigned int y=0; y<ny; y++) {
        for (unsigned int x=0; x<nx; x++) cout<<theta_data[find_i(x,y)]<<' ';
        cout<<endl;
    }
    cout<<endl;
}

void Film::write(string fich_name) {
    fstream fich;
    fich.open(fich_name, ios::out);
    for(unsigned int y=0; y<ny; y++) {
        for(unsigned int x=0; x<nx; x++) fich<<theta_data[find_i(x,y)]<<' ';
        fich<<endl;
    }
    fich.close();
}

//Energies
double Film::E_ech(double A) {   
    double E=0.;
    for (unsigned int y=0; y<ny; y++) {
        for (unsigned int x=0; x<nx; x++) {
            E+=(sin(theta_data[find_i(x+1,y)])-sin(theta_data[find_i(x,y)]))*(sin(theta_data[find_i(x+1,y)])-sin(theta_data[find_i(x,y)]))/(4*lx*lx)
               +(sin(theta_data[find_i(x,y+1)])-sin(theta_data[find_i(x,y)]))*(sin(theta_data[find_i(x,y+1)])-sin(theta_data[find_i(x,y)]))/(4*ly*ly)
                +(cos(theta_data[find_i(x+1,y)])-cos(theta_data[find_i(x,y)]))*(cos(theta_data[find_i(x+1,y)])-cos(theta_data[find_i(x,y)]))/(4*lx*lx)
                 +(cos(theta_data[find_i(x,y+1)])-cos(theta_data[find_i(x,y)]))*(cos(theta_data[find_i(x,y+1)])-cos(theta_data[find_i(x,y)]))/(4*ly*ly);
        }
    }
    return A*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film::E_d(double Nxx, double Nyy, double Nzz) {
    if (Nxx+Nyy+Nzz!=1.) cout<<"Warning: Nxx+Nyy+Nzz must be equal to 1"<<endl;
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) E+=M_norm*M_norm*(Nxx*sin(theta_data[i])*sin(theta_data[i])+Nzz*cos(theta_data[i])*cos(theta_data[i]));
    return (mu_0/2)*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film::E_ani(double K) {
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) E+=cos(theta_data[i])*cos(theta_data[i]);
    return K*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film::E_z(double phi, double H_a_norm) {
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) E+=M_norm*H_a_norm*cos(phi-theta_data[i]);
    return -mu_0*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film::E_tot(double A, double Nxx, double Nyy, double Nzz, double K, double phi, double H_a_norm) {
    double E=E_ech(A)+E_d(Nxx, Nyy, Nzz)+E_ani(K)+E_z(phi, H_a_norm);
    return E;
}

double Film::E_tot(double A, double Nxx, double Nyy, double Nzz, double K) {
    double E=E_ech(A)+E_d(Nxx, Nyy, Nzz)+E_ani(K);
    return E;
}

//Advanced writing
void Film::write(string fich_name, double A, double Nxx, double Nyy, double Nzz, double K, double phi, double H_a_norm) {
    fstream fich;
    fich.open(fich_name, ios::out);
    fich<<M_norm<<endl;
    fich<<A<<endl;
    fich<<Nxx<<endl;
    fich<<Nyy<<endl;
    fich<<Nzz<<endl;
    fich<<K<<endl;
    fich<<phi<<endl;
    fich<<H_a_norm<<endl;
    fich<<E_tot(A, Nxx, Nyy, Nzz, K, phi, H_a_norm)<<endl;
    for(unsigned int y=0; y<ny; y++) {
        for(unsigned int x=0; x<nx; x++) fich<<theta_data[find_i(x,y)]<<' ';
        fich<<endl;
    }
    fich.close();
}
