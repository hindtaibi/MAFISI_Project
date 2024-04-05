#include "Film3D.h"

//Constructors
Film3D::Film3D(string name_, unsigned int nx_, unsigned int ny_, double M_norm_,complex<double> theta_def_) {
    name=name_;
    nx=nx_;
    ny=ny_;
    M_norm=M_norm_;
    theta_def=theta_def_;
    theta_data=new complex<double>[nx*ny];
    for (unsigned int i=0; i<nx*ny; i++) theta_data[i]=theta_def;
    cout<<"Film "<<name<<" created @ "<<theta_data<<endl;
}

Film3D::Film3D(string name_, unsigned int nx_, unsigned int ny_, double M_norm_) {   //Random constructor
    name=name_;
    nx=nx_;
    ny=ny_;
    M_norm=M_norm_;
    theta_data=new complex<double>[nx*ny];
    srand48(time(NULL));
    for (unsigned int i=0; i<nx*ny; i++) {
        double thetaz=drand48()*2*M_PI;
        double thetax=drand48()*M_PI;
        theta_data[i]=thetaz+I*thetax;
    }
    cout<<"Film3D "<<name<<" created @ "<<theta_data<<endl;
}

Film3D::Film3D(const Film3D& other)   //Copy constructor
    : name("copy of "+other.name), ny(other.ny), nx(other.nx), M_norm(other.M_norm), theta_data(nullptr) {
        theta_data=new complex<double>[ny*nx];
        for (unsigned int k = 0; k < ny*nx; k++) theta_data[k]=other.theta_data[k];
        cout<<"Film3D "<<name<<" created @ "<<theta_data<<endl;
    }

//Destructor
Film3D::~Film3D() {
    if (theta_data!=nullptr) delete[] theta_data;
    std::cout <<"Film3D "<<name<<" @ "<<theta_data<<" destroyed"<<endl;
    theta_data = nullptr;
}

//Operators
complex<double> Film3D::operator()(int i) const{  
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);    
    return theta_data[i];
}

complex<double> Film3D::operator()(int x, int y) const {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    return theta_data[nx*(y%ny)+x%nx];
}

complex<double>& Film3D::operator()(int i) {
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);  
	return theta_data[i];
}

complex<double>& Film3D::operator()(int x, int y) {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    return theta_data[nx*(y%ny)+x%nx];
}

//Methods
unsigned int Film3D::find_i(int x, int y) {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    unsigned int i=nx*(y%ny)+x%nx;  
    return i;    
}

array<unsigned int,2> Film3D::find_xy(int i) {
    if (i<0) i+=nx*ny;
    i=i%(nx*ny);
    unsigned int x=i*(1-nx)/(1-nx*ny);
    unsigned int y=i*(1-ny)/(1-nx*ny);
    return {x,y};
}

void Film3D::show() {
    for (unsigned int y=0; y<ny; y++) {
        for (unsigned int x=0; x<nx; x++) cout<<theta_data[find_i(x,y)]<<' ';
        cout<<endl;
    }
    cout<<endl;
}

void Film3D::write(string fichz_name, string fichy_name) {
    fstream fichz;
    fichz.open(fichz_name, ios::out);
    for(unsigned int y=0; y<ny; y++) {
        for(unsigned int x=0; x<nx; x++) fichz<<real(theta_data[find_i(x,y)])<<' ';
        fichz<<endl;
    }
    fichz.close();

    fstream fichy;
    fichy.open(fichy_name, ios::out);
    for(unsigned int y=0; y<ny; y++) {
        for(unsigned int x=0; x<nx; x++) fichy<<imag(theta_data[find_i(x,y)])<<' ';
        fichy<<endl;
    }
    fichy.close();
}

//Energies
double Film3D::E_ech(double A) {   
    double E=0.;
    for (unsigned int y=0; y<ny; y++) {
        for (unsigned int x=0; x<nx; x++) {
            double
                thetaz=real(theta_data[find_i(x,y)]),
                thetazx=real(theta_data[find_i(x+1,y)]),
                thetazy=real(theta_data[find_i(x,y+1)]),
                thetax=imag(theta_data[find_i(x,y)]),
                thetaxx=imag(theta_data[find_i(x+1,y)]),
                thetaxy=imag(theta_data[find_i(x,y+1)]);
            E+=(sin(thetazx)*cos(thetaxx)-sin(thetaz)*cos(thetax))*(sin(thetazx)*cos(thetaxx)-sin(thetaz)*cos(thetax))/(4*lx*lx)
               +(sin(thetazy)*cos(thetaxy)-sin(thetaz)*cos(thetax))*(sin(thetazy)*cos(thetaxy)-sin(thetaz)*cos(thetax))/(4*ly*ly)
                +(sin(thetaxx)-sin(thetax))*(sin(thetaxx)-sin(thetax))/(4*lx*lx)
                 +(sin(thetaxy)-sin(thetax))*(sin(thetaxy)-sin(thetax))/(4*ly*ly)
                  +(cos(thetazx)*cos(thetaxx)-cos(thetaz)*cos(thetax))*(cos(thetazx)*cos(thetaxx)-cos(thetaz)*cos(thetax))/(4*lx*lx)
                   +(cos(thetazy)*cos(thetaxy)-cos(thetaz)*cos(thetax))*(cos(thetazy)*cos(thetaxy)-cos(thetaz)*cos(thetax))/(4*ly*ly);
        }
    }
    return A*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film3D::E_d(double Nxx, double Nyy, double Nzz) {
    if (Nxx+Nyy+Nzz!=1.) cout<<"Warning: Nxx+Nyy+Nzz must be equal to 1"<<endl;
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) {
	double thetaz=real(theta_data[i]);
        double thetax=imag(theta_data[i]);
	E+=M_norm*M_norm*(Nxx*sin(thetaz)*cos(thetax)*sin(thetaz)*cos(thetax)
                          +Nyy*sin(thetax)*sin(thetax)
                           +Nzz*cos(thetaz)*cos(thetax)*cos(thetaz)*cos(thetax));
	}
    return (mu_0/2)*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film3D::E_ani(double K) {
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) {
	double thetaz=real(theta_data[i]);
        double thetax=imag(theta_data[i]);
	E+=cos(thetax)*cos(thetaz)*cos(thetax)*cos(thetaz);
    }
    return K*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film3D::E_z(double phiz, double phix, double H_a_norm) {
    double E=0.;
    for (unsigned int i=0; i<nx*ny; i++) {
	double thetaz=real(theta_data[i]);
        double thetax=imag(theta_data[i]);
	E+=M_norm*H_a_norm*(sin(thetaz)*cos(thetax)*sin(phiz)*cos(phix)
                            +sin(thetax)*sin(phix)
                             +cos(thetaz)*cos(thetax)*cos(phiz)*cos(phix));
	}
    return -mu_0*E*(4/3)*M_PI*lx*ly*lz*nx*ny;
}

double Film3D::E_tot(double A, double Nxx, double Nyy, double Nzz, double K, double phiz, double phix, double H_a_norm) {
    double E=E_ech(A)+E_d(Nxx, Nyy, Nzz)+E_ani(K)+E_z(phiz, phix, H_a_norm);
    return E;
}
