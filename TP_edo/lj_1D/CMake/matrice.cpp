#include "matrice.hpp"


/////////////////////////////////////
// METHODES de la classe vecteur


vec::vec()
// constructeur par defaut
{
    dim = 0;
    coord = NULL;
}


vec::vec(int d)
//construit un vec de taille donnÈe
{
    dim = d;
    coord = new double[d];
}


vec::vec(const vec & v)
// operateur de recopie de v dans le courant
{
    dim = v.dim;
    if (v.coord!=NULL) {
        coord = new double[dim];
        for (int i=0;i<dim;i++) coord[i]=(v.coord[i]);
    } else coord=NULL;
}


vec::~vec()
// destruction
{
    if (coord!=NULL) delete[] coord;
    coord = NULL;
    dim =0;
    
}


vec & vec::operator= (const vec & a)
// affectation surchargÈe
{
    if (this!=&a) {
        dim = a.dim;
        if (coord!=NULL) delete[] coord;
        if (a.coord!=NULL) {
            coord = new double[dim];
            for (int i=0;i<dim;i++) coord[i]=a.coord[i];
        } else {
            coord = NULL;
        }
    }
    return *this;
}


vec vec::operator+ (const vec & a)
// addition surchargÈe
{
    vec b;
    
    if (dim!=a.dim) {
        std::cerr<<"addition de vecs de tailles differentes impossible"<<std::endl;
        exit(1);
    } else {
        if ((coord != NULL) && (a.coord!=NULL)) {
            b.dim = a.dim;
            b.coord = new double[dim];
            for (int i=0;i<dim;i++) b.coord[i] = coord[i] + a.coord[i];
        }
    }
    return b;
}

vec vec::operator- (const vec & a)
// soustraction surchargÈe
{
    vec b;
    
    if (dim!=a.dim) {
        std::cerr<<"soustraction de vecs de tailles differentes impossible"<<std::endl;
        exit(1);
    } else {
        if ((coord != NULL) && (a.coord!=NULL)) {
            b.dim = a.dim;
            b.coord = new double[dim];
            for (int i=0;i<dim;i++) b.coord[i] = coord[i] - a.coord[i];
        }
    }
    return b;
}

vec vec::operator* (const double a)
// multiplication par un scalaire
{
    vec b;
    
    if (coord != NULL) {
        b.dim = dim;
        b.coord = new double[dim];
        for (int i=0;i<dim;i++) b.coord[i] = coord[i]*a;
    }
    
    return b;
}

double & vec::operator() (int a)
// renvoie la composante d'un vec
{
    
    if (dim-1<a) {
        std::cerr<<"on demande une composante inexistante"<<std::endl;
        exit(1);
    } else {
        if (coord!=NULL) {
            return coord[a];
        } else {
            std::cerr<<"on demande la composante d'un vec vide"<<std::endl;
            exit(1);
        }
    }
}

double & vec::operator[] (int a)
// idem que ()
{
    
    if (dim-1<a) {
        std::cerr<<"on demande une composante inexistante"<<std::endl;
        exit(1);
    } else {
        if (coord!=NULL) {
            return coord[a];
        } else {
            std::cerr<<"on demande la composante d'un vec vide"<<std::endl;
            exit(1);
        }
    }
}


void vec::set_size(int d)
// met le vec courant ‡ la longueur souhaitÈe
{
    dim = d;
    if (coord!=NULL) delete[] coord;
    coord = new double[dim];
}

void vec::zeros()
// met le vec courant ‡ 0
{
    if (coord!=NULL) {
        for (int i =0;i<dim;i++) coord[i]=0.0;
    }
}

//ostream& vec::affiche (ostream& os) const
//{
//  for (int i=0;i<dim;i++){
//	os << coord[i]<<" ";
//  }

//  return os;
//}


/////////////////////////////////////////////////
////// METHODES de la classe mat

mat::mat()
// constructeur par defaut
{
    nligne = 0;
    ncol = 0;
    col = NULL;
}


mat::mat(int d1, int d2)
// constructeur avec choix de la dimension
{
    nligne = d1;
    ncol = d2;
    col = new vec[d2];
    for (int i=0;i<d2;i++) {
        // on a construit les vecs: il reste ‡ les mettre
        // a la bonne dimansion et a les initialiser
        (col[i]).set_size(d1);
    }
}


mat::mat(const mat & M)
// operateur de recopie
{
    nligne = M.nligne;
    ncol = M.ncol;
    if (M.col!=NULL) {
        col = new vec[ncol];
        for (int i=0;i<ncol;i++) {
            col[i] = M.col[i]; // utilisation de l'affectattion entre vec
        }
    } else col=NULL;
}


mat::~mat()
// destruction
{
    if (col!=NULL) {
        for (int i=0;i<ncol;i++) col[i].~vec();
    }
    delete[] col;
    col = NULL;
    ncol = 0;
    nligne = 0;
}


mat & mat::operator= (const mat & M)
// affectation surchargÈe
{
    if (this!=&M) {
        nligne = M.nligne;
        ncol = M.ncol;
        if (col!=NULL) for (int i=0;i<ncol;i++) col[i].~vec();
        if (M.col!=NULL) {
            col = new vec[ncol];
            for (int i=0;i<ncol;i++) col[i]=M.col[i]; // affectattion entre vec
        } else {
            col = NULL;
        }
    }
    return *this;
}


mat mat::operator+ (const mat & M)
// addition surchargÈe
{
    mat b;
    
    if ((nligne!=M.nligne) || (ncol!=M.ncol)) {
        std::cerr<<"addition de mats de tailles differentes impossible"<<std::endl;
        exit(1);
    } else {
        if ((col != NULL) && (M.col!=NULL)) {
            b.nligne = nligne;
            b.ncol = ncol;
            b.col = new vec[ncol];
            for (int i=0;i<ncol;i++) b.col[i] = col[i] + M.col[i];
        }
    }
    return b;
}

mat mat::operator* (const double a)
// multiplication par un scalaire
{
    mat b;
    
    if (col != NULL) {
        b.nligne = nligne;
        b.ncol = ncol;
        b.col = new vec[ncol];
        for (int i=0;i<ncol;i++) b.col[i] = col[i]*a;
    }
    return b;
}

double & mat::operator() (int a, int b) const
// renvoie la composante d'une mat
{
    
    if ((nligne-1<a) || (ncol-1<b)) {
        std::cerr<<"on demande une composante inexistante"<<std::endl;
        std::cerr<<"nb ligne="<<nligne<<", num ligne="<<a<<", nb col="<<ncol<<", num col="<<b<<std::endl;
        exit(1);
    } else {
        if (col!=NULL) {
            return(col[b](a));
        } else {
            std::cerr<<"on demande la composante d'une mat vide"<<std::endl;
            exit(1);
        }
    }
}

void mat::set_size(int d1, int d2)
// met le mat courant ‡ la longueur souhaitÈe
{
    nligne = d1;
    ncol = d2;
    
    col = new vec[d2];
    for (int i=0;i<d2;i++) {
        // on a construit les vecs: il reste ‡ les mettre
        // a la bonne dimension et a les initialiser
        (col[i]).set_size(d1);
    }
}

void mat::zeros()
// met le mat courant ‡ 0
{
    if (col!=NULL) {
        for (int i =0;i<ncol;i++) col[i].zeros();
    }
}

vec mat::get_col(int a)
// on renvoie une des colonnes
{  
    vec colonne = col[a];
    return colonne;
}

void mat::set_col(int a , vec v)
// on impose uen colonne dans la matrice
{
    col[a]=v;
}

//ostream& mat::affiche (ostream& os) const
//{
//  for (int i=0;i<nligne;i++){
//	for (int j=0;j<ncol;j++) {
//	  os << col[j](i)<<" ";
//	}
//  }

//  return os;
//}

//ostream& operator <<(ostream& os, const vec & d)
//{
//  return d.affiche(os);
//}

//ostream& operator <<(ostream& os, const mat & d)
//{
//  return d.affiche(os);
//}
