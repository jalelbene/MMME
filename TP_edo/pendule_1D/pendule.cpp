// Etude d'un pendule 1D

// make

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "matrice.hpp"

using namespace std;

// pas de temps
const double dt = 0.2;

// nb de pas de temps
const int pas = 1000;

// frequence d'affichage
const int freq = 1;

// tolerance et nb max d'iterations dans Kunger Kutta implicite
const double tol = 1.e-9;
const int NiterMax= 50;

struct Particule
{
  double q;
  double p;
};


//------------------------
//  potentiels et forces
//------------------------

double V (double pos)
{
  double pot = -cos(pos) + 0.2*sin(2.*pos);
  return pot;
}

double f (double pos)
{
  double force = -sin(pos) - 0.4*cos(2.*pos);
  return force;
}

// hamiltonien
double H (double pos, double imp)
{
  double pot = V(pos);
  pot += 0.5*pow(imp,2);
  return pot;
}

// on ecrit l'edo hamiltonienne sous la forme \dot{x} = Gamma(x)
// avec x=(q,p) et on definit ici le champ de vecteur Gamma(x)
// RHS = Right-Hand Side
vec RHS (vec x)
{
  vec rhs_(2);
  rhs_(0) = x(1);
  rhs_(1) = f(x(0));
  
  return rhs_;
}
  

//--------------------------
//  algorithme Verlet
//-------------------------

Particule Verlet(Particule X)
{
  Particule Y;
  Y.p = X.p + f(X.q)*(0.5*dt);
  Y.q = X.q + Y.p*dt;
  Y.p = Y.p + f(Y.q)*(0.5*dt); 
  return Y;
}

//--------------------------
//  algorithme RK Lobatto III B d'ordre 4
//-------------------------
Particule Lobatto(Particule X)
{
  vec k1(2),k2(2),k3(2);
  vec k_old1(2),k_old2(2),k_old3(2);
  vec current(2);
  
  // Resolution du probleme non lineaire en les k:

  // initialisation
  current(0) = X.q;
  current(1) = X.p;
  k1 = RHS(current);
  k2 = k1;
  k3 = k1;
  double diff = 1.;
  int niter = 0;
  // iteration sur les k
  while ((diff > (tol*tol)) && (niter < NiterMax)) {
	k_old1 = k1;
	k_old2 = k2;
	k_old3 = k3;
      
    vec Z(2);
    Z(0) = X.q + k_old1(0)*(dt/6.) + k_old2(0)*(-dt/6.) + k_old3(0)*(0.*dt);
    Z(1) = X.p + k_old1(1)*(dt/6.) + k_old2(1)*(-dt/6.) + k_old3(1)*(0.*dt);
	k1 = RHS(Z); // a completer
      
    Z(0) = X.q + k_old1(0)*(dt/6.) + k_old2(0)*(dt/3.) + k_old3(0)*(0.*dt);
    Z(1) = X.p + k_old1(1)*(dt/6.) + k_old2(1)*(dt/3.) + k_old3(1)*(0.*dt);
	k2 = RHS(Z); // a completer
      
    Z(0) = X.q + k_old1(0)*(dt/6.) + k_old2(0)*(5.*dt/6.) + k_old3(0)*(0.*dt);
    Z(1) = X.p + k_old1(1)*(dt/6.) + k_old2(1)*(5.*dt/6.) + k_old3(1)*(0.*dt);
	k3 = RHS(Z); // a completer
      
	diff = 0.;
	for (int a=0;a<2;a++) {
	  diff += pow(k1(a)-k_old1(a),2) + pow(k2(a)-k_old2(a),2) + pow(k3(a)-k_old3(a),2);
	}
	niter += 1;
  }
  
  if (niter == NiterMax) {
	cout<<"Nb max d'iterations atteint dans Lobatto"<<endl;
  }
	
  Particule Y;
  Y.q = X.q + k1(0)*(dt/6.) + k2(0)*(2.*dt/3.) + k3(0)*(dt/6.);
  Y.p = X.p + k1(1)*(dt/6.) + k2(1)*(2.*dt/3.) + k3(1)*(dt/6.);

  return Y;
}





//----------------------------
// integration des equations
//---------------------------

int main () 
{
  Particule X;

  // choix des CI
  double pos=0.;
  X.q = pos;

  double imp = 2.2;
  X.p = imp;

  int tracage = freq;
  ofstream energie("energie");
  ofstream position("position");

  double H_init = H(X.q,X.p);
	
  position << 0. << " "<< X.q << " " << X.p << endl;
  energie << 0. << "  " << H(X.q,X.p)-H_init <<endl;

  
  // integration
  for (int i = 0; i < pas; i++) {
//	X = Verlet(X);
	X = Lobatto(X);
	
	
	if (tracage == freq) {
	  tracage = 0;
	  position << (i+1)*dt << " "<< X.q << "  " << X.p << endl;
	  energie << (i+1)*dt << "  " << H(X.q,X.p)-H_init <<endl;
	}
	tracage += 1;
  }  

    return 0;
}
