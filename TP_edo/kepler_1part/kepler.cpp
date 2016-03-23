// une particule dans un potentiel central gravitationnel

// pour compiler: make

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
const double dt = 0.001;

// nb de pas de temps
const int pas = 1000000;

// frequence d'affichage
const int freq = 1;

// dimension
const int dim = 3;

// tolerance pour projection sur la surface d'energie constante
const double tol = 1.e-14;

// nb d'iteration maximum pour projection sur la surface d'energie constante
const int NiterMax= 50;

struct Particule
{
  vec q;
  vec p;
};


//------------------------
//  energie potentielle, force et hamiltonien
//------------------------

double V (vec pos)
{
  double dis = 0.;
  for (int a=0; a<dim; a++) {
	dis += pow(pos(a),2);
  }
  dis = sqrt(dis);
  double pot = -1./dis;
  return pot;
}

vec f (vec pos)
{
  vec force(dim);
  double dis3 = 0.;

  for (int a=0; a<dim; a++) {
	dis3 += pow(pos(a),2);
  }
  dis3 = pow(sqrt(dis3),3);
  for (int a=0; a<dim; a++) {
	force(a) = -pos(a)/dis3;
  }
  
  return force;
}

// hamiltonien
double H (vec pos, vec imp)
{
  double pot = V(pos);
  for (int a=0; a<dim; a++) {
	pot += 0.5*pow(imp(a),2);
  }
  return pot;
}

//--------------------------
//  algorithme Symplectic Euler
//-------------------------

Particule SymplecticEuler(Particule X)
{
  Particule Y;
  Y.p = X.p + f(X.q)*dt;
  Y.q = X.q + Y.p*dt;
  return Y;
}


//--------------------------
//  algorithme Symplectic Euler + energie projection
//-------------------------

// fonction a annuler
double func(double q, double q0, double p0, double E)
{
  double res = 1. + q*q*q0-pow(q,3);
  res = 0.5*p0*p0/(res*res);
  res = res - 1./q - E;
  return res;
}

// derivee de func
double der(double q, double q0, double p0, double E)
{
  double res = 1. + q*q*q0-pow(q,3);
  res = -p0*p0/(pow(res,3));
  res = res*(2.*q*q0-3.*q*q);
  res = res + 1./(q*q);
  return res;
}

// solveur: newton
double solve(double q0, double p0, double E)
{
  
  double q = q0;
  double residu = func(q,q0,p0,E);
  int niter = 0;
  double pente;
  // la fonction abs sur un double semble poser probleme.
  double abs_residu = residu;
  if (residu < 0) {
	abs_residu = -residu;
  }
    
  while ((abs_residu > tol) && (niter < NiterMax)) {
	pente = der(q,q0,p0,E);
	q -= residu/pente;
	residu = func(q,q0,p0,E);
	abs_residu = residu;
	if (residu < 0) {
	  abs_residu = -residu;
	}
	niter += 1;
  }
  
  if (niter == NiterMax) {
	cout<<"Nb max d'iterations atteint dans Newton"<<endl;
  }
  
  return q;
}

Particule SymplecticEuler_proj(Particule X, double E_init)
{
  Particule X0 = SymplecticEuler(X);
  Particule Y = X0;
  vec p0 = X0.p;
  vec q0 = X0.q;
  double p0_mod = 0.;
  double q0_mod = 0.;
  for (int a=0; a< dim; a++) {
	p0_mod += p0(a)*p0(a);
	q0_mod += q0(a)*q0(a);
  }
  p0_mod = sqrt(p0_mod);
  q0_mod = sqrt(q0_mod);
  
  double q_mod = solve(q0_mod,p0_mod,E_init);
  double lambda = (q0_mod-q_mod)*pow(q_mod,2);
  for (int a=0; a< dim; a++) {
	Y.p(a) = p0(a)/(1.+lambda);
	Y.q(a) = q0(a)/(1.+lambda/pow(q_mod,3));
  }
  return Y;
}


  
//----------------------------
// integration des equations
//---------------------------

int main () 
{
    
  Particule X;

  // choix des CI
  vec pos(dim);
  pos(0)=1.;
  for (int a=1;a<dim;a++) {
	pos(a) = 0.;
  }
  X.q = pos;

  vec imp(dim);
  imp(0) = 0.;
  imp(1) = 0.8;
  for (int a=2;a<dim;a++) {
	imp(a) = 0.;
  }
  X.p = imp;

  int tracage = freq;
  ofstream energie("energie");
  energie.setf(ios::scientific);
  energie<<setprecision(10);
  ofstream position("position");

  double H_init = H(X.q,X.p);

  position << 0. << " ";
  for (int a=0;a<dim;a++) {
    position<< X.q(a) << " ";
  }
  for (int a=0;a<dim;a++) {
    position<< X.p(a) << " ";
  }
  position<<  endl;
  energie << 0. << "  " << H(X.q,X.p)-H_init <<endl;
    
  // integration
  for (int i = 0; i < pas; i++) {
//		X = SymplecticEuler(X);
		X = SymplecticEuler_proj(X,H_init);
		
	if (tracage == freq) {
	  tracage = 0;
	  position << (i+1)*dt << " ";
	  for (int a =0;a<dim;a++) {
	    position<< X.q(a) << "  ";
	  }
	  for (int a=0;a<dim;a++) {
	    position<< X.p(a)<<" ";
	  }
	  position<< endl;
	  energie << (i+1)*dt << "  " << H(X.q,X.p)-H_init << endl;
	}
	tracage += 1;
  }  

}
