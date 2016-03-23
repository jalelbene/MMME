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
#include <Imagine/Graphics.h>

using namespace std;
using namespace Imagine;


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


int random(const int &a){
    return rand()%a;
}


//----------------------------
// integration des equations
//---------------------------

int main () 
{
  
  Window w_q = openWindow(700,1000);
  Window w_H = openWindow(700,1000);

  Particule X,Z;
    
  // choix des CI
  double pos=0.;
  X.q = pos;
  Z.q = pos;

  double imp = 1.8;
  X.p = imp;
  Z.p = imp;

  int tracage = freq;
  ofstream energie("energie");
  ofstream position("position");
  ofstream amplitude("amplitude");

  double H_init = H(X.q,X.p);
	
  position << 0. << " "<< X.q << " " << X.p << endl;
  energie << 0. << "  " << H(X.q,X.p)-H_init <<endl;

  double q_min_V = X.q;
  double q_max_V = X.q;
  double H_min_V = H_init;
  double H_max_V = H_init;
    
  double q_min_L = X.q;
  double q_max_L = X.q;
  double H_min_L = H_init;
  double H_max_L = H_init;
  
  Color c1 = RED;
  Color c2 = BLUE;

  // integration
  for (int i = 0; i < pas; i++) {
    
    Particule Y,T;
    Y.q = X.q;
    Y.p = X.p;
    
    T.q = Z.q;
    T.p = Z.p;

    // Recalcul de X et Z
    Z = Verlet(Z);
	X = Lobatto(X);
      
    // Mise à jour des valeurs extremales pour la méthode Verlet (Z)
    if( Z.q < q_min_V ){
        q_min_V =Z.q;
    }
    
    if( Z.q > q_max_V ){
        q_max_V = Z.q;
    }
      
    if( H(Z.q,Z.p) < H_min_V ){
        H_min_V = H(Z.q,Z.p);
    }
    
    if( H(Z.q,Z.p) > H_max_V ){
        H_max_V = H(Z.q,Z.p);
    }
	
      
    // Mise à jour des valeurs extremales pour la méthode Lobatto (X)
    if( X.q < q_min_L ){
        q_min_L =X.q;
    }
      
    if( X.q > q_max_L ){
        q_max_L = X.q;
    }
      
    if( H(X.q,X.p) < H_min_L ){
        H_min_L = H(X.q,X.p);
    }
      
    if( H(X.q,X.p) > H_max_L ){
        H_max_L = H(X.q,X.p);
    }
      
//    Color c = Color(random(256),random(256),random(256));
    setActiveWindow(w_q);
    drawLine(i*dt*5, Y.q*100+350, (i+1)*dt*5, X.q*100+350, c1);
    drawLine(i*dt*5, T.q*100+350, (i+1)*dt*5, Z.q*100+350, c2);
    milliSleep(10);
      
    setActiveWindow(w_H);
    drawLine(i*dt*5, H(T.q,T.p)*10000-6000, (i+1)*dt*5, H(Z.q,Z.p)*10000-6000, c1);
    drawLine(i*dt*5, H(Y.q,Y.p)*10000-6000, (i+1)*dt*5, H(X.q,X.p)*10000-6000, c2);
//    milliSleep(10);
	
	if (tracage == freq) {
	  tracage = 0;
	  position << (i+1)*dt << " "<< X.q << "  " << X.p << endl;
	  energie << (i+1)*dt << "  " << H(X.q,X.p)-H_init <<endl;
	}
	tracage += 1;
  }  

    
    amplitude << " position initiale = " << pos
              << "\n impulsion initiale = " << imp
              << endl;
    
    amplitude << "\n time step = " << dt
              << "\n number of time steps = " << pas
              << endl;
    
    amplitude << "\n t_min = 0"
              << "\n t_max = " << pas*dt
              << endl;
    
    amplitude << "\n q_min_V = " << q_min_V
              << "\n q_max_V = " << q_max_V
              << endl;
    
    amplitude << "\n H_min_V = " << H_min_V
              << "\n H_max_V = " << H_max_V
              << endl;
    
    amplitude << "\n q_min_L = " << q_min_L
              << "\n q_max_L = " << q_max_L
              << endl;
    
    amplitude << "\n H_min_L = " << H_min_L
              << "\n H_max_L = " << H_max_L
              << endl;
    
    endGraphics();
    
    return 0;
}
