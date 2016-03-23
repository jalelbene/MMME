// Etude d'une particule 1D soumise a un potentiel de Lennard-Jones 

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
const int pas = 20000;

// frequence d'affichage
const int freq = 1;

// position du mur de droite (le mur de gauche est en 0). 
const double mur = 10.;

struct Particule
{
  double q;
  double p;
};


//------------------------
//  potentiels, forces et derivees d'ordre superieur
//------------------------

// potentiel de Lennard-Jones elementaire
double W (double pos)
{
  double pot = 1./pow(pos,12) - 2./pow(pos,6);
  return pot;
}

// energie potentielle
double V (double pos)
{
  double pot = W(pos) + W(mur-pos);
  return pot;
}

// force associee au potentiel de Lennard-Jones elementaire
double f_W (double pos)
{
  double force = 12./pow(pos,13) - 12./pow(pos,7);
  return force;
}

// force associee a l'energie potentielle
double f (double pos)
{
  double force = f_W(pos)- f_W(mur-pos);
  return force;
}

// hessien associe au potentiel de Lennard-Jones elementaire
double hessien_W (double pos)
{
  double hess = 12.*13./pow(pos,14) - 12.*7./pow(pos,8);
  return hess;
}

// hessien associe a l'energie potentielle
double hessien_V (double pos)
{
  double hess = hessien_W(pos) + hessien_W(mur-pos);
  return hess;
}

// derivee 3ieme associee au potentiel de Lennard-Jones elementaire
double d3W (double pos)
{
  double res = -12.*13.*14./pow(pos,15) + 12.*7.*8./pow(pos,9);
  return res;
}

// derivee 3ieme associee a l'energie potentielle
double d3V (double pos)
{
  double res = d3W(pos) - d3W(mur-pos);
  return res;
}

// derivee 4ieme associee au potentiel de Lennard-Jones elementaire
double d4W (double pos)
{
//  double res = 0.;

  double res = 12.*13.*14.*15./pow(pos,16) - 12.*7.*8.*9./pow(pos,10);
  
  return res;
}

// derivee 4ieme associee a l'energie potentielle
double d4V (double pos)
{
  double res = d4W(pos) + d4W(mur-pos);
  return res;
}

// hamiltonien
double H (double pos, double imp)
{
  double pot = V(pos);
  pot += 0.5*pow(imp,2);
  return pot;
}

// correction d'ordre 2 pour le hamiltonien modifie
double corr_H2 (double pos, double imp)
{
  double res = 0.;

  double hess = hessien_V(pos);
  double force = f(pos);

  res = 2.*imp*imp*hess;
  res -= force*force;
  
  return res;
}

// correction d'ordre 4 pour le hamiltonien modifie
double corr_H4 (double pos, double imp)
{
  double res = 0.;

  double hess = hessien_V(pos);
  double force = f(pos);

  res += -pow(imp,4)*d4V(pos) - 3*hess*pow(force,2) + 12*pow(imp,2)*pow(hess,2) + 6*pow(imp,2)*d3V(pos)*force;
  
  return res;

}

// hamiltonien modifie d'ordre 2
double H2 (double pos, double imp, double dt)
{
  double res = H(pos,imp) + dt*dt*corr_H2(pos,imp)/24.;
  return res;

}

// hamiltonien modifie d'ordre 4
double H4 (double pos, double imp, double dt)
{
  double res = 0.;
  
  res = H(pos,imp) + dt*dt*corr_H2(pos,imp)/24. + pow(dt,4)/720.*corr_H4(pos,imp);
  
  return res;

}

//--------------------------
//  algorithme Explicte Euler
//-------------------------

Particule EE(Particule X)
{
  Particule Y;
  Y.p = X.p + f(X.q)*dt;
  Y.q = X.q + X.p*dt;
  return Y;
}

//--------------------------
//  algorithme Symplectic Euler
//-------------------------

Particule SE(Particule X)
{
  Particule Y;
  Y.p = X.p + f(X.q)*dt;
  Y.q = X.q + Y.p*dt;
  return Y;
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

//----------------------------
// integration des equations
//---------------------------

int main () 
{
  Particule X,Y,Z;
    
  double delta_H_V, delta_H2_V, delta_H4_V;
  double delta_H_EE, delta_H2_EE, delta_H4_EE;
  double delta_H_SE, delta_H2_SE, delta_H4_SE;

  vec error_max_V(3);
  vec error_max_EE(3);
  vec error_max_SE(3);

  // choix des CI (a choisir entre 0 et mur)
  double pos=5.;
  X.q = pos;
  Y.q = pos;
  Z.q = pos;
    
  double imp = -1;
  X.p = imp;
  Y.p = imp;
  Z.p = imp;

  int tracage = freq;
  ofstream energie("energie");
  energie.setf(ios::scientific);
  energie<<setprecision(10);
  ofstream position("position");
    
  ofstream results("results");
  results.setf(ios::scientific);
  results<<setprecision(10);

  double H_init = H(X.q,X.p);
  double H2_init = H2(X.q,X.p,dt);
  double H4_init = H4(X.q,X.p,dt);
    
  error_max_V(0) = 0;
  error_max_V(1) = 0;
  error_max_V(2) = 0;

  error_max_EE = error_max_V;
  error_max_SE = error_max_V;
    
	
  position << 0. << " "<< X.q << " " << X.p << endl;
  energie << 0. << "  " << H(X.q,X.p)-H_init << " " << H2(X.q,X.p,dt)-H2_init << " "<<H4(X.q,X.p,dt) - H4_init<<endl;

  
  // integration
  for (int i = 0; i < pas; i++) {
	X = Verlet(X);
    Y = EE(Y);
	Z = SE(Z);
      
    delta_H_V = H(X.q,X.p)-H_init;
    delta_H2_V = H2(X.q,X.p,dt)-H2_init;
    delta_H4_V = H4(X.q,X.p,dt)-H4_init;
    
    delta_H_EE = H(Y.q,Y.p)-H_init;
    delta_H2_EE = H2(Y.q,Y.p,dt)-H2_init;
    delta_H4_EE = H4(Y.q,Y.p,dt)-H4_init;
      
    delta_H_SE = H(Z.q,Z.p)-H_init;
    delta_H2_SE = H2(Z.q,Z.p,dt)-H2_init;
    delta_H4_SE = H4(Z.q,Z.p,dt)-H4_init;
      
    if(fabs(error_max_V(0)) < fabs(delta_H_V)){
        error_max_V(0) = delta_H_V;
    }

    if(fabs(error_max_V(1)) < fabs(delta_H2_V)){
        error_max_V(1) = delta_H2_V;
    }
    
    if(fabs(error_max_V(2)) < fabs(delta_H4_V)){
        error_max_V(2) = delta_H4_V;
    }
      
    if(fabs(error_max_EE(0)) < fabs(delta_H_EE)){
        error_max_EE(0) = delta_H_EE;
    }
      
    if(fabs(error_max_EE(1)) < fabs(delta_H2_EE)){
        error_max_EE(1) = delta_H2_EE;
    }
    
    if(fabs(error_max_EE(2)) < fabs(delta_H4_EE)){
        error_max_EE(2) = delta_H4_EE;
    }
      
    if(fabs(error_max_SE(0)) < fabs(delta_H_SE)){
        error_max_SE(0) = delta_H_SE;
    }
    
    if(fabs(error_max_SE(1)) < fabs(delta_H2_SE)){
        error_max_SE(1) = delta_H2_SE;
    }
    
    if(fabs(error_max_SE(2)) < fabs(delta_H4_SE)){
        error_max_SE(2) = delta_H4_SE;
    }
      
	if (tracage == freq) {
	  tracage = 0;
	  position << (i+1)*dt << " "<< X.q << "  " << X.p << endl;
	  energie << (i+1)*dt << "  " << H(X.q,X.p)-H_init << " " << H2(X.q,X.p,dt)-H2_init << " "<<H4(X.q,X.p,dt) - H4_init<<endl;
	}
	tracage += 1;
  }  

    
  results << " time step = " << dt
          << "\n number of time steps = " << pas
          << endl;
    
  results << "\n delta_H_EE = " << error_max_EE(0)
          << "\n delta_H2_EE = " << error_max_EE(1)
          << "\n delta_H4_EE = " << error_max_EE(2)
          << endl;
    
  results << "\n delta_H_SE = " << error_max_SE(0)
          << "\n delta_H2_SE = " << error_max_SE(1)
          << "\n delta_H4_SE = " << error_max_SE(2)
          << endl;
    
  results << "\n delta_H_V = " << error_max_V(0)
          << "\n delta_H2_V = " << error_max_V(1)
          << "\n delta_H4_V = " << error_max_V(2)
          << endl;

    
    return 0;
}
