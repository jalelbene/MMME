//--------------------------------------------
//
//           Un fluide de particules interagissant via LJ
//
//   compilation : make
//
//--------------------------------------------

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "matrice.hpp"
#include <Imagine/Graphics.h>

using namespace std;
using namespace Imagine;


// pas de temps
const double dt = 0.001;

// nb de pas de temps
const int pas = 1000;

// frequence d'affichage
const int freq = 1;

// nb de particules par dimension
const int NbPart_dim = 6;
// nb total de particules
const int NbPart = 36; // doit etre egal a NbPart_dim a la puissance dim
// dimension
const int dim = 2;

struct Particule
{
  mat q;       
  mat p;
};


//------------------------
//  forces et potentiels
//------------------------

double V (Particule X)
{

  double pot = 0.;
  double dis2 = 0.;

  for (int i=0; i<(NbPart-1); i++) {
	for (int j=(i+1); j<NbPart; j++){
	  dis2 = 0.;
	  for (int a=0; a<dim; a++) {
		dis2 += pow(X.q(i,a)-X.q(j,a),2);
	  }
	  pot += 1./pow(dis2,6) - 2./pow(dis2,3);
	}
  }

  for (int i=0; i<NbPart; i++) {
	dis2 = 0.;
	for (int a=0; a<dim; a++) {
	  dis2 += pow(X.q(i,a),2);
	}
	pot += 0.5*dis2;
  }
  
  return pot;
}

mat f (Particule X)
{  
  mat force(NbPart,dim);
  force.zeros();
  double dis2 = 0.;
  double dis_[dim];

  for (int i=0; i<(NbPart-1); i++) {
	for (int j=(i+1); j<NbPart; j++){
        dis2 = 0.;
        for (int a=0; a<dim; a++) {
            dis_[a] = X.q(i,a)-X.q(j,a);
            dis2 += pow(dis_[a],2);
        }
        for (int a=0; a<dim; a++) {
            force(i,a) += 12*pow(dis2,5)*(dis_[a])/pow(dis2,12) - 12*pow(dis2,2)*(dis_[a])/pow(dis2,6);
        }
	  // a completer
	}
  }

  for (int i=0; i<NbPart; i++) {
      for(int a=0; a<dim; a++) {
          force(i,a) -= X.q(i,a);
      }
	// a completer
  }
    
  return force;
}

mat hessien_V (Particule X)
{
  mat hess(dim*NbPart,dim*NbPart);
  hess.zeros();
  int index_1,index_2,index_3,index_4;
  double over_d16,over_d14,over_d10,over_d8;
  double pref1,pref2;
  double dis2;
  double dis_[dim];
  
  for (int i=0; i < NbPart; i++) {
	for (int a=0;a<dim;a++) {
	  index_1 = a+i*dim;
	  for (int j=0; j < NbPart; j++) {
		if (i != j) {
		  index_4 = a + j*dim;
		  for (int b=0; b< dim; b++) {
			index_2 = b+j*dim;
			index_3 = b + i*dim;
			
			dis2 = 0.;
			for (int c=0; c<dim; c++) {
			  dis_[c] = X.q(i,c)-X.q(j,c);
			  dis2 += pow(dis_[c],2);
			}
			over_d16 = 1./pow(dis2,8);
			over_d10 = 1./pow(dis2,5);
			over_d14 = 1./pow(dis2,7);
			over_d8 = 1./pow(dis2,4);
			pref1 = 12.*14.*over_d16 - 12.*8.*over_d10;
			pref2 = -12.*over_d14 + 12.*over_d8;
			// on complete (i,i)
			hess(index_1,index_3) += pref1*dis_[a]*dis_[b];
            
			// on complete (i,j)
			hess(index_1,index_2) -= pref1*dis_[a]*dis_[b];
			
		  }
		  hess(index_1,index_1) += pref2;
		  hess(index_1,index_4) -= pref2;
		}
	  }
	}
  }

  for (int i=0; i < NbPart; i++) {
	for (int a=0;a<dim;a++) {
	  index_1 = a+i*dim;
	  hess(index_1,index_1) += 1.;
	}
  }

  return hess;
}

// distance minimale
double dist_min (Particule X)
{  
  double dis = 0.;
  double dis_[dim];
  double min = 100.;

  for (int i=0; i<(NbPart-1); i++) {
	for (int j=(i+1); j<NbPart; j++){
	  dis = 0.;
	  for (int a=0; a<dim; a++) {
		dis_[a] = X.q(i,a)-X.q(j,a);
		dis += pow(dis_[a],2);
	  }
	  dis = sqrt(dis);
	  if (dis < min) {
		min = dis;
	  }
	}
  }
  
  return min;
}


double H (Particule X)
{
  double pot = V(X);
  for (int i=0; i< NbPart; i++) {
	for (int a=0; a<dim; a++) {
	  pot += 0.5*pow(X.p(i,a),2);
	}
  }
  return pot;
}

// correction d'ordre 2 du hamiltonien modifie
double corr_H2 (Particule X)
{
  double res = 0.;

  mat hess = hessien_V(X);
  mat force = f(X);
  int index_1, index_2;

  for (int i=0; i < NbPart; i++) {
	for (int a=0;a<dim;a++) {
	  index_1 = a+i*dim;
	  for (int j=0; j < NbPart; j++) {
		for (int b=0; b< dim; b++) {
		  index_2 = b+j*dim;
		  res += 2.*X.p(i,a)*X.p(j,b)*hess(index_1,index_2);
		}
	  }
	  res -= force(i,a)*force(i,a);
	}
  }
  
  return res;

}

// hamiltonien modifie d'ordre 2
double H2 (Particule X, double dt)
{
  double res = H(X) + dt*dt*corr_H2(X)/24.;
  return res;

}


//--------------------------
//  algorithme Verlet
//-------------------------

Particule Verlet(Particule X)
{
  Particule Y;
  Y.p = X.p + f(X)*(0.5*dt);
  Y.q = X.q + Y.p*dt;
  Y.p = Y.p + f(Y)*(0.5*dt);

  return Y;
}

//----------------------------
// Choix des conditions initiales
//---------------------------

void choix_ci(mat & position, mat & impulsion)
{

  int index;
  
  // les positions
  for (int a=0; a<NbPart_dim; a++) {
	for (int b=0; b<NbPart_dim; b++) {
	  index = b + a*NbPart_dim;
	  position(index,0) = b;
	  position(index,1) = a;
	}
  }

  double alea = 0.;
  double pi = 4.*atan(1.);
  double ampli = 5.;
  // les vitesses
  for (int a=0; a<NbPart_dim; a++) {
	for (int b=0; b<NbPart_dim; b++) {
	  alea = 2.*pi*rand()/RAND_MAX;
	  index = b + a*NbPart_dim;
	  impulsion(index,0) = ampli*cos(alea);
	  impulsion(index,1) = ampli*sin(alea);
	}
  }
 
}


//----------------------------
// integration des equations
//---------------------------

int main () 
{
  //  srand(time(NULL));
    
  Window w_H = openWindow(1000,800);
  Window w_H2 = openWindow(1000,800);

  Particule X;
  mat position0(NbPart,dim);
  mat impulsion0(NbPart,dim);
  choix_ci(position0,impulsion0);
  X.q = position0;
  X.p = impulsion0;

  int tracage = freq;
  ofstream energie("energie");
//  ofstream position("position");
    ofstream amplitude("amplitude");

  double H_init = H(X);
  double H2_init = H2(X,dt);
  
  double H_min = H_init;
  double H_max = H_init;
  double H2_min = H2_init;
  double H2_max = H2_init;
//  position << 0. << " "<< X.q << " " << X.p << endl;
  energie << 0. << " " << H(X)-H_init << " "<< H2(X,dt)-H2_init << endl;

    
  Color c1 = RED;
  Color c2 = BLUE;
    
  // integration
  for (int i = 0; i < pas; i++) {
    Particule Y;
    Y.q = X.q;
    Y.p = X.p;
      
	X = Verlet(X);
      
    if (H(X) < H_min) {
        H_min = H(X);
    }

    if (H(X) > H_max) {
        H_max = H(X);
    }
      
    if (H2(X,dt) < H2_min) {
        H2_min = H2(X,dt);
    }
      
    if (H2(X,dt) > H_max) {
        H2_max = H2(X,dt);
    }
      
    setActiveWindow(w_H);
    drawLine(i*dt*1000, H(Y)*0.03, (i+1)*dt*1000, H(X)*0.03, c1);
//    milliSleep(10);
      
    setActiveWindow(w_H2);
    drawLine(i*dt*1000, H2(Y,dt)*0.03, (i+1)*dt*1000, H2(X,dt)*0.03, c2);
//    milliSleep(10);
		
	if (tracage == freq) {
	  tracage = 0;
//	  position << (i+1)*dt << " "<< X.q << " " << X.p << endl;
	  energie << (i+1)*dt << " " << H(X)-H_init << " "<< H2(X,dt)-H2_init << " "<<dist_min(X)<<endl;
	}
	tracage += 1;
  }
  
  
  amplitude << " H_min = " << H_min
            << "\n H_max = " << H_max
            << "\n\n H2_min = " << H2_min
            << "\n H2_max = " << H2_max;
    
//  endGraphics();
    
  return 0;
}
