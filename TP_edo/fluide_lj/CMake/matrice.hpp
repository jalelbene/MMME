#ifndef MATRICE_HPP
#define MATRICE_HPP
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <math.h>
#include <string>


///////////////////////////////////////////////
//// CLASSE vecteur

class vec
{
public:

  int dim;
  double* coord;

  vec();
  // constructeur par defaut

  vec(int);
  //construit un vec de taille donnée

  vec(const vec &);
  // operateur de recopie

  ~vec();
  // destruction

  vec & operator= (const vec &);
  // affectation

  vec operator+ (const vec &);
  // addition

  vec operator- (const vec &);
  // soustraction

  vec operator* (const double a);
  // multiplication par un scalaire

  double & operator() (int);
  // renvoie la composante d'un vec

  double & operator[] (int);
  // idem que ()
  
  void set_size(int);
  // met le vec courant à la longueur souhaitée

  void zeros();
  // met un vecteur a 0

  //ostream& affiche (ostream& os) const;
  // pour afficher un vecteur
  
};


/////////////////////////////////////////////
/////////// CLASSE matrice

class mat
{

public:

  int nligne;
  int ncol;
  vec * col;

  mat();
  // constructeur par defaut

  mat(int,int);
  // constructeur avec choix de la dimension
	
  mat(const mat &);
  // operateur de recopie

  ~mat();
  // destruction

  mat & operator= (const mat &);
  // affectation

  mat operator+ (const mat &);
  // addition
  
  mat operator* (const double a);
  // multiplication par un scalaire

  double & operator() (int, int) const;
  // renvoie une composante de mat

  void set_size(int,int);
  // met la mat courant à la longueur souhaitée

  void zeros();
  // met une mat a 0

  vec get_col(int a);
  // renvoie une colonne

  void set_col(int a, vec v);
  // met a jour une colonne
  
  //ostream& affiche(ostream& os) const;
  // pour afficher une matrice
  
};

//ostream& operator <<(ostream& os, const mat & d);
// pour afficher une matrice

//ostream& operator <<(ostream& os, const vec & d);
// pour afficher un vecteur

#endif

