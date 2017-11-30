/*=========================================================================

Program:   Numerical Analysis tools
Module:    NewtonSolver.h

Copyright (c) Guilbert Pierre
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
* @brief   Numerical Analysis
*
*/
#ifndef GRADIENT_DESCENT_SOLVER_H
#define GRADIENT_DESCENT_SOLVER_H

// STD
#include <iostream>
#include <stdlib.h>
#include <vector>

template<typename VectorType, typename Function, typename Gradient>
class GradientDescentSolver
{
public:
  // default constructor
  GradientDescentSolver();

  // Constructor with function and starting point specification
  GradientDescentSolver(VectorType argX0, Function argF, Gradient argG);

  // set the initial point
  void SetX0(VectorType argX0);

  // Set the max iteration stop condition
  void SetMaxIteration(unsigned int input);

  // Set the Epsilon stop condition
  void SetEpsilon(double input);

  // Solve the minimisation problem
  void Solve();
protected:
  // Initial point
  VectorType X0;

  // solution
  VectorType Xf;

  // Function to minimize
  Function F;

  // Gradient of the function to minimize
  Gradient G;

  // step iteration
  double stepGrad;

  // Condition to stop the solver algorithm
  unsigned int MaxIteration;
  double Epsilon;
};

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
GradientDescentSolver<VectorType, Function, Gradient>::GradientDescentSolver()
{
  // Initialization of the initial point to zero
  for (int k = 0; k < X0.rows(); ++k)
  {
    X0(k) = 0;
  }

  this->MaxIteration = 100;
  this->Epsilon = 1e-8;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
GradientDescentSolver<VectorType, Function, Gradient>::GradientDescentSolver(VectorType argX0, Function argF, Gradient argG)
{
  this->X0 = argX0;
  this->F = argF;
  this->G = argG;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
void GradientDescentSolver<VectorType, Function, Gradient>::SetX0(VectorType argX0)
{
  this->X0 = argX0;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
void GradientDescentSolver<VectorType, Function, Gradient>::SetEpsilon(double input)
{
  this->Epsilon = input;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
void GradientDescentSolver<VectorType, Function, Gradient>::SetMaxIteration(unsigned int input)
{
  this->MaxIteration = input;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Gradient>
void GradientDescentSolver<VectorType, Function, Gradient>::Solve()
{
  unsigned int nbrIte = 0;
  this->Xf = this->X0;
  double lambda = 2.0 / (2.0 + nbrIte);
  // First computation of the gradient
  VectorType Gx = this->G(this->X0);

  while (Gx.norm() > this->Epsilon && nbrIte < this->MaxIteration)
  {
    this->Xf = this->Xf - lambda * Gx;

    // update variables
    Gx = this->G(this->Xf);
    lambda = 2.0 / (2.0 + nbrIte);
    nbrIte++;
  }

  std::cout << "iteration made : " << nbrIte << std::endl;
  std::cout << "solution : " << this->Xf << std::endl;
  std::cout << "value : " << this->F(this->Xf) << std::endl;
  std::cout << "Norme of gradient : " << Gx.norm() << std::endl;
}
#endif // GRADIENT_DESCENT_SOLVER_H