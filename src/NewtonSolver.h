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
#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

// STD
#include <iostream>
#include <stdlib.h>
#include <vector>

template<typename VectorType, typename Function, typename Jacobian>
class NewtonSolver
{
public:
  // Default constructor
  NewtonSolver();

  // Constructor with function and starting point specification
  NewtonSolver(VectorType argX0, Function argF, Jacobian argJ);

  // Set the starting point of the newton algorithm
  void SetX0(VectorType argX0);

  // Solve the non-linear equation
  void Solve();

  // return the solution point
  VectorType GetXf();

  // Set the epsilon stop condition
  void SetEpsilon(double argEspilon);

  // Set the max iteration stop condition
  void SetMaxIteration(unsigned int argMaxIt);

  // Set if the jacobian should be estimated
  void SetShouldEstimateJacobian(bool input);

  // Set the delta X to estimate Jacobian
  void SetdX(VectorType argdX);

  // Set verbal / not verbal mode
  void SetVerbal(bool input);

protected:
  // Starter point to initialize the
  // newton solver algorithm
  VectorType X0;

  // Solution point
  VectorType Xf;

  // Function to solve
  Function F;

  // Gradient of the function to solve
  Jacobian J;

  // If true, the jacobian is estimated
  // the estimation is in O(N^2) where N is
  // the dimension of the vector space
  bool ShouldEstimateJacobian;

  // delta X used to estimate the jacobian
  // if ShouldEstimateJacobian is set to true
  VectorType dX;

  // Stop condition
  double Epsilon;
  unsigned int MaxIteration;

  // Estimate the Jacobian at X
  auto EstimateJacobian(VectorType X);

  // compute the jacobian
  auto ComputeJacobian(VectorType X);

  // Verbal mode
  bool Verbal;
};

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
NewtonSolver<VectorType, Function, Jacobian>::NewtonSolver()
{
  // Startting point set at zero
  for (int k = 0; k < X0.rows(); ++k)
  {
    this->X0(k) = 0;
  }

  this->MaxIteration = 25;
  this->Epsilon = 1e-8;
  this->ShouldEstimateJacobian = false;
  this->Verbal = false;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
NewtonSolver<VectorType, Function, Jacobian>::NewtonSolver(VectorType argX0, Function F, Jacobian J)
{
  this->F = argF;
  this->J = argJ;
  this->X0 = argX0;
  this->MaxIteration = 25;
  this->Epsilon = 1e-8;
  this->ShouldEstimateJacobian = false;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetX0(VectorType argX0)
{
  this->X0 = argX0;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetdX(VectorType argdX)
{
  this->dX = argdX;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetShouldEstimateJacobian(bool input)
{
  this->ShouldEstimateJacobian = input;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
VectorType NewtonSolver<VectorType, Function, Jacobian>::GetXf()
{
  return this->Xf;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetEpsilon(double argEpsilon)
{
  this->Epsilon = argEpsilon;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetMaxIteration(unsigned int argMaxIt)
{
  this->MaxIteration = argMaxIt;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::SetVerbal(bool input)
{
  this->Verbal = input;
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
void NewtonSolver<VectorType, Function, Jacobian>::Solve()
{
  // Initialization of the algorithm parameters
  this->Xf = this->X0;
  unsigned int nbrIt = 0;

  // first computation of F(x)
  // and Jacobian(X)
  auto Y = this->F(X0);
  auto Jx = ComputeJacobian(X0);

  // Iterate while we are far from the solution
  while (Y.norm() > this->Epsilon && nbrIt < this->MaxIteration)
  {
    // compute one step
    this->Xf = this->Xf - Jx.inverse() * Y;

    // update jacobian, norme and value
    Y = this->F(this->Xf);
    Jx = ComputeJacobian(Xf);

    nbrIt++;
  }

  if (this->Verbal)
  {
    std::cout << "iteration made : " << nbrIt << std::endl;
    std::cout << "solution : " << this->Xf << std::endl;
    std::cout << "value : " << Y << std::endl;
    std::cout << "Norme : " << Y.norm() << std::endl;
  }
}

//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
auto NewtonSolver<VectorType, Function, Jacobian>::ComputeJacobian(VectorType X)
{
  auto J1 = this->EstimateJacobian(X);
  auto J2 = this->J(X);

  if (this->ShouldEstimateJacobian)
  {
    return this->EstimateJacobian(X);
  }

  return this->J(X);
}


//-----------------------------------------------------------
template<typename VectorType, typename Function, typename Jacobian>
auto NewtonSolver<VectorType, Function, Jacobian>::EstimateJacobian(VectorType X)
{
  auto estimatedJacobian = this->J(X);

  const unsigned int h = estimatedJacobian.rows();
  const unsigned int w = estimatedJacobian.cols();

  for (int j = 0; j < w; ++j)
  {
    // estimation of dF / dxj
    VectorType dxj;
    for (int i = 0; i < w; ++i)
    {
      dxj(i) = 0;
    }
    dxj(j) = this->dX(j);

    VectorType dFdXj = (this->F(X + dxj) - this->F(X - dxj)) / (2.0 * this->dX(j));
    for (int i = 0; i < h; ++i)
    {
      estimatedJacobian(i, j) = dFdXj(i);
    }
  }

  return estimatedJacobian;
}
#endif // NEWTON_SOLVER_H