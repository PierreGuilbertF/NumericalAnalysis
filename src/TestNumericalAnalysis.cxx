/*=========================================================================

Program:   Numerical Analysis tools
Module:    TestNumericalAnalysis.cxx

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

// EIGEN
#include <Eigen/Dense>

// STD
#include <iostream>

// LOCAL
#include "NumericalAnalysis.h"
#include "TestNewtonSolver.h"
#include "TestGradientDescentSolver.h"


int main()
{
  int nbrError = 0;
  // Test Newton Solver
  nbrError += TestNewtonSolverAnalytical();
  nbrError += TestNewtonSolverApproximated();

  // Test Gradient Descent Solver
  nbrError += TestGradientDescentAnalytical();

  std::cout << "Number of error in this test : " << nbrError << std::endl;
  system("PAUSE");
}