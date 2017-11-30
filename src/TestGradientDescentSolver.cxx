/*=========================================================================

Program:   Numerical Analysis tools
Module:    TestGradientDescentSolver.cxx

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
#include "TestGradientDescentSolver.h"

struct F
{
  double operator()(Eigen::Matrix<double, 3, 1> X)
  {
    double y = std::pow(X(0) - 3.0, 2) + std::pow(X(1) - 6.4, 2) + std::pow(X(2) + 4.28, 2);
    return y;
  }
};

struct G
{
  Eigen::Matrix<double, 3, 1> operator()(Eigen::Matrix<double, 3, 1> X)
  {
    Eigen::Matrix<double, 3, 1> Gradient;
    Gradient(0) = 2.0 * std::pow(X(0) - 3.0, 1);
    Gradient(1) = 2.0 * std::pow(X(1) - 6.4, 1);
    Gradient(2) = 2.0 * std::pow(X(2) + 4.28, 1);
    return Gradient;
  }
};

int TestGradientDescentAnalytical()
{
  GradientDescentSolver<Eigen::Matrix<double, 3, 1>, F, G> Solver;
  Solver.SetEpsilon(1e-8);
  Solver.SetMaxIteration(150);
  Solver.Solve();
  return 0;
}