/*=========================================================================

Program:   Numerical Analysis tools
Module:    TestNewtonSolver.cxx

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
#include "TestNewtonSolver.h";

//-----------------------------------------------------------
struct F
{
  Eigen::Matrix<double, 2, 1> operator()(Eigen::Matrix<double, 2, 1> X)
  {
    Eigen::Matrix<double, 2, 1> Y;
    Y(0) = (1.0 * std::pow(X(0), 3) - 54.0 * std::pow(X(0), 2)
      + 609.0 * X(0) - 1960) * X(1);
    Y(1) = X(0) * std::pow(X(1), 2) - 8.0 * std::pow(X(1), 2) - 4 * X(0) * X(1)
      + 4 * X(0) + 32 * X(1) - 32;
    return Y;
  }
};

//-----------------------------------------------------------
struct JF
{
  Eigen::Matrix<double, 2, 2> operator()(Eigen::Matrix<double, 2, 1> X)
  {
    Eigen::Matrix<double, 2, 2> Jacobian;
    Jacobian(0, 0) = (3.0 * std::pow(X(0), 2) - 108 * X(0) + 609) * X(1);
    Jacobian(0, 1) = (1.0 * std::pow(X(0), 3) - 54.0 * std::pow(X(0), 2)
      + 609.0 * X(0) - 1960);
    Jacobian(1, 0) = std::pow(X(1), 2) - 4 * X(1) + 4;
    Jacobian(1, 1) = 2 * X(0) * X(1) - 16 * X(1) - 4 * X(0) + 32;
    return Jacobian;
  }
};

//-----------------------------------------------------------
int TestNewtonSolverAnalytical()
{
  std::cout << "-----TEST NEWTON SOLVER ANALYTICAL-----" << std::endl;

  NewtonSolver<Eigen::Matrix<double, 2, 1>, F, JF> Solver;

  // real solution is (7, 2)
  Solver.SetX0(Eigen::Matrix<double, 2, 1>(8.49, 3.65));
  Solver.SetMaxIteration(150);
  Solver.SetEpsilon(1e-16);
  Solver.SetVerbal(true);
  Solver.Solve();

  Eigen::Matrix<double, 2, 1> Solution = Solver.GetXf();

  double d = std::max(std::abs(Solution(0) - 7.0), std::abs(Solution(1) - 2.0));

  if (d > 1e-5)
  {
    std::cout << "TEST FAILED : " << d << std::endl;
    return 1;
  }

  std::cout << "TEST SUCCEEDED" << std::endl;
  return 0;
}

//-----------------------------------------------------------
int TestNewtonSolverApproximated()
{
  std::cout << "-----TEST NEWTON SOLVER APPROXIMATED-----" << std::endl;

  NewtonSolver<Eigen::Matrix<double, 2, 1>, F, JF> Solver;

  // real solution is (7, 2)
  Solver.SetX0(Eigen::Matrix<double, 2, 1>(8.49, 3.65));
  Solver.SetMaxIteration(150);
  Solver.SetEpsilon(1e-16);
  Solver.SetdX(Eigen::Matrix<double, 2, 1>(1e-3, 1e-3));
  Solver.SetShouldEstimateJacobian(true);
  Solver.SetVerbal(true);
  Solver.Solve();

  Eigen::Matrix<double, 2, 1> Solution = Solver.GetXf();

  double d = std::max(std::abs(Solution(0) - 7.0), std::abs(Solution(1) - 2.0));

  if (d > 1e-5)
  {
    std::cout << "TEST FAILED : " << d << std::endl;
    return 1;
  }

  std::cout << "TEST SUCCEEDED" << std::endl;
  return 0;
}