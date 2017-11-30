cmake_minimum_required (VERSION 3.0)

project(TestNumericalAnalysis CXX)

if (MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11")
    set(TestNumericalAnalysisSRC src)
    set(EIGEN_INCLUDE_DIR D:/eigen3)

    find_package(Eigen3 REQUIRED)

    include_directories(${EIGEN_INCLUDE_DIR})

    add_executable(${PROJECT_NAME}
                   ${TestNumericalAnalysisSRC}/TestNumericalAnalysis.cxx
                   ${TestNumericalAnalysisSRC}/TestNewtonSolver.cxx
                   ${TestNumericalAnalysisSRC}/TestGradientDescentSolver.cxx
                   ${TestNumericalAnalysisSRC}/NumericalAnalysis.cxx
                   )


    target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
    
else()

endif()
