#include <iostream>
#include "MatrixSolver.hpp" // Assuming the MatrixSolver class is defined in MatrixSolver.h

int main() {
    // Example usage of MatrixSolver

    // Create a coefficient matrix A
    ArshiaMath::Matrix A = {
        {4, 1, 0},
        {1, 4, 1},
        {0, 1, 4}
    };

    // Create a right-hand side vector b
    std::vector<double> b = {5, 6, 5};

    // Create a MatrixSolver object
    MatrixSolver solver(A, b);

    // Solve the system of equations using Gauss-Seidel method
    bool success = solver.solveGaussSeidel();

    if (success) {
        std::cout << "Solution (Gauss-Seidel method): ";
        std::vector<double> solution = solver.getSolution();
        for (double x : solution) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Failed to solve the system using Gauss-Seidel method." << std::endl;
    }

    // Solve the system of equations using TDMA (Thomas algorithm) if it satisfies the requirements
    if (solver.checkTDMA()) {
        success = solver.solveTDMA();

        if (success) {
            std::cout << "Solution (TDMA): ";
            std::vector<double> solution = solver.getSolution();
            for (double x : solution) {
                std::cout << x << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Failed to solve the system using TDMA." << std::endl;
        }
    } else {
        std::cout << "The coefficient matrix does not satisfy the requirements for TDMA." << std::endl;
    }

    return 0;
}