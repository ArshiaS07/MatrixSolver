#include "../EquationMatrix/EquationMatrix.hpp"
#include <vector>
#include <algorithm>

namespace ArshiaMath{

// Class for solving matrices
class MatrixSolver {
public:
    // Constructor
    MatrixSolver(EquationMatrix equationMatrix, double guess = 20.0, double tolerance = 1e-6, size_t maxIterations = 1e7);

    // Methods for solving the matrix
    bool solveGaussSeidel();
    bool solveTDMA();

    // Data members
    std::vector<double> m_x; // Solution vector
    bool checkTDMA(void); // Method for checking if TDMA method is applicable
protected:
private:
    const EquationMatrix m_equationMatrix; // Equation matrix
    const double m_tolerance; // Tolerance for convergence
    const size_t m_maxIterations; // Maximum number of iterations
    double m_guess; // Initial guess for the solution
};


// Class for solving matrices
class MatrixSolver {
public:
    // Constructor
    MatrixSolver(EquationMatrix equationMatrix, double guess = 20.0, double tolerance = 1e-6, size_t maxIterations = 1e7);

    // Methods for solving the matrix
    bool solveGaussSeidel();
    bool solveTDMA();

    // Data members
    std::vector<double> m_x; // Solution vector
    bool checkTDMA(void); // Method for checking if TDMA method is applicable
};

}