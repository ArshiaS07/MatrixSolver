#include "MatrixSolver.hpp"

MatrixSolver::MatrixSolver(EquationMatrix equationMatrix, double guess, double tolerance, size_t maxIterations) :
	m_equationMatrix(equationMatrix),
	m_tolerance(tolerance),
	m_maxIterations(maxIterations),
	m_guess(guess)
{
	m_x = std::vector<double>(m_equationMatrix.getDimension(), m_guess);
};


bool MatrixSolver::solveGaussSeidel() {
	const size_t n = m_equationMatrix.getDimension();
	size_t k = 0;
	double InfNorm = std::numeric_limits<double>::max()
		, AbsNorm = std::numeric_limits<double>::max()
		, EucNorm = std::numeric_limits<double>::max()
		, Sum = std::numeric_limits<double>::max();

	std::vector<double> oldX(n);
	std::vector<double> sumVec(n - 1);
	while (InfNorm > m_tolerance && k < m_maxIterations)
	{

		for (size_t i = 0; i < n; i++) {
			Sum = 0;
			for (size_t j = 0; j < i; j++)
				Sum += m_equationMatrix.A[i][j] * m_x[j];
			for (size_t j = i + 1; j < n; j++)
				Sum += m_equationMatrix.A[i][j] * m_x[j];

			oldX[i] = m_x[i];
			m_x[i] = (m_equationMatrix.b[i] - Sum) / m_equationMatrix.A[i][i];
		}

		InfNorm = 0;
		AbsNorm = 0;
		EucNorm = 0;

		for (size_t i = 0; i < n; i++) {
			AbsNorm += abs(oldX[i] - m_x[i]);
			EucNorm += pow((oldX[i] - m_x[i]), 2);
			InfNorm = std::max(InfNorm, abs(oldX[i] - m_x[i]));
		}
		EucNorm = pow(EucNorm, 1 / 2);
		k++;
	}
	if (k = m_maxIterations)
		return false;
	return true;

}

bool MatrixSolver::checkTDMA(void) {
	bool isTDMA = true;
	for (size_t j = 2; j < m_equationMatrix.A.getCols(); j++) {
		isTDMA = isTDMA
			&& m_equationMatrix.A[0][j] == 0;
	}
	for (size_t i = 1; i < m_equationMatrix.A.getRows() - 1; i++) {
		for (size_t j = 0; j < m_equationMatrix.A.getCols(); j++) {
			isTDMA = isTDMA
				&& (((j == i - 1) || (j == i) || (j == i + 1)) ?
					true : m_equationMatrix.A[i][j] == 0);
		}
	}
	for (size_t j = 0; j < m_equationMatrix.A.getCols() - 1; j++) {
		isTDMA = isTDMA
			&& m_equationMatrix.A[m_equationMatrix.A.getRows()-1][j] == 0;
	}

	return isTDMA;
}

bool MatrixSolver::solveTDMA() {
	size_t n = m_equationMatrix.getDimension();
	ArshiaMath::Matrix Mat_A = m_equationMatrix.A;
	std::vector<double> Vec_b = m_equationMatrix.b;
	std::vector<double> Vec_x(n, 0);


	Mat_A[0][1] = Mat_A[0][1] / Mat_A[0][0];
	Vec_b[0] /= Mat_A[0][0];
	Mat_A[0][0] = 1;

	for (size_t i = 1; i < n - 1; i++) {
		//j   = i-1   i   i+1    i
		//i-1 =  1    c    0 ,   d
		//
		//j   = i-1   i   i+1    i
		//i   =  a    b    c ,   d
		//   b    = b -     a          *    c
		Mat_A[i][i] -= Mat_A[i][i - 1] * Mat_A[i - 1][i];
		//   d    = d -     a          *    d
		Vec_b[i] -= Mat_A[i][i - 1] * Vec_b[i - 1];
		//   a    = a -     a          *    1
		Mat_A[i][i - 1] = 0;
		//j   = i-1   i   i+1    i
		//i   =  0    b'   c',   d'
		Mat_A[i][i + 1] /= Mat_A[i][i];
		Vec_b[i] /= Mat_A[i][i];
		Mat_A[i][i] = 1;
		//j   = i-1   i   i+1    i
		//i   =  0    1    c ,   d
	}
	//last line
	Mat_A[n - 1][n - 1] -= Mat_A[n - 1][n - 2] * Mat_A[n - 2][n - 1];
	Vec_b[n - 1] -= Mat_A[n - 1][n - 2] * Vec_b[n - 2];
	Mat_A[n - 1][n - 2] = 0;
	Vec_b[n - 1] /= Mat_A[n - 1][n - 1];
	Mat_A[n - 1][n - 1] = 1;
	//Back-Substitution
	Vec_x[n - 1] = Vec_b[n - 1];
	for (size_t i = n - 2; i != -1; i--) {
		Vec_x[i] = Vec_b[i] - Mat_A[i][i + 1] * Vec_x[i + 1];
	}
	m_x = Vec_x;
	return true;

}