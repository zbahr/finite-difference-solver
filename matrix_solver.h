/*! \file
 *  Matrix_Solver class definition/declaration.
 */

#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H
#include "gaussian_solver.h"
#include "cholesky_solver.h"
#include "generators.hpp"

/*! Matrix solver class */
template <class T>
class Matrix_Solver { 
    private:
        // Matrix-vector containers
        int m_size;
        Base_Matrix<T>* m_matrix;
        Vector<T> m_vec;

        // Solving strategy 
        Solver_Strategy<T>* m_method;

    public:
        /*! Constructor for a given matrix-vector pair
          *
          * \param matrix the matrix to solve
          * \param vec the solution vector that is paired with `matrix`
          * 
          * \pre none
          * \post solver class is constructed with given matrix-vector pair
        */
        Matrix_Solver(const Base_Matrix<T>& matrix, const Vector<T>& vec);

        /*! Constructor to generate a matrix-vector pair via the finite difference
          * method
          *
          * \param lower_bound the lower bound of the mesh
          * \param upper_bound the upper bound of the mesh
          * \param mesh_length the mesh length
          * \param gauss_override flag to force gaussian elimination strategy
          * \param upper a pointer to the upper boundary function
          * \param lower a pointer to the lower boundary function
          * \param right a pointer to the right boundary function
          * \param left a pointer to the left boundary function
          * 
          * \pre upper_bound > lower_bound
          * \pre mesh_length > 0
          * \post solver class is constructed with matrix-vector pair utilizing
          *       finite difference method
          * \throws domain_error thrown if pre-conditions broken
        */
        Matrix_Solver(const double& lower_bound, const double& upper_bound, const int& mesh_length, const bool& gauss_override, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double));

        /* Destructor */
        ~Matrix_Solver() { if (m_matrix != nullptr) delete m_matrix; }

        /*! Driver function to select solver strategy and begin solving process
          * 
          * \return a vector containing the solution to the matrix-vector members
          * 
          * \pre pre-conditions for auxiliary functions should be met
          * \pre parameter matrix.get_size() == parameter vector.get_size()
          * \post (see return)
          * \throws domain_error thrown if size of parameters are not equal
        */
        Vector<T> solve();
};

#include "matrix_solver.hpp"
#endif