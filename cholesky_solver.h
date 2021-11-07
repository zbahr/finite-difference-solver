/*! \file
 *  Cholesky_Solver class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef CHOLESKY_SOLVER_H
#define CHOLESKY_SOLVER_H
#include "solver_strategy.h"
#include "l_triangle_matrix.h"
#include "u_triangle_matrix.h"

/*! Chomsky solver class */
template <class T>
class Cholesky_Solver : public virtual Solver_Strategy<T> { 
    private:
        int m_size;

    public:
        /*! Function operator overload implementing Cholesky decomposition followed by substitution
          *
          * \param matrix the matrix to perform Cholesky decomposition
          * \param vec the solution vector that is paired with `matrix` in the Cholesky decomposition process
          * \return a vector x representing the solution of matrix * x = vec
          * 
          * \pre pre-conditions for auxiliary functions should be met
          * \pre parameter matrix.get_size() == parameter vector.get_size()
          * \pre T(0) is defined
          * \pre dividend values in the algorithm should be non-zero; values in square root should be >= 0
          * \post (see return)
          * \throws domain_error thrown if preconditions 2 or 4 are broken
        */
        Vector<T> solve(const Base_Matrix<T>& matrix, Vector<T> vec);
};

#include "cholesky_solver.hpp"
#endif