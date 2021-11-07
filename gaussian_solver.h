/*! \file
 *  Gaussian_Solver class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef GAUSSIAN_SOLVER_H
#define GAUSSIAN_SOLVER_H
#include "solver_strategy.h"

/*! Gaussian solver class */
template <class T>
class Gaussian_Solver : public virtual Solver_Strategy<T> { 
    private:
        int m_size;
        Vector<Vector<T>> m_matrix_data;
        Vector<T> m_scales;
        Vector<T> m_ratios;

    public:
        /*! Function operator overload implementing Gaussian elimination with scaled partial pivoting
          *
          * \param matrix the matrix to perform gaussian elimination on
          * \param vec the solution vector that is paired with `matrix` in the gaussian elimination process
          * \return a vector x representing the solution of matrix * x = vec
          * 
          * \pre pre-conditions for auxiliary functions should be met
          * \pre parameter matrix.get_size() == parameter vector.get_size()
          * \post (see return)
          * \throws domain_error thrown if size of parameters are not equal
        */
        virtual Vector<T> solve(const Base_Matrix<T>& matrix, Vector<T> vec);

        // Class specific functions
        /*! Computes the scaling vector (max absolute elements of each row) (auxiliary function used for scaled partial pivoting)
          * 
          * \pre none
          * \post scaling vector calculated and stored in `m_scales`
        */
        void calculate_scales();

        /*! Computes the ratio vector (leading entry in column `col` / associated entry in scaling vector) (auxiliary function used for scaled partial pivoting)
          * 
          * \param col the column we are about to evaluate
          * 
          * \pre (/) operator defined for type T
          * \pre scaling vector `m_scales` should be populated
          * \pre zero should not be the max element in a row and also the leading entry in column `col`
          * \post ratio vector is calculated and stored in `m_ratios`
          * 
          * \throws std::domain_error is thrown if scaling vector is not populated
          * \throws std::domain_error is thrown if 0/0 division occurs
        */
        void calculate_ratios(const int& col);

        /*! Swaps the row with highest ratio value with row `current_row`  (auxiliary function used for scaled partial pivoting)
          * 
          * \param current_row the row we are currently evaluating 
          * \param vec the vector that is being used in tandem with the matrix for gaussian elimination
          * 
          * \pre ratio vector `m_ratios` should be populated
          * \post row with highest ratio is swapped with row `current_row` if these rows are different
          * 
          * \throws std::domain_error is thrown if ratio vector is not populated
        */
        void rearrange(const int& current_row, Vector<T>& vec);

        /*! Reduce below row/column `row_col`  (auxiliary function)
          * 
          * \param row_col the entry in the diagonal that we are about to row reduce 
          * \param vec the vector that is being used in tandem with the matrix for gaussian elimination
          * 
          * \pre (==) operator defined for type T and numeric value '0'
          * \pre (/) operator defined for type T
          * \post all elements in matrix below diagonal entry `row_col` are "zeroed" out and the rows are updated by the rules of gaussian elimination 
        */
        void row_reduce(const int& row_col, Vector<T>& vec);
};

#include "gaussian_solver.hpp"
#endif