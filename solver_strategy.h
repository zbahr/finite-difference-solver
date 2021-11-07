/*! \file
 *  Solver_Strategy class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef SOLVER_STRATEGY_H
#define SOLVER_STRATEGY_H
#include "base_matrix.h"

/*! Generic solver class */
template <class T>
class Solver_Strategy { 
    public:
        /*! Pure virtual function for specific solver methods to implement */
        virtual Vector<T> solve(const Base_Matrix<T>& matrix, Vector<T> vec) = 0;

        /*! Virtual destructor */
        virtual ~Solver_Strategy() {}
};

#endif