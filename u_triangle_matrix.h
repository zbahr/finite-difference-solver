/*! \file
 *  U_Triangle_Matrix class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef U_TRIANGLE_MATRIX_H
#define U_TRIANGLE_MATRIX_H
#include "general_matrix.h"

/*! Upper triangular square matrix class */
template <class T>
class U_Triangle_Matrix : public virtual General_Matrix<T> {
    private:
        Status m_state;

    public:
        /*! Constructs an empty upper triangular matrix
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of size 3
         */
        U_Triangle_Matrix();

        /*! Constructs an empty upper triangular matrix with assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of max size `size`
         */
        U_Triangle_Matrix(const int& size);

        /*! Copy constructor
         *
         * \pre T has a default constructor
         * \post a new upper triangular matrix is created with (deep) copies of all elements from source,
         *       and storage size
         */
        U_Triangle_Matrix(const U_Triangle_Matrix<T>& source);

        /*! Copy constructor (from base)
         * 
         * \pre T has a default constructor
         * \post a new upper triangular matrix is created with (deep) copies of all elements from source,
         *       and storage size
         * \post data that is not apart of a upper triangular matrix is discarded
         */
        U_Triangle_Matrix(const Base_Matrix<T>& source);

        /*! Move constructor
         * 
         * \pre T has a default constructor
         * \post the "moved-from" object is valid and destructible
         */
        U_Triangle_Matrix(U_Triangle_Matrix&& other);

        //////////////////////////////////////// Matrix Operators ////////////////////////////////////////
        /*! Copies the size and elements of source into *this
         *
         *  \param source the upper triangular matrix to copy elements from
         *  \return a reference to the upper triangular matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        U_Triangle_Matrix<T>& operator=(const U_Triangle_Matrix<T>& source);

        /*! Adds the respective elements of this container and m2's container (more efficient than base funct.)
         *
         *  \param m2 the upper triangular matrix to add to *this
         *  \return a new upper triangular matrix containing sum of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        U_Triangle_Matrix<T> operator+(const U_Triangle_Matrix& m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator+;

        /*! Subtracts the respective elements of m2 from this container (more efficient than base funct.)
         *
         *  \param m2 the upper triangular matrix to subtract from *this
         *  \return a new upper triangular matrix containing difference of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */

        U_Triangle_Matrix<T> operator-(const U_Triangle_Matrix<T>& m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator-;

        /*! Multiplies the respective elements of this container and m2's container using the simplified approach for upper triangular matrices
         *
         *  \param m2 the upper triangular matrix to multiply by
         *  \return a new upper triangular matrix of size `size`, containing product of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \pre  T has a paramterized constructor for numeric value '0'
         *  \pre  (+=) operator defined for type T
         *  \pre  (*) operator defined for type T
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        U_Triangle_Matrix<T> operator*(const U_Triangle_Matrix& m2) const;
        
        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator*;

        /*! Multiplies all elements of this container with the value of a scalar quantitiy using the simplified approach for upper triangular matrices
         *
         *  \param scalar the scalar value to multiply by
         *  \return a new upper triangular matrix with all elements having been multiplied by `scalar`
         * 
         *  \pre  (*) operator is defined for type T
         *  \pre  (=) operator is defined for type T
         *  \post (see return)
         */
        U_Triangle_Matrix operator*(const T& scalar) const;

        /*! Multiplies the respective elements of this container and vec's container using the simplified rules of upper triangular matrix-vector multiplication
         *
         *  \param vec the vector to multiply by
         *  \return a vector of size m_size that is the result of the product of `this` matrix and `vec`
         * 
         *  \pre  m_size == vec.m_size
         *  \pre  T has a paramterized constructor for numeric value '0'
         *  \pre  (+=) operator defined for type T
         *  \pre  (*) operator defined for type T
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != vec.m_size
         */
        virtual Vector<T> operator*(const Vector<T>& vec) const;

        /*! Calculates the mathematical transpose of the underlying container
         *
         *  \return a transpose of *this's matrix
         * 
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         */
        virtual General_Matrix<T> transpose() const;

        ///////////////////////////////////////// Setters ////////////////////////////////////////////
        /*! Sets value at row `row` and column `col` of the upper triangular matrix to `val`
         *
         *  \pre `col` >= `row`
         *  \pre `row` >= 0
         *  \pre `row` < m_size
         *  \pre `col` >= 0
         *  \pre `col` < m_size
         *  \post element at row `row` and column `col` is updated with `val`
         *  \throws std::out_of_range if pre-conditions are broken
         * 
         */
        virtual void set_element(const int& row, const int& col, const T& val);

        /////////////////////////////////////// Accessors /////////////////////////////////////////
        /*! Gets value at row `row` and column `col` of the matrix
         *
         *  \pre `col` >= `row`
         *  \pre `row` >= 0
         *  \pre `row` < m_size
         *  \pre `col` >= 0
         *  \pre `col` < m_size
         *  \post none
         *  \throws std::out_of_range if pre-conditions 2-4 are broken
         * 
         *  \return 0 if pre-condition 1 is broken; otherwise, value at row `row` and column `col`
         */
        virtual T get_element(const int& row, const int& col) const;

        /*! Accessor for the general matrix representation of the data
         *
         *  \pre  T(0) is defined
         *  \post (see return)
         *  \return a `Vector` of rows (`Vector`'s) representing the `General_Matrix` representation of the matrix
         * 
         */
        virtual Vector<Vector<T> > get_gen_elements() const;

        /*! Gets the current status of the matrix (none, row-reduced, symmetric, etc.)
         *
         *  \pre none
         *  \post (see return)
         *  \return the status of the matrix
         * 
         */
        virtual Status get_status() const { return m_state; }

        ///////////////////////////// Gaussian auxiliary functions //////////////////////////////////
        /*! Performs back-substitution on the already row-reduced `m_elements`
          * 
          * \return a new vector `x` that solves the equation `m_elements` * `x` = `vec`
          * 
          * \pre (=) operator defined for type T
          * \pre (*) operator defined for type T
          * \pre (+=) operator defined for type T
          * \pre (-, unary) operator defined for type T
          * \pre (/) operator defined for type T
          * \pre (==) operator defined for type T with numeric types (0)
          * \post (see return) 
          * 
          * \throws std::domain_error is thrown if division by zero occurs
        */
        virtual Vector<T> back_sub(const Vector<Vector<T>>& matrix_data, const Vector<T>& vec) const; 
};

/*! Stream extraction operator for `U_Triangle_Matrix`. Data is read in as if it is in upper-triangular matrix form. Any matrix members that would be "zeroes" are discarded.
 * 
 * \pre stream extraction operator is defined for `T`
 * \post the contents of the matrix are extracted from the input stream into matrix_in
 * \return the modified input stream
 *
 * \relatesalso U_Triangle_Matrix
 */
template <typename T>
istream& operator>>(istream& in, U_Triangle_Matrix<T>& vec_in);

#include "u_triangle_matrix.hpp"
#endif