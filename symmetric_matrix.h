/*! \file
 *  Symmetric_Matrix class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef SYMMETRIC_MATRIX_H
#define SYMMETRIC_MATRIX_H
#include "base_matrix.h"

/*! Symmetric square matrix class */
template <class T>
class Symmetric_Matrix : public virtual General_Matrix<T> {
    private:
        Status m_state;
    
    public:
        /*! Constructs an empty symmetric matrix
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of size 3
         */
        Symmetric_Matrix();

        /*! Constructs an empty symmetric matrix with assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of size `size`
         */
        Symmetric_Matrix(const int& size);

        /*! Constructs a symmetric matrix with default value and assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of maximum capacity `size`
         */
        Symmetric_Matrix(const int& size, const T& default_val);

        /*! Copy constructor
         * 
         * \pre T has a default constructor
         * \post a new symmetric matrix is created with (deep) copies of all elements from source,
         *       and storage size
         */
        Symmetric_Matrix(const Symmetric_Matrix& source);

        /*! Copy constructor (from base)
         * 
         * \pre T has a default constructor
         * \post a new symmetric matrix is created with (deep) copies of all elements from source,
         *       and storage size
         * \post data that is not apart of a symmetric matrix is discarded
         */
        Symmetric_Matrix(const Base_Matrix<T>& source);

        /*! Move constructor
         * 
         * \pre T has a default constructor
         * \post the "moved-from" object is valid and destructible
         */
        Symmetric_Matrix(Symmetric_Matrix&& other);

        //////////////////////////////////////// Matrix Operators ////////////////////////////////////////
        /*! Copies the size and elements of source into *this
         *
         *  \param source the symmetric matrix to copy elements from
         *  \return a reference to the symmetric matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        Symmetric_Matrix& operator=(const Symmetric_Matrix& source);

        /*! Adds the respective elements of this container and m2's container (more efficient than base funct.)
         *
         *  \param m2 the symmetric matrix to add to *this
         *  \return a new symmetric matrix containing sum of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        Symmetric_Matrix<T> operator+(Symmetric_Matrix<T> m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator+;

        /*! Subtracts the respective elements of m2 from this container (more efficient than base funct.)
         *
         *  \param m2 the symmetric matrix to subtract from *this
         *  \return a new symmetric matrix containing difference of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */

        Symmetric_Matrix<T> operator-(Symmetric_Matrix m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator-;

        /*! Multiplies all elements of this container with the value of a scalar quantitiy using simplified approach for symmetric matrices
         *
         *  \param scalar the scalar value to multiply by
         *  \return a new diagonal matrix with all elements having been multiplied by `scalar`
         * 
         *  \pre  (*) operator is defined for type T
         *  \pre  (=) operator is defined for type T
         *  \post (see return)
         */
        Symmetric_Matrix<T> operator*(const T& scalar) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator*;

        /*! Calculates the mathematical transpose of the underlying container
         *
         *  \return a transpose of *this's matrix
         * 
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         */
        virtual General_Matrix<T> transpose() const { return *this; }

        ///////////////////////////////////////// Setters ////////////////////////////////////////////
        /*! Inserts a vector into the square matrix (if not full) (auxiliary function for operator>>)
         * 
         *  \pre  m_size != m_max
         *  \post vector is inserted into matrix and size increases by one
         *  \throws std::domain_error is thrown if m_size == m_max
         */
        virtual void insert_vector(const Vector<T>& vec);

        /*! Sets value at row `row` and column `col` of the symmetric matrix to `val`
         *
         *  \pre `col` <= `row`
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
         *  \pre `col` <= `row`
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
};

/*! Stream extraction operator for `Symmetric_Matrix`. Data is read in as if it is in lower-triangular matrix form. Any matrix members that would be "zeroes" are discarded.
 * \pre stream extraction operator is defined for `T`
 * \post the contents of the matrix are extracted from the input stream into matrix_in
 * \return the modified input stream
 *
 * \relatesalso Symmetric_Matrix
 */
template <typename T>
istream& operator>>(istream& in, Symmetric_Matrix<T>& vec_in);

#include "symmetric_matrix.hpp"
#endif