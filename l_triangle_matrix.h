/*! \file
 *  L_Triangle_Matrix class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef L_TRIANGLE_MATRIX_H
#define L_TRIANGLE_MATRIX_H
#include "general_matrix.h"

/*! Lower triangular square matrix class */
template <class T>
class L_Triangle_Matrix : public virtual General_Matrix<T> {
    private:
        Status m_state;

    public:
        /*! Constructs an empty lower triangular matrix
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of size 3
         */
        L_Triangle_Matrix();

        /*! Constructs an empty lower triangular matrix with assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of max size `size`
         */
        L_Triangle_Matrix(const int& size);

        /*! Copy constructor
         *
         * \pre T has a default constructor
         * \post a new lower triangular matrix is created with (deep) copies of all elements from source,
         *       and storage size
         */
        L_Triangle_Matrix(const L_Triangle_Matrix<T>& source);

        /*! Copy constructor (from base)
         * 
         * \pre T has a default constructor
         * \post a new lower triangular matrix is created with (deep) copies of all elements from source,
         *       and storage size
         * \post data that is not apart of a lower triangular matrix is discarded
         */
        L_Triangle_Matrix(const Base_Matrix<T>& source);

        /*! Move constructor
         * 
         * \pre T has a default constructor
         * \post the "moved-from" object is valid and destructible
         */
        L_Triangle_Matrix(L_Triangle_Matrix&& other);

        //////////////////////////////////////// Matrix Operators ////////////////////////////////////////
        /*! Copies the size and elements of source into *this
         *
         *  \param source the lower triangular matrix to copy elements from
         *  \return a reference to the lower triangular matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        L_Triangle_Matrix<T>& operator=(const L_Triangle_Matrix<T>& source);

        /*! Adds the respective elements of this container and m2's container (more efficient than base funct.)
         *
         *  \param m2 the lower triangular matrix to add to *this
         *  \return a new lower triangular matrix containing sum of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        L_Triangle_Matrix<T> operator+(L_Triangle_Matrix m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator+;

        /*! Subtracts the respective elements of m2 from this container (more efficient than base funct.)
         *
         *  \param m2 the lower triangular matrix to subtract from *this
         *  \return a new lower triangular matrix containing difference of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */

        L_Triangle_Matrix<T> operator-(L_Triangle_Matrix m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator-;

        /*! Multiplies the respective elements of this container and m2's container using the simplified approach for lower triangular matrices
         *
         *  \param m2 the lower triangular matrix to multiply by
         *  \return a new lower triangular matrix of size `size`, containing product of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \pre  T has a paramterized constructor for numeric value '0'
         *  \pre  (+=) operator defined for type T
         *  \pre  (*) operator defined for type T
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        L_Triangle_Matrix<T> operator*(const L_Triangle_Matrix& m2) const;

        /* Redeclare the hidden parent function */
        using General_Matrix<T>::operator*;

        /*! Multiplies all elements of this container with the value of a scalar quantitiy using simplified approach for lower triangular matrices
         *
         *  \param scalar the scalar value to multiply by
         *  \return a new lower triangular matrix with all elements having been multiplied by `scalar`
         * 
         *  \pre  (*) operator is defined for type T
         *  \pre  (=) operator is defined for type T
         *  \post (see return)
         */
        L_Triangle_Matrix<T> operator*(const T& scalar) const;

        /*! Multiplies the respective elements of this container and vec's container using simplified rules of lower triangular matrix-vector multiplication
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
        /*! Sets value at row `row` and column `col` of the lower triangular matrix to `val`
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

        ///////////////////////////// Gaussian auxiliary functions //////////////////////////////////
        /*! Performs back-substitution in reverse on the already row-reduced `m_elements`
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
        virtual Vector<T> back_sub(const Vector<Vector<T> >& matrix_data, const Vector<T>& vec) const; 
};

/*! Stream extraction operator for `L_Triangle_Matrix`. Data is read in as if it is in lower-triangular matrix form. Any matrix members that would be "zeroes" are discarded.
 * \pre stream extraction operator is defined for `T`
 * \post the contents of the matrix are extracted from the input stream into matrix_in
 * \return the modified input stream
 *
 * \relatesalso L_Triangle_Matrix
 */
template <typename T>
istream& operator>>(istream& in, L_Triangle_Matrix<T>& vec_in);

#include "l_triangle_matrix.hpp"
#endif