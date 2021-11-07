/*! \file
 *  General_Matrix class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H
#include "base_matrix.h"

/*! General square matrix class */
template <class T>
class General_Matrix : public Base_Matrix<T> {
    private:
        /*! State of the matrix (none, row-reduced, symmetric, or tridiagonal) */
        Status m_state;
        
    protected:
        /*! Size of the matrix */
        int m_size;

        /*! Maximum size of the underlying container */
        int m_max;

        /*! Vector of vectors containing matrix data */
        Vector<Vector<T>> m_elements;
    
    public:
        /*! Constructs an empty general matrix
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of size 3
         */
        General_Matrix();

        /*! Constructs an empty general matrix with assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of maximum capacity `size`
         */
        General_Matrix(const int& size);

        /*! Constructs a general matrix with default value and assigned allocation storage
         * 
         * \pre  T has a default constructor
         * \post the matrix is square of maximum capacity `size`
         */
        General_Matrix(const int& size, const T& default_val);

        /*! Copy constructor
         *
         *  \pre T has a default constructor
         *  \post a new general matrix is created with (deep) copies of all elements from source,
         *       and storage size
         */
        General_Matrix(const General_Matrix<T>& source);

        /*! Copy constructor (from base)
         * 
         * \pre T has a default constructor
         * \post a new general matrix is created with (deep) copies of all elements from source,
         *       and storage size
         */
        General_Matrix(const Base_Matrix<T>& source);

        /*! Move constructor
         * 
         * \pre T has a default constructor
         * \post the "moved-from" object is valid and destructible
         */
        General_Matrix(General_Matrix&& other);

        //////////////////////////////////////// Matrix Operators ////////////////////////////////////////
        /*! Copies the size and elements of source into *this (copy and swap idiom)
         *
         *  \param source the general matrix to copy elements from
         *  \return a reference to the general matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        General_Matrix& operator=(General_Matrix<T>& source);

        /*! Copies the size and elements of source into *this
         *
         *  \param source the base matrix to copy elements from
         *  \return a reference to the general matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        virtual General_Matrix& operator=(const Base_Matrix<T>& source);

        /*! Adds the respective elements of this container and m2's container
         *
         *  \param m2 the base matrix to add to *this
         *  \return a new general matrix containing sum of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        virtual General_Matrix operator+(const Base_Matrix<T>& m2) const;

        /*! Subtracts the respective elements of m2 from this container
         *
         *  \param m2 the base matrix to subtract from *this
         *  \return a new general matrix containing difference of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != m2.m_size
         */
        virtual General_Matrix operator-(const Base_Matrix<T>& m2) const;

        /*! Multiplies the respective elements of this container and m2's container following the rules of matrix multiplication
         *
         *  \param m2 the base matrix to multiply by
         *  \return a new general matrix of size `size`, containing product of the two matrices (*this and m2)
         * 
         *  \pre  m_size == m2.m_size
         *  \pre  T has a paramterized constructor for numeric value '0'
         *  \pre  (+=) operator defined for type T
         *  \pre  (*) operator defined for type T
         *  \pre  (=) operator defined for type T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != vm.m_size
         */
        virtual General_Matrix operator*(const Base_Matrix<T>& m2) const;

        /*! Multiplies all elements of this container with the value of a scalar quantitiy
         *
         *  \param scalar the scalar value to multiply by
         *  \return a new general matrix with all elements having been multiplied by `scalar`
         * 
         *  \pre  (*) operator is defined for type T
         *  \pre  (=) operator is defined for type T
         *  \post (see return)
         */
        General_Matrix operator*(const T& scalar) const;

        /*! Multiplies the respective elements of this container and vec's container following the rules of matrix-vector multiplication
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
        /*! Inserts a vector into the square matrix (if not full) (auxiliary function for operator>>)
         * 
         *  \pre  m_max == vec.m_size
         *  \pre  m_size != m_max
         *  \post vector is inserted into matrix and size increases by one
         *  \throws std::domain_error is thrown if m_size != vec.m_size
         *  \throws std::domain_error is thrown if m_size == m_max
         */
        virtual void insert_vector(const Vector<T>& vec);

        /*! Sets value at row `row` and column `col` of the matrix to `val`
         *
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
        /* Replacement for operator[] 
         *  
         *  \pre `row` >= 0
         *  \pre `row` < m_size
         *  \post (see return)
         *  \returns a reference to a `Vector` corresponding to the row of the matrix
         * 
        */
        virtual Vector<T>& get_row_ref(const int& row) { return m_elements[row]; }

        /* Replacement for operator[] const
         *
         *  \pre `row` >= 0
         *  \pre `row` < m_size
         *  \post (see return)
         *  \returns a const reference to a `Vector` corresponding to the row of the matrix
         */
        virtual const Vector<T>& get_row_element(const int& row) const { return m_elements[row];}

        /*! Gets value at row `row` and column `col` of the matrix
         *
         *  \pre `row` >= 0
         *  \pre `row` < m_size
         *  \pre `col` >= 0
         *  \pre `col` < m_size
         *  \post none
         *  \throws std::out_of_range if pre-conditions are broken
         * 
         */
        virtual T get_element(const int& row, const int& col) const;

        /*! Getter for number of elements (rows) in the matrix 
         * 
         *  \return how many rows are in the matrix
         */
        virtual int get_size() const { return m_size; }

        /*! Getter for maximum capacity of the matrix
         * 
         *  \return maximum capacity of matrix
         */
        virtual int get_max() const { return m_max; }

        /*! Accessor for the underlying matrix representation of the data 
         *
         *  \pre  none
         *  \post (see return)
         *  \return a `Vector` of rows (`Vector`'s) representing the matrix is returned
         * 
         */
        virtual Vector<Vector<T> > get_elements() const { return m_elements; }

        /*! Accessor for the general matrix representation of the data
         *
         *  \pre  none
         *  \post (see return)
         *  \return a `Vector` of rows (`Vector`'s) representing the `General_Matrix` representation of the matrix
         * 
         */
        virtual Vector<Vector<T> > get_gen_elements() const { return m_elements; }

        ///////////////////////////// Solver auxiliary functions //////////////////////////////////
        /*! Performs back-substitution on the row-reduced `m_gauss_result`
          * 
          * \return a new vector `x` that solves the equation `m_elements` * `x` = `vec`
          * 
          * \pre (=) operator defined for type T
          * \pre (*) operator defined for type T
          * \pre (+=) operator defined for type T
          * \pre (-, unary) operator defined for type T
          * \pre (/) operator defined for type T
          * \pre (==) operator defined for type T with numeric types (0)
          * \pre gaussify() has been called prior to this function's call
          * \post (see return) 
          * 
          * \throws std::domain_error is thrown if gaussian elimination has not been called
          * \throws std::domain_error is thrown if division by zero occurs
        */
        virtual Vector<T> back_sub(const Vector<Vector<T>>& reduced_result, const Vector<T>& vec) const; 

        //////////////////////////////////// State auxiliary functions //////////////////////////////////////
        /*! Gets the current status of the matrix (none, row-reduced, symmetric, etc.)
         *
         *  \pre none
         *  \post (see return)
         *  \return the status of the matrix
         * 
         */
        virtual Status get_status() const { return m_state; }

        /*! Determines whether this matrix is symmetric (auxiliary function)
         * 
         *  \pre T is comparable to numeric value '0'
         *  \post (see return)
         *  \return true if matrix is symmetric; false otherwise
         * 
         */
        bool is_symmetric() const;

        /*! Determines whether this matrix is tridiagonal (auxiliary function)
         * 
         *  \pre T is comparable to numeric value '0'
         *  \post (see return)
         *  \return true if matrix is tridiagonal; false otherwise
         * 
         */
        bool is_tridiagonal() const;

        /*! Function to enable operator<< for use from a base class pointer
         *
         *  \pre operator<< is defined for `General_Matrix`
         *  \post *this is passed to the output stream
         * 
         */
        virtual void print() const { cout << *this; }

        //////////////////////////////////// Friend functions ///////////////////////////////////////
        #ifndef DOXYGEN
        template<typename U>
        friend void swap(General_Matrix<U>& v1, General_Matrix<U>& v2);
        #endif
};

/*! Auxiliary swap function for `General_Matrix` (copy and swap idiom)
  *
  * \pre none
  * \post the contents of m1 and m2 have been swapped
  *
  * \relatesalso General_Matrix
  */
template <typename T>
void swap(General_Matrix<T>& v1, General_Matrix<T>& v2);

/*! Stream insertion operator for `General_Matrix`.
 *
 * \pre stream insertion operator is defined for `T`
 * \post the contents of the matrix are printed to the ouptut stream
 * \return the modified output stream
 *
 * \relatesalso General_Matrix
 */
template <typename T>
ostream& operator<<(ostream& out, const General_Matrix<T>& matrix_out);

/*! Stream extraction operator for `General_Matrix`.
 *
 * \pre stream extraction operator is defined for `T`
 * \post the contents of the matrix are extracted from the input stream into matrix_in
 * \return the modified input stream
 *
 * \relatesalso General_Matrix
 */
template <typename T>
istream& operator>>(istream& in, General_Matrix<T>& vec_in);

#include "general_matrix.hpp"
#endif