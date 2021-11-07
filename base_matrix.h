/*! \file
 *  Base_Matrix class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H
#include "vector.h"

// Forward declaration of general matrix class
#ifndef DOXYGEN
template <class T>
class General_Matrix;
#endif

/*! Flag to be used by the matrix solver */
enum Status { none, row_reduced, symmetric, tridiagonal };

/*! Abstract `Base_Matrix` class definition. See child classes for detailed function documentation.*/
template <class T>
class Base_Matrix {
    public:
        /////////////////////////////////////// Matrix Operations /////////////////////////////////////////
        /*! Copies the size and elements of source into *this
         *
         *  \param source the matrix to copy elements from
         *  \return a reference to the general matrix that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        virtual General_Matrix<T>& operator=(const Base_Matrix<T>& source) = 0;

        /*! Adds the respective elements of this container and m2's underlying container */
        virtual General_Matrix<T> operator+(const Base_Matrix<T>& m2) const = 0;

        /*! Subtracts the respective elements of m2 from the underlying container */
        virtual General_Matrix<T> operator-(const Base_Matrix<T>& m2) const = 0;

        /*! Multiplies the respective elements of this container and m2's container following the rules of matrix multiplication */
        virtual General_Matrix<T> operator*(const Base_Matrix<T>& m2) const = 0;

        /*! General purpose matrix-vector multiplication function */
        virtual Vector<T> operator*(const Vector<T>& vec) const = 0;

        /*! Returns a general matrix that is the transpose of the underlying matrix */
        virtual General_Matrix<T> transpose() const = 0;

        //////////////////////////////////////////// Setters //////////////////////////////////////////////
        /*! Replacement for operator[] (used to access a reference to a vector of a matrix) */
        virtual Vector<T>& get_row_ref(const int& row) = 0;

        /*! Inserts a vector (row) into a matrix. Used as an auxiliary function for extraction operator.*/
        virtual void insert_vector(const Vector<T>& vec) = 0;

        /*! Used to set an individual element of a matrix.*/
        virtual void set_element(const int& row, const int& col, const T& val) = 0;

        ////////////////////////////////////////// Accessors //////////////////////////////////////////////
        /*! Replacement for const operator[] (used to access a vector of a matrix) */
        virtual const Vector<T>& get_row_element(const int& row) const = 0;

        /*! Used to get an individual element of a matrix. */
        virtual T get_element(const int& row, const int& col) const = 0;

        /*! Gets the size of a matrix. */
        virtual int get_size() const = 0;

        /*! Gets the maximum allotted size of a matrix. */
        virtual int get_max() const = 0;

        /*! Gets all of the elements of the underlying matrix. */
        virtual Vector<Vector<T> > get_elements() const = 0;

        /*! Gets the general matrix representation of the underlying matrix */
        virtual Vector<Vector<T> > get_gen_elements() const = 0;

        ///////////////////////////////// Universal solving functions /////////////////////////////////////
        /*! Performs appropriate substitution method with a row-reduced matrix and an associated vector.*/
        virtual Vector<T> back_sub(const Vector<Vector<T>>&, const Vector<T>& vec) const = 0;

        //////////////////////////////////////// State function //////////////////////////////////////////
        /*! Gets the current state of a matrix (row-reduced, symmetric, etc.) */
        virtual Status get_status() const = 0;

        /*! Prints the matrix to std::out */
        virtual void print() const = 0;

        /*! Virtual destructor */
        virtual ~Base_Matrix() {}
};

#endif