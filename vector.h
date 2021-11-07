/*! \file
 *  Vector class definition/declaration.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef VECTOR_H
#define VECTOR_H
#include "libraries.h"

/*! Vector class adapted from Barton and Nackman's Array class.
 */
template <class T>
class Vector {
    private: 
        int m_size;
        int m_max;
        T* ptr_to_data;

        /*! Auxiliary function to copy elements from source
         * 
         *  \pre  (=) operator defined for T type
         *  \post this's underlying container contains elements source's container
         */
        void vector_copy(const Vector& source);

    public:
        /*! Constructs an empty vector
         * 
         * \pre  T(T) is implemented
         * \post the container is allocated storage for 4 elements
         */
        Vector();

        /*! Constructs an empty vector with assigned allocation storage
         * 
         * \pre  T(T) is implemented
         * \post the container is allocated storage for max_size elements
         */
        Vector(const int& max_size); 

        /*! Constructs a default value vector with assigned allocation storage
         * 
         * \pre  T(T) is implemented
         * \post the container is allocated storage for max_size elements
         */
        Vector(const int& max_size, const T& default_val);
        
        /*! Copy Constructor
         * 
         * \pre T(T) is implemented
         * \post a new vector is created with copies of all elements from source,
         *       and storage size
         */
        Vector(const Vector& source);

        /*! Move constructor
         * 
         * \pre T has a default constructor
         * \post the "moved-from" object is valid and destructible
         */
        Vector(Vector&& other);

        /*! Destructor
         *
         *  \pre none
         *  \post ptr_to_data is de-allocated
         */
        ~Vector();

        // Overloaded operators
        /*! Reads an element of the vector
         *
         *  \param i the index of the element to read
         *  \return the element at position i is returned
         * 
         *  \pre i is within the range [0, m_size]
         *  \post (see return)
         *  \throws std::out_of_range is thrown if i < 0 or i >= m_size
         */
        T& operator[](const int& i);

        /*! Reads an element of the vector
         *
         *  \param i the index of the element to read
         *  \return the element at position i is returned (const)
         * 
         *  \pre i is within the range [0, m_size]
         *  \post (see return)
         *  \throws std::out_of_range is thrown if i < 0 or i >= m_size
         */
        const T& operator[](const int& i) const;

        /*! Negates all elements in the vector
         *
         *  \return a new vector is returned with all elements negated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post None
         */
        Vector operator-() const;

        /*! Copies the size and elements of source into *this
         *
         *  \param source the vector to copy elements from
         *  \return a reference to the vector that has been updated
         * 
         *  \pre (-) unary operator is defined for type T
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        Vector& operator=(Vector source);

        /*! Sets all elements in container
         *
         *  \param val the value to set all elements in the container to
         *  \return a reference to the vector that has been updated
         * 
         *  \pre (=) operator is defined for type T
         *  \post (see return)
         */
        Vector& operator=(const T& val);

        /*! Adds the respective elements of this container and v2's container
         *
         *  \param v2 the vector to add to *this
         *  \return a new vector containing sum of the two vectors (*this and v2)
         * 
         *  \pre  m_size == v2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != v2.m_size
         */
        Vector<T> operator+(const Vector<T>& v2) const;

        /*! Adds the respective elements of this container and v2's container
         *
         *  \param v2 the vector to add to *this
         *  \return a reference to `this` vector with updated entries
         * 
         *  \pre  m_size == v2.m_size
         *  \pre  (+=) operator is defined for type T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != v2.m_size
         */
        Vector<T>& operator+=(const Vector<T>& v2);

        /*! Subtracts the respective elements of v2's container from this container
         *
         *  \param v2 the vector to subtract from *this
         *  \return a new vector containing difference of the two vectors (*this and v2)
         * 
         *  \pre  m_size == v2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != v2.m_size
         */
        Vector<T> operator-(const Vector<T>& v2) const;

        /*! Subtracts the respective elements of v2's container from this container
         *
         *  \param v2 the vector to subtract from *this
         *  \return a reference to `this` vector with updated entries
         * 
         *  \pre  m_size == v2.m_size
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != v2.m_size
         */
        Vector<T>& operator-=(const Vector<T>& v2);

        /*! Calculates the dot product of this container and v2's container
         *
         *  \param v2 the vector to compute the dot product with *this
         *  \return the dot product value
         * 
         *  \pre  m_size == v2.m_size
         *  \pre (+=) operator defined for T
         *  \post (see return)
         *  \throws std::domain_error is thrown if m_size != v2.m_size
         */
        T operator*(const Vector<T>& v2) const;

        /*! Multiplies all entries in `this` container with the value `scalar`
         *
         *  \param scalar the value to multiply entries by
         *  \return a new scaled vector
         * 
         *  \pre  m_size == v2.m_size
         *  \pre (*=) operator defined for T
         *  \post (see return)
         */
        Vector<T> operator*(const T& scalar) const;

        // Setters
        /*! Resets the vector capacity and sets a default value
         * 
         *  \param size the size to reset the vector to
         *  \pre  T() has been defined and is a numeric type
         *  \post the container is now of of size size, every 
         *        element in the container has a default value,
         *        and m_max is double the size
         */
        void reset_vector(const int& size);

        /*! Resets the vector capacity
         * 
         *  \pre  T() has been defined and is a numeric type
         *  \pre  T is not a pointer type (could lead to memory leaks)
         *  \post the container is now of of size zero
         */
        void clear();

        /*! Appends value val to the end of the container
         * 
         *  \param val the value to append
         *  \pre  (=) defined for T
         *  \pre  destructor defined for T
         *  \post container's size has incremented by one, val
         *        has been appended to the end of the container,
         *        new memory is created if vector is full
         */
        void push_back(const T& val);

        /*! Decrements the container size (effectively removing last element)
         * 
         *  \pre  m_size > 0 
         *  \pre  T is not a pointer type (could lead to memory leaks)
         *  \post container's size has been decremented by one
         *  \throws std::out_of_range if m_size <= 0
         */
        void pop_back();

        /*! Determines the vector element with highest absolute value
         * 
         *  \return the element in the container with maximum absolute valuevalue
         *  
         *  \pre  T is defined for the abs() function
         *  \pre  T is defined for the (<) operator
         *  \post (see return)
         */
        T abs_max_element() const;

        // Getters
        /*! Getter for number of elements in the vector
         * 
         *  \return how many elements are in the vector
         */
        int get_size() const { return m_size; }

        /*! Getter for maximum capacity of the vector
         * 
         *  \return maximum capacity of vector
         */
        int get_max() const { return m_max; }

        /*! Getter for pointer to container
         * 
         *  \return the pointer to the underlying container
         */
        T* get_ptr() const { return ptr_to_data; }

        //////////////////////////////////// Friend functions ///////////////////////////////////////
        #ifndef DOXYGEN
        template <typename U>
        friend void swap(Vector<U>& v1, Vector<U>& v2);
        #endif
};

/*! Auxiliary swap function for `Vector` (copy and swap idiom)
 *
 * \pre none
 * \post the contents of v1 and v2 have been swapped
 *
 * \relatesalso Vector
 */
template <typename T>
void swap(Vector<T>& v1, Vector<T>& v2);

/*! Stream insertion operator for `Vector`.
 *
 * \pre stream insertion operator is defined for `T`
 * \post the contents of the vector are printed to the ouptut stream
 * \return the modified output stream
 *
 * \relatesalso Vector
 */
template <typename T>
ostream& operator<<(ostream& out, const Vector<T>& vec_out);

/*! Stream extraction operator for `Vector`.
 *
 * \pre stream extraction operator is defined for `T`
 * \post the contents of the vector are extracted from the input stream into vec_in
 * \return the modified input stream
 *
 * \relatesalso Vector
 */
template <typename T>
istream& operator>>(istream& in, Vector<T>& vec_in);

#include "vector.hpp"
#endif