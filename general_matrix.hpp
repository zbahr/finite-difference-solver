/*! \file
 *  Function definitions for the `General_Matrix` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
General_Matrix<T>::General_Matrix() : m_state(none), m_size(0), m_max(3), m_elements(3) {}

template <typename T>
General_Matrix<T>::General_Matrix(const int& size) : m_state(none), m_size(0), m_max(size), m_elements(size) {}

template <typename T>
General_Matrix<T>::General_Matrix(const int& size, const T& default_val) : m_state(none), m_size(0), m_max(size), m_elements(size) {
    for (int i = 0; i < size; i++) {
        insert_vector(Vector<T>(size, default_val));
    }
}

template <typename T>
General_Matrix<T>::General_Matrix(const General_Matrix<T>& source) : m_state(none), m_size(source.m_size), m_max(source.m_max), m_elements(source.get_gen_elements()) {}

template <typename T>
General_Matrix<T>::General_Matrix(const Base_Matrix<T>& source): m_state(none), m_size(source.get_size()), m_max(source.get_max()), m_elements(source.get_gen_elements()) {}

template <typename T>
General_Matrix<T>::General_Matrix(General_Matrix<T>&& other) : m_state(none), m_size(other.m_size), m_max(other.m_max), m_elements(other.get_gen_elements()) {} 

template <typename T>
General_Matrix<T>& General_Matrix<T>::operator=(General_Matrix<T>& source) {
    swap(*this, source);
    return (*this);
}

template <typename T>
General_Matrix<T>& General_Matrix<T>::operator=(const Base_Matrix<T>& source) {
    m_size = source.get_size();
    m_max = source.get_max();
    m_elements = source.get_gen_elements();
    m_state = none;
    
    return (*this);
}

template <typename T>
General_Matrix<T> General_Matrix<T>::operator+(const Base_Matrix<T>& m2) const {
    if (m_size != m2.get_size()) { throw domain_error("Error: Matrices to be added must be of same size."); }

    General_Matrix<T> result(m2);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            result.set_element(row, col, get_element(row, col) + m2.get_element(row, col));
        }
    }

    return result;
}

template <typename T>
General_Matrix<T> General_Matrix<T>::operator-(const Base_Matrix<T>& m2) const {
    if (m_size != m2.get_size()) { throw domain_error("Error: Matrices to be subtracted must be of same size."); }

    General_Matrix<T> result(m2);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            result.set_element(row, col, get_element(row, col) - m2.get_element(row, col));
        }
    }

    return result;
}

template <typename T>
General_Matrix<T> General_Matrix<T>::operator*(const Base_Matrix<T>& m2) const {
    if (m_size != m2.get_size()) { throw domain_error("Error: Square matrices to be multiplied must be of same size."); }

    General_Matrix<T> result_matrix(*this);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            T result(0);
            for (int runner = 0; runner < m_size; runner++) {
                result += get_element(row, runner) * m2.get_element(runner, col);
            }

            result_matrix.set_element(row, col, result);
        }
    }

    return result_matrix;
}

template <typename T>
General_Matrix<T> General_Matrix<T>::operator*(const T& scalar) const {
    General_Matrix<T> result_matrix(*this);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            result_matrix.set_element(row, col, result_matrix.get_element(row, col) * scalar);
        }
    }

    return result_matrix;
}

template <typename T>
Vector<T> General_Matrix<T>::operator*(const Vector<T>& vec) const {
    if (m_size != vec.get_size()) { throw domain_error("Error: Square matrix to be multiplied by vector must have same number of rows"); }

    Vector<T> result_vector(vec);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            T result(0);
            for (int runner = 0; runner < m_size; runner++) {
                result += get_element(row, runner) * vec[runner];
            }
            result_vector[row] = result;
        }
    }

    return result_vector;
}

template <typename T>
General_Matrix<T> General_Matrix<T>::transpose() const {
    General_Matrix<T> result_matrix(*this);

    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            result_matrix.set_element(col, row, m_elements[row][col]);
        }
    }

    return result_matrix;
}

template <typename T>
void General_Matrix<T>::insert_vector(const Vector<T>& vec) {
    if (m_max != vec.get_size()) { throw domain_error("Error: Vector to be added to square matrix must be of same dimension as other matrix elements"); }
    if (m_size == m_max) { throw domain_error("Error: Square matrix is already full."); }

    m_elements.push_back(vec);
    m_size++;

    return;
}

template <typename T>
void General_Matrix<T>::set_element(const int& row, const int& col, const T& val) {
    if (row < 0 || row >= m_size || col < 0 || col >= m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }
    else { m_elements[row][col] = val; }

    return;
}

template <typename T>
T General_Matrix<T>::get_element(const int& row, const int& col) const {
    if (row < 0 || row >= m_size || col < 0 || col >= m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }

    return m_elements[row][col];
}

template <typename T>
Vector<T> General_Matrix<T>::back_sub(const Vector<Vector<T>>& reduced_matrix, const Vector<T>& vec) const {
    // Indicator to ensure gaussian elimination was called prior to back substitution
    if (reduced_matrix.get_size() != m_max) { throw domain_error("Error: Back substitution should be applied on a row-reduced matrix."); }

    Vector<T> result(vec);

    for (int row = vec.get_size() - 1; row >= 0; row--) {
        result[row] = vec[row];

        for (int col = row + 1; col < vec.get_size(); col++) {
            result[row] +=  (-reduced_matrix[row][col]) * result[col];
        }

        if (reduced_matrix[row][row] == 0) { throw domain_error("Error: Division by zero during back substitution."); }

        result[row] = result[row] / reduced_matrix[row][row];
    }

    return result;
}

template <typename T>
bool General_Matrix<T>::is_symmetric() const {
    for (int row = 0; row < m_size; row++) {
        for (int col = 0; col < m_size; col++) {
            if (m_elements[row][col] != m_elements[col][row]) return false;
        }
    }

    return true;
}

template <typename T>
bool General_Matrix<T>::is_tridiagonal() const {
    for (int row = 0; row < m_size; row++) {
        if (row == 0) {
            for (int col = 2; col < m_size; col++) {
                if (m_elements[row][col] != 0) return false;
            }
        }
        else if (row == m_size - 1) {
            for (int col = 0; col < m_size - 2; col++) {
                if (m_elements[row][col] != 0) return false;
            }
        }
        else {
            for (int col = 0; col < row - 1; col++) {
                if (m_elements[row][col] != 0) return false;
            }
        }
    }

    return true;
}

template <typename T>
ostream& operator<<(ostream& out, const General_Matrix<T>& matrix_out) {
    for (int i = 0; i < matrix_out.get_size(); i++) {
        for (int j = 0; j < matrix_out.get_size(); j++) {
            out << setprecision(PRECISION) << matrix_out.get_element(i, j) << " ";
        }
        out << endl;
    }

    return out;
}

template <typename T>
istream& operator>>(istream& in, General_Matrix<T>& matrix_in) {
    for (int i = 0; i < matrix_in.get_max(); i++) {
        Vector<T> temp(matrix_in.get_max());
        in >> temp;
        matrix_in.insert_vector(temp);

        // Ensure there is enough data (only true if whole rows are missing)
        if (i + 1 == matrix_in.get_max() && in.peek() == ifstream::traits_type::eof()) {
            throw runtime_error("Error: Insufficient amount of data provided.");
        }
    }

    return in;
}

template <typename T>
void swap(General_Matrix<T>& m1, General_Matrix<T>& m2) {
    std::swap(m1.m_size, m2.m_size);
    std::swap(m1.m_max, m2.m_max);
    std::swap(m1.m_elements, m2.m_elements);
    std::swap(m1.m_state, m2.m_state);

    return;
}