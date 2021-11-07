/*! \file
 *  Function definitions for the `U_Triangle_Matrix` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
U_Triangle_Matrix<T>::U_Triangle_Matrix() : General_Matrix<T>(), m_state(row_reduced) {}

template <typename T>
U_Triangle_Matrix<T>::U_Triangle_Matrix(const int& size) : General_Matrix<T>(size), m_state(row_reduced) {}

template <typename T>
U_Triangle_Matrix<T>::U_Triangle_Matrix(const U_Triangle_Matrix<T>& source) : General_Matrix<T>(source.m_size), m_state(row_reduced) {
    this->m_elements = source.m_elements;
    this->m_size = source.m_size;
}

template <typename T>
U_Triangle_Matrix<T>::U_Triangle_Matrix(const Base_Matrix<T>& source) : General_Matrix<T>(source.get_size()), m_state(row_reduced) {
    for (int row = 0; row < this->m_max; row++) {
        Vector<T> temp(this->m_max - row);

        for (int col = row; col < this->m_max; col++) {
            temp.push_back(source.get_element(row, col));
        }

        this->m_elements.push_back(temp);
        this->m_size++;
    }
} 

template <typename T>
U_Triangle_Matrix<T>::U_Triangle_Matrix(U_Triangle_Matrix<T>&& other) : General_Matrix<T>(other), m_state(other.m_state) {}

template <typename T>
U_Triangle_Matrix<T>& U_Triangle_Matrix<T>::operator=(const U_Triangle_Matrix<T>& source) {
    this->m_size = source.get_size();
    this->m_max = source.get_max();
    this->m_elements = source.get_elements();
    m_state = row_reduced;
    
    return (*this);
}

template <typename T>
U_Triangle_Matrix<T> U_Triangle_Matrix<T>::operator+(const U_Triangle_Matrix& m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be added must be of same size."); }

    U_Triangle_Matrix result(m2);

    for (int row = 0; row < this->m_size; row++) {
        result.get_row_ref(row) += this->m_elements[row];
    }

    return result;
}

template <typename T>
U_Triangle_Matrix<T> U_Triangle_Matrix<T>::operator-(const U_Triangle_Matrix<T>& m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be subtracted must be of same size."); }

    U_Triangle_Matrix<T> result(m2);

    for (int row = 0; row < this->m_size; row++) {
        result.get_row_ref(row) = this->m_elements[row] - result.get_row_element(row);
    }

    return result;
}

template <typename T>
U_Triangle_Matrix<T> U_Triangle_Matrix<T>::operator*(const U_Triangle_Matrix<T>& m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Square matrices to be multiplied must be of same size."); }

    U_Triangle_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = row; col < this->m_size ; col++) {
            T result(0);
            for (int runner = row; runner < col + 1; runner++) {
                result += get_element(runner, col) * m2.get_element(row, runner);
            }
            result_matrix.set_element(row, col, result);
        }
    }

    return result_matrix;
}

template <typename T>
U_Triangle_Matrix<T> U_Triangle_Matrix<T>::operator*(const T& scalar) const {
    U_Triangle_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = row; col < this->m_size; col++) {
            result_matrix.set_element(row, col, get_element(row, col) * scalar);
        }
    }

    return result_matrix;
}

template <typename T>
Vector<T> U_Triangle_Matrix<T>::operator*(const Vector<T>& vec) const {
    if (this->m_size != vec.get_size()) { throw domain_error("Error: Square matrix to be multiplied by vector must have same number of rows"); }

    Vector<T> result_vector(vec);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = row; col < this->m_size; col++) {
            T result(0);
            for (int runner = row; runner < col + 1; runner++) {
                result += get_element(row, runner) * vec[runner];
            }
            result_vector[row] = result;
        }
    }

    return result_vector;
}

template <typename T>
General_Matrix<T> U_Triangle_Matrix<T>::transpose() const {
    General_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = 0; col < row; col++) {
            result_matrix.set_element(col, row, 0);
        }
        for (int col = row; col < this->m_size; col++) {
            result_matrix.set_element(col, row, get_element(row, col));
        }
    }

    return result_matrix;
}

template <typename T>
void U_Triangle_Matrix<T>::set_element(const int& row, const int& col, const T& val) {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range.");
    }
    if (col < row) { 
        throw out_of_range("Error: Invalid matrix location for upper triangular matrix."); 
    }
    else { this->m_elements[row][col - row] = val; }

    return;
}

template <typename T>
T U_Triangle_Matrix<T>::get_element(const int& row, const int& col) const {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { throw out_of_range("Error: Attempt to access matrix member out of range."); }

    return (col < row ? 0 : this->m_elements[row][col - row]);
}

template <typename T>
Vector<Vector<T> > U_Triangle_Matrix<T>::get_gen_elements() const {
    Vector<Vector<T> > elements(this->m_size);

    for (int row = 0; row < this->m_size; row++) {
        Vector<T> temp(this->m_size);

        for (int col = 0; col < row; col++) {
            temp.push_back(0);
        }

        for (int col = row; col < this->m_size; col++) {
            temp.push_back(get_element(row, col));
        }
    
        elements.push_back(temp);
    }

    return elements;
}

template <typename T>
Vector<T> U_Triangle_Matrix<T>::back_sub(const Vector<Vector<T>>& matrix_data, const Vector<T>& vec) const {
    Vector<T> result(vec);

    for (int row = vec.get_size() - 1; row >= 0; row--) {
        result[row] = vec[row];

        for (int col = row + 1; col < vec.get_size(); col++) {
            result[row] +=  (-matrix_data[row][col - row]) * result[col];
        }

        if (matrix_data[row][0] == 0) {
            throw domain_error("Error: Division by zero during back substitution.");
        }

        result[row] = result[row] / matrix_data[row][0];
    }

    return result;
}

template <typename T>
istream& operator>>(istream& in, U_Triangle_Matrix<T>& matrix_in) {
    for (int i = 0; i < matrix_in.get_max(); i++) {
        // Discard the zeroes
        T zeroes;
        for (int j = i; j > 0; j--) {
            in >> zeroes;
        }

        Vector<T> temp(matrix_in.get_max() - i);
        in >> temp;
        matrix_in.insert_vector(temp);

        // Ensure there is enough data (only true if whole rows are missing)
        if (i + 1 == matrix_in.get_max() && in.peek() == ifstream::traits_type::eof()) {
            throw runtime_error("Error: Insufficient amount of data provided.");
        }
    }

    return in;
}