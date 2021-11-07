/*! \file
 *  Function definitions for the `L_Triangle_Matrix` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
L_Triangle_Matrix<T>::L_Triangle_Matrix() : General_Matrix<T>(), m_state(row_reduced) {}

template <typename T>
L_Triangle_Matrix<T>::L_Triangle_Matrix(const int& size) : General_Matrix<T>(size), m_state(row_reduced) {}

template <typename T>
L_Triangle_Matrix<T>::L_Triangle_Matrix(const L_Triangle_Matrix<T>& source) : General_Matrix<T>(source.m_size), m_state(row_reduced) {
    this->m_elements = source.m_elements;
    this->m_size = source.m_size;
}

template <typename T>
L_Triangle_Matrix<T>::L_Triangle_Matrix(const Base_Matrix<T>& source) : General_Matrix<T>(source.get_size()), m_state(row_reduced) {
    for (int row = 0; row < this->m_max; row++) {
        Vector<T> temp(row + 1);

        for (int col = 0; col < temp.get_max(); col++) {
            temp.push_back(source.get_element(row, col));
        }

        this->m_elements.push_back(temp);
        this->m_size++;
    }
} 

template <typename T>
L_Triangle_Matrix<T>::L_Triangle_Matrix(L_Triangle_Matrix<T>&& other) : General_Matrix<T>(other), m_state(other.m_state) {}

template <typename T>
L_Triangle_Matrix<T>& L_Triangle_Matrix<T>::operator=(const L_Triangle_Matrix<T>& source) {
    this->m_size = source.get_size();
    this->m_max = source.get_max();
    this->m_elements = source.get_elements();
    m_state = row_reduced;
    
    return (*this);
}

template <typename T>
L_Triangle_Matrix<T> L_Triangle_Matrix<T>::operator+(L_Triangle_Matrix<T> m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be added must be of same size."); }

    for (int row = 0; row < this->m_size; row++) {
        m2.get_row_ref(row) += this->m_elements[row];
    }

    return m2;
}

template <typename T>
L_Triangle_Matrix<T> L_Triangle_Matrix<T>::operator-(L_Triangle_Matrix<T> m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be subtracted must be of same size."); }

    for (int i = 0; i < this->m_size; i++) {
        m2.get_row_ref(i) = this->m_elements[i] - m2.get_row_element(i);
    }

    return m2;
}

template <typename T>
L_Triangle_Matrix<T> L_Triangle_Matrix<T>::operator*(const L_Triangle_Matrix<T>& m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Square matrices to be multiplied must be of same size."); }

    L_Triangle_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = 0; col < row + 1; col++) {
            T result(0);
            for (int runner = col; runner < row + 1; runner++) {
                result += get_element(row, runner) * m2.get_element(runner, col);
            }
            result_matrix.set_element(row, col, result);
        }
    }

    return result_matrix;
}

template <typename T>
L_Triangle_Matrix<T> L_Triangle_Matrix<T>::operator*(const T& scalar) const {
    L_Triangle_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = 0; col < row + 1; col++) {
            result_matrix.set_element(row, col, get_element(row, col) * scalar);
        }
    }

    return result_matrix;
}

template <typename T>
Vector<T> L_Triangle_Matrix<T>::operator*(const Vector<T>& vec) const {
    if (this->m_size != vec.get_size()) { throw domain_error("Error: Square matrix to be multiplied by vector must have same number of rows"); }

    Vector<T> result_vector(vec);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = 0; col < row + 1; col++) {
            T result(0);
            for (int runner = 0; runner < row + 1; runner++) {
                result += get_element(row, runner) * vec[runner];
            }
            result_vector[row] = result;
        }
    }

    return result_vector;
}

template <typename T>
General_Matrix<T> L_Triangle_Matrix<T>::transpose() const {
    General_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        int row_length = this->m_elements[row].get_size();
        for (int col = 0; col < row_length; col++) {
            result_matrix.set_element(col, row, this->m_elements[row][col]);
        }
        for (int col = row_length; col < this->m_size; col++) {
            result_matrix.set_element(col, row, 0);
        }
    }

    return result_matrix;
}

template <typename T>
void L_Triangle_Matrix<T>::set_element(const int& row, const int& col, const T& val) {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }
    if (col > row) { 
        throw out_of_range("Error: Invalid matrix location for lower triangular matrix."); 
    }
    else { this->m_elements[row][col] = val; }

    return;
}

template <typename T>
T L_Triangle_Matrix<T>::get_element(const int& row, const int& col) const {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }

    return (col > row ? 0 : this->m_elements[row][col]);
}

template <typename T>
Vector<Vector<T> > L_Triangle_Matrix<T>::get_gen_elements() const {
    Vector<Vector<T> > elements(this->m_elements);

    for (int i = 0; i < this->m_size; i++) {
        for (int j = i + 1; j < this->m_size; j++) {
            elements[i].push_back(0);
        }
    }

    return elements;
}

template <typename T>
Vector<T> L_Triangle_Matrix<T>::back_sub(const Vector<Vector<T> >& matrix_data, const Vector<T>& vec) const {
    Vector<T> result(vec);

    for (int row = 0; row < this->m_size; row++) {
        result[row] = vec[row];

        for (int col = 0; col < row; col++) {
            result[row] +=  (-matrix_data[row][col]) * result[col];
        }

        if (get_element(row, row) == 0) {
            throw domain_error("Error: Division by zero during back substitution.");
        }

        result[row] = result[row] / matrix_data[row][row];
    }

    return result;
}

template <typename T>
istream& operator>>(istream& in, L_Triangle_Matrix<T>& matrix_in) {
    for (int i = 0; i < matrix_in.get_max(); i++) {
        Vector<T> temp(i + 1);
        in >> temp;
        matrix_in.insert_vector(temp);

        // Discard the zeroes
        T zeroes;
        for (int j = matrix_in.get_max() - (i + 1); j > 0; j--) {
            in >> zeroes;
        }

        // Ensure there is enough data (only true if whole rows are missing)
        if (i + 1 == matrix_in.get_max() && in.peek() == ifstream::traits_type::eof()) {
            throw runtime_error("Error: Insufficient amount of data provided.");
        }
    }

    return in;
}
