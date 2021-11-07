/*! \file
 *  Function definitions for the `Symmetric_Matrix` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
Symmetric_Matrix<T>::Symmetric_Matrix() : General_Matrix<T>(), m_state(symmetric) {}

template <typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(const int& size) : General_Matrix<T>(size), m_state(symmetric) {}

template <typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(const int& size, const T& default_val) : General_Matrix<T>(size), m_state(symmetric) {
    for (int i = 0; i < size; i++) {
        insert_vector(Vector<T>(i + 1, default_val));
    }
}

template <typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(const Symmetric_Matrix<T>& source) : General_Matrix<T>(source.m_size), m_state(symmetric) {
    this->m_elements = source.m_elements;
    this->m_size = source.m_size;
}

template <typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(const Base_Matrix<T>& source) : General_Matrix<T>(source.get_size()), m_state(symmetric) {
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
Symmetric_Matrix<T>::Symmetric_Matrix(Symmetric_Matrix<T>&& other) : General_Matrix<T>(other), m_state(other.m_state) {} 

template <typename T>
Symmetric_Matrix<T>& Symmetric_Matrix<T>::operator=(const Symmetric_Matrix<T>& source) {
    this->m_size = source.get_size();
    this->m_max = source.get_max();
    this->m_elements = source.get_elements();
    m_state = symmetric;
    
    return (*this);
}

template <typename T>
Symmetric_Matrix<T> Symmetric_Matrix<T>::operator+(Symmetric_Matrix<T> m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be added must be of same size."); }

    for (int row = 0; row < this->m_size; row++) {
        m2.get_row_ref(row) += this->m_elements[row];
    }

    return m2;
}

template <typename T>
Symmetric_Matrix<T> Symmetric_Matrix<T>::operator-(Symmetric_Matrix<T> m2) const {
    if (this->m_size != m2.get_size()) { throw domain_error("Error: Matrices to be added must be of same size."); }

    for (int row = 0; row < this->m_size; row++) {
        m2.get_row_ref(row) = this->m_elements[row] - m2.get_row_element(row);
    }

    return m2;
}

template <typename T>
Symmetric_Matrix<T> Symmetric_Matrix<T>::operator*(const T& scalar) const {
    Symmetric_Matrix<T> result_matrix(*this);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = 0; col < row + 1; col++) {
            result_matrix.set_element(row, col, get_element(row, col) * scalar);
        }
    }

    return result_matrix;
}

template <typename T>
void Symmetric_Matrix<T>::insert_vector(const Vector<T>& vec) {
    if (this->m_size == this->m_max) { throw domain_error("Error: Square matrix is already full."); }

    this->m_elements.push_back(vec);
    this->m_size++;

    return;
}

template <typename T>
void Symmetric_Matrix<T>::set_element(const int& row, const int& col, const T& val) {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }
    if (col > row) { 
        this->m_elements[col][row] = val;
    }
    else { 
        this->m_elements[row][col] = val; 
    }

    return;
}

template <typename T>
T Symmetric_Matrix<T>::get_element(const int& row, const int& col) const {
    if (row < 0 || row >= this->m_size || col < 0 || col >= this->m_size) { 
        throw out_of_range("Error: Attempt to access matrix member out of range."); 
    }

    return (col > row ? this->m_elements[col][row] : this->m_elements[row][col]);
}

template <typename T>
Vector<Vector<T> > Symmetric_Matrix<T>::get_gen_elements() const {
    Vector<Vector<T> > elements(this->m_elements);

    for (int row = 0; row < this->m_size; row++) {
        for (int col = row + 1; col < this->m_size; col++) {
            elements[row].push_back(get_element(col, row));
        }
    }

    return elements;
}

template <typename T>
istream& operator>>(istream& in, Symmetric_Matrix<T>& matrix_in) {
    for (int i = 0; i < matrix_in.get_max(); i++) {
        Vector<T> temp(i + 1);
        in >> temp;
        matrix_in.insert_vector(temp);

        // Discard the symmetric entries
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