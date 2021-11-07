/*! \file
 *  Function definitions for the `Vector` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
Vector<T>::Vector() : m_size(0), m_max(DEFAULT_MAX), ptr_to_data(new T[DEFAULT_MAX]) {}

template <typename T>
Vector<T>::Vector(const int& max_size) : m_size(0), m_max(max_size), ptr_to_data(new T[max_size]) { }

template <typename T>
Vector<T>::Vector(const int& max_size, const T& default_val) : m_size(max_size), m_max(max_size), ptr_to_data(new T[max_size]) {
    for (int i = 0; i < max_size; i++) {
        ptr_to_data[i] = default_val;
    }
}

template <typename T>
Vector<T>::Vector(const Vector<T>& source) {
    m_size = source.get_size();
    m_max = source.get_max();
    ptr_to_data = new T[m_max];
    vector_copy(source);
}

template <typename T>
Vector<T>::Vector(Vector<T>&& other) : m_size(other.m_size), m_max(other.m_max), ptr_to_data(other.ptr_to_data) {
    other.ptr_to_data = nullptr;
}

template <typename T>
void Vector<T>::vector_copy(const Vector<T>& source) {
    T* p = ptr_to_data + m_size;
    T* q = source.get_ptr() + m_size;
    while (p > ptr_to_data) { *--p = *--q; }
    return;
}

template <typename T>
Vector<T>::~Vector() {
    delete [] ptr_to_data;
}

template <typename T>
Vector<T> Vector<T>::operator-() const {
    Vector<T> temp(*this);
    for (int i = 0; i < m_size; i++) {
        temp[i] = -temp[i];
    }

    return temp;
}

template <typename T>
T& Vector<T>::operator[](const int& i) {
    if (i < 0 || i >= m_size) { throw out_of_range("Error: Attempt to access vector member out of range."); }

    return ptr_to_data[i];
}

template <typename T>
const T& Vector<T>::operator[](const int& i) const {
    if (i < 0 || i >= m_size) { throw out_of_range("Error: Attempt to access vector member out of range."); }

    return ptr_to_data[i];
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T> source) {
    swap(*this, source);
    return (*this);
}

template <typename T>
Vector<T>& Vector<T>::operator=(const T& val) {
    for (int i = 0; i < m_size; i++) {
        ptr_to_data[i] = val;
    }

    cout << "test2" << endl;

    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& v2) const {
    if (m_size != v2.get_size()) { throw domain_error("Error: Vectors to be added must be of same size."); }
    Vector<T> added_vector(*this);

    return (added_vector += v2);
}

template <typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& v2) {
    if (m_size != v2.get_size()) { throw domain_error("Error: Vectors to be added must be of same size."); }

    for (int i = 0; i < m_size; i++) {
        ptr_to_data[i] += v2[i];
    }

    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& v2) const {
    if (m_size != v2.get_size()) { throw domain_error("Error: Vectors to be subtracted must be of same size."); }
    Vector<T> subtracted_vector(*this);

    return (subtracted_vector -= v2);
}

template <typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& v2) {
    if (m_size != v2.get_size()) { throw domain_error("Error: Vectors to be subtracted must be of same size."); }

    return (*this += (-v2));
}

template <typename T>
T Vector<T>::operator*(const Vector<T>& v2) const {
    if (m_size != v2.get_size()) { throw domain_error("Error: Vectors dot product is applied must be of same size."); }
    T dot_product = 0;

    for (int i = 0; i < m_size; i++) {
        dot_product += ptr_to_data[i] * v2[i];
    }

    return dot_product;
}

template <typename T>
Vector<T> Vector<T>::operator*(const T& scalar) const {
    Vector<T> result(*this);

    for (int i = 0; i < m_size; i++) {
        result[i] *= scalar;
    }

    return result;
}

template <typename T>
void Vector<T>::reset_vector(const int& size) {
    if (size != m_size) {
        delete [] ptr_to_data;
        m_size = size;
        ptr_to_data = new T[size];
    }

    m_max = size * 2;

    return;
}

template <typename T>
void Vector<T>::clear() {
    m_size = 0;

    return;
}

template <typename T>
void Vector<T>::push_back(const T& val) {
    // Container is not full but not empty
    if (m_size < m_max) { 
        ptr_to_data[m_size] = val;
        m_size++;
    }
    // Container is completely empty
    else if (ptr_to_data == nullptr) {
        ptr_to_data = new T[DEFAULT_MAX];
        ptr_to_data[0] = val;
        m_size++;
        m_max = DEFAULT_MAX;
    }
    // Container is full (m_size == m_max)
    else {
        // Create new container
        m_max *= 2;
        T* new_ptr_to_data = new T[m_max];

        // Copy data from old container to new container
        for (int i = 0; i < m_size; i++) {
            new_ptr_to_data[i] = ptr_to_data[i];
        }
        new_ptr_to_data[m_size] = val;
        m_size++;

        // Memory housekeeping
        delete [] ptr_to_data;
        ptr_to_data = new_ptr_to_data;
    }

    return;
}

template <typename T>
void Vector<T>::pop_back() {
    if (m_size <= 0) { throw out_of_range("Error: Attempted pop from empty vector."); }
    m_size--;
    
    return;
}

template <typename T>
T Vector<T>::abs_max_element() const {
    T result = abs(ptr_to_data[0]);
    for (int col = 0; col < m_size; col++) {
        result = (result < abs(ptr_to_data[col])) ? abs(ptr_to_data[col]) : result;
    }

    return result;
}

template <typename T>
ostream& operator<<(ostream& out, const Vector<T>& vec_out) {
    for (int i = 0; i < vec_out.get_size(); i++) {
        if (abs(vec_out[i]) < ZERO_LIMIT) { out << setprecision(PRECISION) << 0; }
        else { out << setprecision(PRECISION) << vec_out[i]; }
        if (i != vec_out.get_size() - 1) { out << endl; }
    }

    out << endl;

    return out;
}

template <typename T>
istream& operator>>(istream& in, Vector<T>& vec_in) {
    for (int i = 0; i < vec_in.get_max(); i++) {
        T temp;
        in >> temp;
        vec_in.push_back(temp);

        // Ensure that program halts if not given enough data
        if (i + 1 != vec_in.get_max() && in.peek() == '\n') {
            throw runtime_error("Error: Insufficient amount of data provided."); 
        }
    }

    return in;
}

template <typename T>
void swap(Vector<T>& v1, Vector<T>& v2) {
    std::swap(v1.m_size, v2.m_size);
    std::swap(v1.m_max, v2.m_max);
    std::swap(v1.ptr_to_data, v2.ptr_to_data);

    return;
}