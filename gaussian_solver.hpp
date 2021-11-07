/*! \file
 *  Function definitions for the `Gaussian_Solver` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
Vector<T> Gaussian_Solver<T>::solve(const Base_Matrix<T>& matrix, Vector<T> vec) {
    m_size = matrix.get_size();
    m_matrix_data = matrix.get_elements();

    // Calculate scaling vector
    calculate_scales();

    // Perform Gaussian elimination with scaled partial pivoting
    for (int row_col = 0; row_col < m_size; row_col++) {
        //calculate_ratios(row_col);

        // Rearrange matrix so first row is the row with highest ratio
        //rearrange(row_col, vec);

        // Zero out all elements in column below diagonal row_col
        row_reduce(row_col, vec);

    }

    // Debug
    /*int j = 0;
    for (int i = 0; i < m_size; i++) {
        cout << m_matrix_data[i][j] << endl;
        j++;
    }*/

    return matrix.back_sub(m_matrix_data, vec);
}

template <typename T>
void Gaussian_Solver<T>::calculate_scales() {
    for (int row = 0; row < m_size; row++) {
        m_scales.push_back(m_matrix_data[row].abs_max_element());
    }

    return;
}

template <typename T>
void Gaussian_Solver<T>::calculate_ratios(const int& col) {
    if (m_scales.get_size() != m_size) { throw domain_error("Error: Scaling vector is not populated for gaussian elimination."); }

    m_ratios.clear();

    for (int row = 0; row < m_size; row++) {
        if (row < col) { m_ratios.push_back(0); }
        else {  
            if (m_matrix_data[row][col] == 0 && m_scales[row] == 0) {
                throw domain_error("Error: 0/0 division during ratio vector calculation.");
            }

            m_ratios.push_back(abs(m_matrix_data[row][col]) / m_scales[row]); 
        }
    }

    return;
}

template <typename T>
void Gaussian_Solver<T>::rearrange(const int& current_row, Vector<T>& vec) {
    if (m_ratios.get_size() != m_size) { throw domain_error("Error: Ratio vector is not populated for gaussian elimination."); }

    T highest_ratio = m_ratios.abs_max_element();
    int row_to_swap = 0;
    while (m_ratios[row_to_swap] != highest_ratio && row_to_swap < m_ratios.get_size()) { row_to_swap++; }

    if (row_to_swap != current_row) {
        std::swap(m_matrix_data[row_to_swap], m_matrix_data[current_row]);
        std::swap(vec[row_to_swap], vec[current_row]);
    }

    return;
}

template <typename T>
void Gaussian_Solver<T>::row_reduce(const int& row_col, Vector<T>& vec) {
    for (int runner = row_col + 1; runner < m_size; runner++) {
        if (m_matrix_data[runner][row_col] == 0) { continue; }
        T common_factor = m_matrix_data[row_col][row_col] / m_matrix_data[runner][row_col];
        m_matrix_data[runner] = m_matrix_data[runner] * (-common_factor) + m_matrix_data[row_col];
        vec[runner] = vec[runner] * (-common_factor) + vec[row_col];
    }

    return;
}