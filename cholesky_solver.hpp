/*! \file
 *  Function definitions for the `Cholesky_Solver` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
Vector<T> Cholesky_Solver<T>::solve(const Base_Matrix<T>& matrix, Vector<T> vec) {
    m_size = matrix.get_size();

    // Storage for final lower triangular matrix produced from Cholesky decomposition
    L_Triangle_Matrix<T> l_matrix(matrix);

    // Used to not alter original matrix
    Vector<Vector<T>> matrix_data = matrix.get_elements();

    // Checking for division by zero
    T dividend;

    // Decompose
    for (int row = 0; row < m_size; row++) {
        // Compute values under the diagonal for the given row
        for (int col = 0; col < row; col++) {
            T sum(0);
            for (int runner = 0; runner < col; runner++) { 
                sum += l_matrix.get_element(row, runner) * l_matrix.get_element(col, runner);
            }

            dividend = l_matrix.get_element(col, col);
            if (dividend == 0) { throw domain_error("Error: Division by zero while solving symmetric matrix."); }

            T val = (1 / dividend) * (matrix_data[row][col] - sum);
            l_matrix.set_element(row, col, val);
        }

        // Compute value on the diagonal of this row
        T sum_squared_row(0);
        for (int col = 0; col < row; col++) {
            sum_squared_row += pow(l_matrix.get_element(row, col), 2);
        }

        dividend = matrix_data[row][row] - sum_squared_row;
        if (dividend < 0) { throw domain_error("Error: Imaginary numbers are about to run amok while solving a symmetric matrix."); }

        l_matrix.set_element(row, row, sqrt(dividend));
    }

    // Perform back substitutions for L and L*
    U_Triangle_Matrix<T> u_matrix(l_matrix.transpose());
    Vector<T> temp = l_matrix.back_sub(l_matrix.get_elements(), vec);

    return u_matrix.back_sub(u_matrix.get_elements(), temp);
}