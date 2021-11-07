/*! \file
 *  Function definitions for the `Matrix_Solver` class.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

template <typename T>
Matrix_Solver<T>::Matrix_Solver(const Base_Matrix<T>& matrix, const Vector<T>& vec) : m_size(matrix.get_size()), m_matrix(&matrix), m_vec(vec), m_method(nullptr) {}

template <typename T>
Matrix_Solver<T>::Matrix_Solver(const double& lower_bound, const double& upper_bound, const int& mesh_length, const bool& gauss_override, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double)) {
    if (upper_bound <= lower_bound) { throw domain_error("Error: Upper bound should be greater than lower bound."); }
    if (mesh_length <= 1) { throw domain_error("Error: Mesh length should be greater than 1."); }
    
    m_size = (mesh_length - 1) * (mesh_length - 1);
    m_matrix = gauss_override ? (new General_Matrix<T>(gen_coefficient_matrix<T>(mesh_length))) : (new Symmetric_Matrix<T>(gen_coefficient_matrix<T>(mesh_length)));
    m_vec = gen_callback_vec<T>(lower_bound, upper_bound, mesh_length, upper, lower, right, left);
    m_method = nullptr;
}

template <typename T>
Vector<T> Matrix_Solver<T>::solve() {
    Vector<T> result;

    // Ensure size constraints are not violated
    if (m_size != m_vec.get_size()) { throw domain_error("Solving strategies must be performed on a dimensionally consistent matrix/vector pair."); }

    // If already row reduced ...
    if (m_matrix->get_status() == row_reduced) {
        result = m_matrix->back_sub(m_matrix->get_elements(), m_vec);
    } 
    // Perform Cholesky decomposition
    else if (m_matrix->get_status() == symmetric) {
        m_method = new Cholesky_Solver<T>();
        result = m_method->solve(*m_matrix, m_vec);
        delete m_method;
    }
    // Perform Gaussian elimination
    else {
        m_method = new Gaussian_Solver<T>();
        result = m_method->solve(*m_matrix, m_vec);
        delete m_method;
    }

    return result;
}