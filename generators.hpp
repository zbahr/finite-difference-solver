/*! \file
 *   Various matrix/vector generators
 */

//Programmers: Zachary Bahr and Jacob LeGrand

/*! Determines which boundary function to be called when generating the b vector
*
*  \param x the x value to be calculated
*  \param y the y value to be calculated
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param *upper a pointer to the upper boundary function
*  \param *lower a pointer to the lower boundary function
*  \param *right a pointer to the right boundary function
*  \param *left a pointer to the left boundary function
*  \return a value calculated from plugging x and y into the appropriate boundary function, or 0 if no boundary function is applicable
* 
*  \pre upper_bound > lower_bound
*  \post (see return)
*  \relates Matrix_Solver
*/
double boundary_func(const double x, const double y, const double lower_bound, const double upper_bound, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double)) {
    double result;

    if (x == lower_bound) result = left(x, y);
    else if (y == lower_bound) result = lower(x, y);
    else if (x == upper_bound) result = right(x, y);
    else if (y == upper_bound) result = upper(x, y);
    else result = 0;

    return result;
}

/*! Sum up the U (phi) values for all 4 points a distance of delta away from (x, y) in the mesh
*
*  \param x the x value to be calculated
*  \param y the y value to be calculated
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param delta the value of upper_bound/mesh_length
*  \param *upper a pointer to the upper boundary function
*  \param *lower a pointer to the lower boundary function
*  \param *right a pointer to the right boundary function
*  \param *left a pointer to the left boundary function
*  \return the sum of all the U (phi) values for all 4 points a distance away from (x, y) in the mesh
* 
*  \pre upper_bound > lower_bound
*  \pre T_function is defined for (double, double, double, double, double (double, double), double (double, double), double (double, double), double (double, double))
*  \post (see return)
*  \relates Matrix_Solver
*/
template <double T_Function(double, double, double, double, double (double, double), double (double, double), double (double, double), double (double, double))>
double callback(const double x, const double y, const double delta, const double lower_bound, const double upper_bound, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double)) {
    double sum = T_Function(x + delta, y, lower_bound, upper_bound, upper, lower, right, left) 
               + T_Function(x - delta, y, lower_bound, upper_bound, upper, lower, right, left) 
               + T_Function(x, y + delta, lower_bound, upper_bound, upper, lower, right, left) 
               + T_Function(x, y - delta, lower_bound, upper_bound, upper, lower, right, left);
    return sum;
}

/*! Generate the b vector to be used in Ax=b
*
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param mesh_length the length of the mesh
*  \param upper a pointer to the upper boundary function
*  \param lower a pointer to the lower boundary function
*  \param right a pointer to the right boundary function
*  \param left a pointer to the left boundary function
*  \return the generated b vector
* 
*  \pre upper_bound > lower_bound
*  \pre mesh_length > 0
*  \pre push_back() is defined for Vector<T>
*  \pre callback<boundary_function> returns a T value
*  \post (see return)
*  \relates Matrix_Solver
*/
template <typename T>
Vector<T> gen_callback_vec(const double lower_bound, const double upper_bound, const int mesh_length, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double)) {
    Vector<T> vec((mesh_length - 1) * (mesh_length - 1));
    double delta = upper_bound / mesh_length;

    /* The ZERO_LIMIT (defined in libraries.h) is used to differentiate between y's and x's that are very, 
       very close to the boundary. Without it, the vector is not generated with the correct size */
    for (double y = lower_bound + delta; y < (upper_bound - ZERO_LIMIT); y += delta) {
        for (double x = lower_bound + delta; x < (upper_bound - ZERO_LIMIT); x += delta) {
            vec.push_back(callback<boundary_func>(x, y, delta, lower_bound, upper_bound, upper, lower, right, left));
        }
    }

    return vec * 0.25;
}

/*! Generate the coefficient matrix A to be used in Ax=b
*
*  \param mesh_length the length of the mesh
*  \return a Symmetric_Matrix object that represents the coefficient A matrix
* 
*  \pre mesh_length > 0
*  \pre set_element() is defined for Symmetric_Matrix<T>
*  \pre callback<boundary_function> returns a T value
*  \post (see return)
*  \relates Matrix_Solver
*/
template <typename T>
Symmetric_Matrix<T> gen_coefficient_matrix(const int& mesh_length)
{
    int matrix_size = (mesh_length - 1) * (mesh_length - 1);
    Symmetric_Matrix<T> a(matrix_size, 0);

    // Set top left corner before starting loop
    a.set_element(0, 0, 1);

    for (int i = 1; i < matrix_size; i++)
    {
        // Set the diagonal to 1
        a.set_element(i, i, 1);

        // Alternate between -upperlimit/mesh_length and 0 every certain number of elements
        if (i % (mesh_length - 1) != 0)
            a.set_element(i, i - 1, -0.25);
        
        // Banded diagonals
        if (i >= mesh_length-1)
            a.set_element(i, i - (mesh_length - 1), -0.25);
    }

    return a;
}

/*! Write the x, y, and z values of our solution to a file, to be used with graphing
*
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param mesh_length the length of the mesh
*  \param result a Vector representing the solution vector x in Ax=b
*  \param file_name the name of the file to be written to
* 
*  \pre upper_bound > lower_bound
*  \pre mesh_length > 0
*  \pre result is not an empty Vector
*  \post the x, y, and z values of our soltuion are written to a file
*  \relates Matrix_Solver
*/
template <typename T>
void output_to_file(const double& lower_bound, const double& upper_bound, const int& mesh_length, const Vector<T>& result, const string& file_name) {
    ofstream fout;
    fout.open(file_name);

    // Output data;
    double delta = upper_bound / mesh_length;
    int i = 0;
    /* The ZERO_LIMIT (defined in libraries.h) is used to differentiate between y's and x's that are very, 
    very close to the boundary. Without it, the vector is not generated with the correct size */
    for (double y = lower_bound + delta; y < (upper_bound - ZERO_LIMIT); y += delta) {
        for (double x = lower_bound + delta; x < (upper_bound - ZERO_LIMIT); x += delta) {
            fout << x << "\t" << y << "\t" << result[i] << endl;
            i++;
        }
    }

    fout.close();

    return;
}

/*! Write the x, y, and z values of our solution to a file, to be used with graphing
*
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param mesh_length the length of the mesh
*  \param *exact_eqn a pointer to a function that returns the exact solution for the points we are approximating
*  \return a solution Vector that has been populated using the same x and y points as our approximation, but plugged into the exact equation 
*
*  \pre upper_bound > lower_bound
*  \pre mesh_length > 0
*  \pre push_back() is defined for Vector
*  \pre exact_eqn returns a T value
*  \post (see return)
*  \relates Matrix_Solver
*/
template <typename T>
Vector<T> gen_exact_sol(const double lower_bound, const double upper_bound, const int mesh_length, long double (*exact_eqn)(long double, long double)) {
    Vector<T> vec((mesh_length - 1) * (mesh_length - 1));
    double delta = upper_bound / mesh_length;

    /* The ZERO_LIMIT (defined in libraries.h) is used to differentiate between y's and x's that are very, 
       very close to the boundary. Without it, the vector is not generated with the correct size */
    for (double y = lower_bound + delta; y < (upper_bound - ZERO_LIMIT); y += delta) {
        for (double x = lower_bound + delta; x < (upper_bound - ZERO_LIMIT); x += delta) {
            vec.push_back(exact_eqn(x, y));
        }
    }

    return vec;
}