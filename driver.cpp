/*! driver.cpp
 * 
 * Desc:
 *      This file is essentially a driver to test the functionality of the
 *      Vector, Matrix_Solver, and matrix classes.
 */

// Programmers: Zachary Bahr and Jacob LeGrand

#include "u_triangle_matrix.h"
#include "l_triangle_matrix.h"
#include "symmetric_matrix.h"
#include "matrix_solver.h"
#include "speed_test.hpp"

/* Upper boundary function */
double upper(double x, double y)
{
    return ((0 * x) + (0 * y));
}

/* Lower boundary function */
double lower(double x, double y)
{
    return (sin(x) + (0 * y));
}

/* Right boundary function */
double right(double x, double y)
{
    return ((0 * x) + (0 * y));
}

/* Left boundary function */
double left(double x, double y)
{
    return ((0 * x) + sin(y));
}

/* Exact function used to test error */
long double exact_eqn(long double x, long double y) {
    return (1/sinh(M_PI)) * ((sin(x) * sinh(M_PI-y)) + (sin(y) * sinh(M_PI-x)));
}

int main() {

    try {
        double lower_bound = 0;
        double upper_bound = M_PI;
        int mesh_length;
        bool gauss_override;

        // Input mesh length
        cout << "Mesh length: "; 
        cin >> mesh_length;

        // Gaussian over-ride option
        cout << "Solver Technique (1 for Gaussian, 0 for Cholesky): ";
        cin >> gauss_override;

        // Solving method
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        Matrix_Solver<long double> solver(lower_bound, upper_bound, mesh_length, gauss_override, upper, lower, right, left);
        Vector<long double> result;

        // Solver test
        result = solver.solve();
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        double elapsed = duration_cast<microseconds>(t2-t1).count() / 1000000.0;
        cout << "Elapsed time: " << elapsed <<endl;

        // Find exact solution
        Vector<long double> exact_solution;
        exact_solution = gen_exact_sol<long double>(lower_bound, upper_bound, mesh_length, &exact_eqn);

        // Find norm
        long double norm;
        long double sum = 0;
        int size = result.get_size();
        for (int i=0; i<size; i++)
            sum += (result[i] - exact_solution[i]) * (result[i] - exact_solution[i]);
        norm = sqrt(((upper_bound/mesh_length)*(upper_bound/mesh_length)) * sum);
        cout << "Norm: " << norm <<endl;

        // Output data to file
        cout << "Outputting solution to file..." << endl;
        output_to_file(lower_bound, upper_bound, mesh_length, result, gauss_override ? "gauss_points_" + to_string(mesh_length) + ".txt" : "cholesky_points_" + to_string(mesh_length) + ".txt");
    }
    catch(const invalid_argument& err) { cerr << err.what() << endl; }
    catch(const runtime_error& err) { cerr << err.what() << endl; }
    catch(const domain_error& err) { cerr << err.what() << endl; }
    catch(const out_of_range& err) { cerr << err.what() << endl; }

    return 0;
}