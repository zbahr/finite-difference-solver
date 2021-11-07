/*! \file
 *   Testing the speed of gaussian and cholesky
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#include <chrono>
using namespace std::chrono;

/*! Test the speed of various mesh lengths using Gaussian elimination and Cholesky decomposition
*
*  \param lower_bound the lower bound of the mesh
*  \param upper_bound the upper bound of the mesh
*  \param max_mesh the maximum mesh length to be calculated and timed
*  \param *upper a pointer to the upper boundary function
*  \param *lower a pointer to the lower boundary function
*  \param *right a pointer to the right boundary function
*  \param *left a pointer to the left boundary function
* 
*  \pre upper_bound > lower_bound
*  \pre max_mesh > 0
*  \pre Matrix_Solver<long double> is defined
*  \post The speeds, in seconds, of Gaussian elimination and Cholesky decomposition for mesh lengths starting at 5 and increasing by 5 until reaching max_mesh are output to files
*  \relates Matrix_Solver
*/
void speed_test(const double lower_bound, const double upper_bound, const double max_mesh, double (*upper)(double, double), double (*lower)(double, double), double (*right)(double, double), double (*left)(double, double)) {
    ofstream fout;

    for (int i=0; i<2; i++)
    {
        fout.open(i ? "gaussian_times.txt" : "cholesky_times.txt");
        for (int mesh = 5; mesh <= max_mesh; mesh+=5)
        {
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            cout << "Calculating mesh length " << mesh << " using " << (i ? "Gaussian" : "Cholesky") <<endl;
            Matrix_Solver<long double> solver(lower_bound, upper_bound, mesh, i, upper, lower, right, left);
            Vector<long double> result;
            result = solver.solve();
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            double elapsed = duration_cast<microseconds>(t2-t1).count() / 1000000.0;
            fout << mesh << " " << elapsed << "\n";
        }
        fout.close();
    }
}