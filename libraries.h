/*! \file
 *  Library includes, global constants, and function prototypes.
 */

//Programmers: Zachary Bahr and Jacob LeGrand

#ifndef LIBRARIES_H
#define LIBRARIES_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
using namespace std;

/*!  The default container maximum size for the Vector class
*/
const int DEFAULT_MAX = 1;

/*!  The precision of numeric data
*/
const int PRECISION = 8;

/*! The maximum number of iterations for the Gauss Seidel method
*/
const double ZERO_LIMIT = 0.000000001;

#endif