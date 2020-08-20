/**
 * @file BenchmarkFunctions.h
 * @author  Al Timofeyev
 * @date    April 17, 2019
 * @brief   A library of benchmark functions.
 */

#ifndef BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H
#define BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H



#define _USE_MATH_DEFINES // Uncomment if cmath constants are desirable, like M_PI.

#include <vector>
#include <math.h>   // Sine and Cosine, square root
#include <cmath>    // Absolute value of doubles

using namespace std;

/** Performs the Schefel's Function on a vector of elements. */
double schefelsFunc(vector<double> &vect, int size);

/** Performs the 1st De Jong’s Function on a vector of elements. */
double deJongsFunc(vector<double> &vect, int size);

/** Performs the Rosenbrock Function on a vector of elements. */
double rosenbrockFunc(vector<double> &vect, int size);

/** Performs the Rastrigin Function on a vector of elements. */
double rastriginFunc(vector<double> &vect, int size);

/** Performs the Griewangk Function on a vector of elements. */
double griewangkFunc(vector<double> &vect, int size);

/** Performs the Sine Envelope Sine Wave Function on a vector of elements. */
double sineEnvelopeSineWaveFunc(vector<double> &vect, int size);

/** Performs the Stretched V Sine Wave Function on a vector of elements. */
double stretchedVSineWaveFunc(vector<double> &vect, int size);

/** Performs the Ackley’s One Function on a vector of elements. */
double ackleysOneFunc(vector<double> &vect, int size);

/** Performs the Ackley’s Two Function on a vector of elements. */
double ackleysTwoFunc(vector<double> &vect, int size);

/** Performs the Egg Holder Function on a vector of elements. */
double eggHolderFunc(vector<double> &vect, int size);

/** Performs the Rana Function on a vector of elements. */
double ranaFunc(vector<double> &vect, int size);

/** Performs the Pathological Function on a vector of elements. */
double pathologicalFunc(vector<double> &vect, int size);

/** Performs the Michalewicz Function on a vector of elements. */
double michalewiczFunc(vector<double> &vect, int size);

/** Performs the Masters Cosine Wave Function on a vector of elements. */
double mastersCosWaveFunc(vector<double> &vect, int size);

/** Performs the Quartic Function on a vector of elements. */
double quarticFunc(vector<double> &vect, int size);

/** Performs the Levy Function on a vector of elements. */
double levyFunc(vector<double> &vect, int size);

/** Performs the Step Function on a vector of elements. */
double stepFunc(vector<double> &vect, int size);

/** Performs the Alpine Function on a vector of elements. */
double alpineFunc(vector<double> &vect, int size);

#endif //BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H