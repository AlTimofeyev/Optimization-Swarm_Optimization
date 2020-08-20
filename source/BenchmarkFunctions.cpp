/**
 * @file BenchmarkFunctions.cpp
 * @author  Al Timofeyev
 * @date    April 17, 2019
 * @brief   A library of benchmark functions.
 */

#include "BenchmarkFunctions.h"

using namespace std;

// **********************************************************************************
// **************************** Benchmark Functions Below ***************************
// **********************************************************************************
/**
 * @brief Performs the Schefel's Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double schefelsFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        summedUp += (-vect[i]) * sin(sqrt(abs(vect[i])));

    answer = (418.9829 * size) - summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the 1st De Jong’s Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double deJongsFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        answer += pow(vect[i], 2);

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Rosenbrock Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double rosenbrockFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += 100 * pow((pow(vect[i], 2) - vect[i+1]), 2) + pow((1-vect[i]), 2);

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Rastrigin Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double rastriginFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        summedUp += pow(vect[i], 2) - (10 * cos(2*M_PI*vect[i]));

    answer = 10 * size * summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Griewangk Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double griewangkFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;
    double productUp = 1;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        summedUp += pow(vect[i], 2) / 4000;

    for(int i = 0; i < size; ++i)
        productUp *= cos(vect[i] / sqrt(i+1));

    answer = 1 + summedUp - productUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Sine Envelope Sine Wave Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double sineEnvelopeSineWaveFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        summedUp += 0.5 + pow(sin(pow(vect[i], 2) + pow(vect[i+1], 2) - 0.5), 2) / pow((1 + 0.001*(pow(vect[i], 2) + pow(vect[i+1], 2))), 2);

    answer = -summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Stretched V Sine Wave Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double stretchedVSineWaveFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += pow((pow(vect[i], 2) + pow(vect[i+1], 2)), 1.0/4) * pow(sin(50 * pow((pow(vect[i], 2) + pow(vect[i+1], 2)), 1.0/10)), 2) + 1;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Ackley’s One Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double ackleysOneFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += (1 / pow(exp(1.0), 0.2)) * sqrt(pow(vect[i] ,2) + pow(vect[i+1], 2)) + 3 * (cos(2*vect[i]) + sin(2*vect[i+1]));

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Ackley’s Two Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double ackleysTwoFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += 20 + exp(1.0) - (20 / pow(exp(1.0), (0.2 * sqrt((pow(vect[i], 2) + pow(vect[i+1], 2)) / 2))))
                - pow(exp(1.0), (0.5 * (cos(2*M_PI*vect[i]) + cos(2*M_PI*vect[i+1]))));

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Egg Holder Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double eggHolderFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += (-vect[i]) * sin(sqrt(abs(vect[i] - vect[i+1] - 47)))
                - (vect[i+1] + 47) * sin(sqrt(abs(vect[i+1] + 47 + (vect[i] / 2))));

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Rana Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double ranaFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += vect[i] * sin(sqrt(abs(vect[i+1] - vect[i] + 1))) * cos(sqrt(abs(vect[i+1] + vect[i] + 1)))
                + (vect[i+1] + 1) * cos(sqrt(abs(vect[i+1] - vect[i] + 1))) * sin(sqrt(abs(vect[i+1] + vect[i] + 1)));

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Pathological Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double pathologicalFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        answer += 0.5 + (pow(sin(sqrt(100 * pow(vect[i], 2) + pow(vect[i+1], 2))), 2) - 0.5)
                / (1 + 0.001*pow((pow(vect[i], 2) - 2*vect[i] * vect[i+1] + pow(vect[i+1], 2)), 2));

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Michalewicz Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double michalewiczFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        summedUp += sin(vect[i]) * pow(sin(((i+1) * pow(vect[i], 2)) / M_PI), 20);

    answer = -summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Masters Cosine Wave Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double mastersCosWaveFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
        summedUp += pow(exp(1.0), ((-1.0/8.0) * (pow(vect[i], 2) + pow(vect[i+1], 2) + 0.5 * vect[i+1] * vect[i])))
                * cos(4 * sqrt(pow(vect[i], 2) + pow(vect[i+1], 2) + 0.5 * vect[i] * vect[i+1]));

    answer = -summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Quartic Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double quarticFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        answer += (i+1) * pow(vect[i], 4);

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Levy Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double levyFunc(vector<double> &vect, int size)
{
    double answer;
    double summedUp = 0;

    double w1 = vect[0];
    double wi;
    double wn = vect[size-1];

    // perform function calculations.
    for(int i = 0; i < size-1; ++i)
    {
        wi = 1 + (vect[i] - 1) / 4;
        summedUp += pow((wi - 1), 2) * (1 + 10 * pow(sin(M_PI*wi + 1), 2))
                + pow((wn - 1), 2) * (1 + pow(sin(2*M_PI*wn), 2));
    }

    answer = pow(sin(M_PI * w1), 2) + summedUp;

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Step Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double stepFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        answer += pow((abs(vect[i]) + 0.5), 2);

    return answer;
}

// ----------------------------------------------------------------------------------
/**
 * @brief Performs the Alpine Function on a vector of elements.
 *
 * @param vect The vector of elements on which to perform calculations.
 * @param size The number of elements in vector.
 * @return The results of the calculations (fitness).
 */
double alpineFunc(vector<double> &vect, int size)
{
    double answer = 0;

    // perform function calculations.
    for(int i = 0; i < size; ++i)
        answer += abs(vect[i] * sin(vect[i]) + 0.1 * vect[i]);

    return answer;
}