/**
 * @file SwarmUtilities.cpp
 * @author  Al Timofeyev
 * @date    May 10, 2019
 * @brief   Utilities library for swarm algorithms.
 */

#include "SwarmUtilities.h"


/**
 * @brief Prints all the possible Function IDs to the screen.
 *
 * Prints all possible Function ID, as well as the functions they
 * reference, to the screen.
 */
void printAllFunctionIDs()
{
    cout << "\n********************************************************\n";
    cout << "All Possible Function IDs and Their Respective Functions";
    cout << "\n--------------------------------------------------------\n";
    cout << "Function ID: 1\tFunction Name: Schwefels function\n";
    cout << "Function ID: 2\tFunction Name: 1st De Jongs function\n";
    cout << "Function ID: 3\tFunction Name: Rosenbrock function\n";
    cout << "Function ID: 4\tFunction Name: Rastrigin function\n";
    cout << "Function ID: 5\tFunction Name: Griewangk function\n";
    cout << "Function ID: 6\tFunction Name: Sine Envelope Sine Wave function\n";
    cout << "Function ID: 7\tFunction Name: Stretched V Sine Wave function\n";
    cout << "Function ID: 8\tFunction Name: Ackleys One function\n";
    cout << "Function ID: 9\tFunction Name: Ackleys Two function\n";
    cout << "Function ID: 10\tFunction Name: Egg Holder function\n";
    cout << "Function ID: 11\tFunction Name: Rana function\n";
    cout << "Function ID: 12\tFunction Name: Pathological function\n";
    cout << "Function ID: 13\tFunction Name: Michalewicz function\n";
    cout << "Function ID: 14\tFunction Name: Masters Cosine Wave function\n";
    cout << "Function ID: 15\tFunction Name: Quartic function\n";
    cout << "Function ID: 16\tFunction Name: Levy function\n";
    cout << "Function ID: 17\tFunction Name: Step function\n";
    cout << "Function ID: 18\tFunction Name: Alpine function\n";
    cout << "********************************************************\n\n";
}

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------

/**
 * @brief Creates a matrix of doubles using Mersenne Twister.
 *
 * A matrix is constructed using the Mersenne Twister in the <random> library
 * with the user-specified min/max boundaries.
 *
 * @param rows The number of vectors in the matrix.
 * @param columns The number of elements in each vector of the matrix.
 * @param minBound, maxBound The max/min boundaries are the range
 *                           in which to generate numbers.
 *
 * @return The fully constructed matrix of doubles.
 */
vector<vector<double>> createMatrix(int rows, int columns, double minBound, double maxBound)
{
    // Create a Mersenne Twister pseudo-random number generator and a distribution.
    mt19937 randGenerator(time(NULL));
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Declare the matrix
    vector<vector<double>> matrix;

    // Declare the necessary variables to hold temporary numbers.
    double num, randNum;

    // Create all the rows for the matrix.
    for(int row = 0; row < rows; row++)
    {
        vector<double> vectorOfDoubles;

        // Generate all the elements of the vector.
        for (int col = 0; col < columns; col++)
        {
            // Generate a random number using Mersenne Twister
            randNum = dis(randGenerator);

            // Normalize the random number to the bounds.
            num = minBound + randNum*(maxBound - minBound);

            // Add value to vector.
            vectorOfDoubles.push_back(num);
        }

        // Add the vector of doubles to the matrix.
        matrix.push_back(vectorOfDoubles);
    }

    // Return the matrix.
    return matrix;
}

/**
 * @brief Creates a matrix of doubles using Mersenne Twister.
 *
 * A matrix is constructed using the Mersenne Twister in the <random> library
 * with the user-specified min/max boundaries.
 *
 * @param rows The number of vectors in the matrix.
 * @param columns The number of elements in each vector of the matrix.
 * @param minBound, maxBound The max/min boundaries are the range
 *                           in which to generate numbers.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 *
 * @return The fully constructed matrix of doubles.
 */
vector<vector<double>> createMatrixMT(int rows, int columns, double minBound, double maxBound, mt19937 &randGenerator)
{
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Declare the matrix
    vector<vector<double>> matrix;

    // Declare the necessary variables to hold temporary numbers.
    double num, randNum;

    // Create all the rows for the matrix.
    for(int row = 0; row < rows; row++)
    {
        vector<double> vectorOfDoubles;

        // Generate all the elements of the vector.
        for (int col = 0; col < columns; col++)
        {
            // Generate a random number using Mersenne Twister
            randNum = dis(randGenerator);

            // Normalize the random number to the bounds.
            num = minBound + randNum*(maxBound - minBound);

            // Add value to vector.
            vectorOfDoubles.push_back(num);
        }

        // Add the vector of doubles to the matrix.
        matrix.push_back(vectorOfDoubles);
    }

    // Return the matrix.
    return matrix;
}

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------

/**
 * @brief Calculates the fitness of a vector.
 *
 * The fitness of a vector is calculated by the Benchmark Function
 * referenced by the functionID.
 *
 * @note This function makes a call to BenchmarkFunctions.h.
 *
 * @param vect The vector of elements on which the Benchmark Functions operate.
 * @param functionID The ID that references which Benchmark Function to use.
 *
 * @return The fitness of the vector.
 */
double calculateFitnessOfVector(vector<double> &vect, int functionID)
{
    switch(functionID)
    {
        case 1:
            return schefelsFunc(vect, vect.size());
        case 2:
            return deJongsFunc(vect, vect.size());
        case 3:
            return rosenbrockFunc(vect, vect.size());
        case 4:
            return rastriginFunc(vect, vect.size());
        case 5:
            return griewangkFunc(vect, vect.size());
        case 6:
            return sineEnvelopeSineWaveFunc(vect, vect.size());
        case 7:
            return stretchedVSineWaveFunc(vect, vect.size());
        case 8:
            return ackleysOneFunc(vect, vect.size());
        case 9:
            return ackleysTwoFunc(vect, vect.size());
        case 10:
            return eggHolderFunc(vect, vect.size());
        case 11:
            return ranaFunc(vect, vect.size());
        case 12:
            return pathologicalFunc(vect, vect.size());
        case 13:
            return michalewiczFunc(vect, vect.size());
        case 14:
            return mastersCosWaveFunc(vect, vect.size());
        case 15:
            return quarticFunc(vect, vect.size());
        case 16:
            return levyFunc(vect, vect.size());
        case 17:
            return stepFunc(vect, vect.size());
        case 18:
            return alpineFunc(vect, vect.size());

        default:
            cout << "\n********************************************************\n";
            cout << "Fitness Process Failed for Function ID: " << functionID << endl;
            printAllFunctionIDs();
            cout << "************ TERMINATING PROGRAM EXECUTION *************\n\n";
            exit(1);
    }
}

/**
 * @brief Calculates the fitness of all vectors of a matrix.
 *
 * Calculates the fitness of all the vectors of the matrix stored
 * All the fitness results are stored in the fitness vector variable.
 *
 * @param matrix The matrix that holds all the vectors for calculating the fitness.
 * @param functionID The ID of the function to use for calculating the fitness.
 *
 * @return A vector of fitness values.
 */
vector<double> calculateFitnessOfMatrix(vector<vector<double>> matrix, int functionID)
{
    // Variables to hold the fitness of each vector.
    double fitness;
    vector<double> fitnessList;

    // Calculate the fitness of all rows in matrix.
    for(int row = 0; row < matrix.size(); row++)
    {
        fitness = calculateFitnessOfVector(matrix[row], functionID);
        fitnessList.push_back(fitness);
    }

    // Return the fitness list.
    return fitnessList;
}

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------

/**
 * @brief Calculates the average value of a vector of doubles.
 * @param vect The vector of doubles.
 * @return The average value of the vector.
 */
double calculateAverage(vector<double> vect)
{
    double average;
    double summedUp = 0;

    // Sum up all the fitness values.
    for(int row = 0; row < vect.size(); row++)
        summedUp += vect[row];

    // Calculate average.
    average = summedUp / vect.size();

    // Return the average.
    return average;
}

/**
 * @brief Calculates the standard deviation value of a vector of doubles.
 * @param vect The vector of doubles.
 * @return The standard deviation value of the vector.
 */
double calculateStandardDeviation(vector<double> vect)
{
    double stdDeviation;
    double summedUp = 0;
    double average = calculateAverage(vect);
    int size = vect.size();

    for(int row = 0; row < size; row++)
        summedUp += pow((vect[row] - average), 2);

    stdDeviation = sqrt((1.0/size) * summedUp);

    return stdDeviation;
}

//*******************************************************************************************************
//*******************************************************************************************************

/**
 * @brief Sorts a matrix and its fitness vector based on the fitness.
 *
 * @note Sorted in Ascending Order.
 * @note Smallest (minimum) fitness gets moved to index 0, along with its vector from matrix.
 * @note Largest (maximum) fitness gets moved to the last index, along with its vector from matrix.
 *
 * @param fitnessList The list of fitness values that correspond to each row of the matrix.
 * @param matrix A matrix of double values.
 * @param L The starting index for the quicksort (inclusive).
 * @param R The ending index for the quicksort (inclusive).
 */
void quicksort(vector<double> &fitnessList, vector<vector<double>> &matrix, int L, int R)
{
    int i, j, mid;
    double piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = fitnessList[mid];

    while (i<R || j>L)
    {
        while (fitnessList[i] < piv)
            i++;
        while (fitnessList[j] > piv)
            j--;

        if (i <= j)
        {
            swap(fitnessList, matrix, i, j);
            i++;
            j--;
        }
        else
        {
            if (i < R)
                quicksort(fitnessList, matrix, i, R);
            if (j > L)
                quicksort(fitnessList, matrix, L, j);
            return;
        }
    }
}

/**
 * @brief Swaps the fitness' and their corresponding vectors in the matrix.
 *
 * @param fitnessList The list of fitness values that correspond to each row of the matrix.
 * @param matrix A matrix of double values.
 * @param x The 1st index of the fitness/vector for the swap.
 * @param y The 2nd index of the fitness/vector for the swap.
 */
void swap(vector<double> &fitnessList, vector<vector<double>> &matrix, int x, int y)
{
    // Swap fitness values.
    double fitTemp = fitnessList[x];
    fitnessList[x] = fitnessList[y];
    fitnessList[y] = fitTemp;

    // Swap vector values.
    vector<double> vectTemp = matrix[x];
    matrix[x] = matrix[y];
    matrix[y] = vectTemp;
}

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------

/**
 * @brief A normal Quicksort implementation for vector arrays of doubles.
 *
 * @note Sorted in Ascending Order.
 * @note Smallest value gets moved to index 0.
 * @note Largest value gets moved to the last index.
 *
 * @param vec Vector array of doubles.
 * @param L The starting index for the quicksort (inclusive).
 * @param R The ending index for the quicksort (inclusive).
 */
void quicksort(vector<double> &vec, int L, int R) {
    int i, j, mid;
    double piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = vec[mid];

    while (i<R || j>L)
    {
        while (vec[i] < piv)
            i++;
        while (vec[j] > piv)
            j--;

        if (i <= j)
        {
            swap(vec, i, j); //error=swap function doesnt take 3 arguments
            i++;
            j--;
        }
        else
        {
            if (i < R)
                quicksort(vec, i, R);
            if (j > L)
                quicksort(vec, L, j);
            return;
        }
    }
}

/**
 * @brief Swaps two values of a vector array of doubles.
 *
 * @param v The vector in which values are swapped.
 * @param x The 1st index of the fitness/vector for the swap.
 * @param y The 2nd index of the fitness/vector for the swap.
 */
void swap(vector<double> &v, int x, int y)
{
    double temp = v[x];
    v[x] = v[y];
    v[y] = temp;
}