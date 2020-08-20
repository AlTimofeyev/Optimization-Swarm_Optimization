/**
 * @file FireflyAlgorithm.cpp
 * @author  Al Timofeyev
 * @date    May 14, 2019
 * @brief   This is an implementation of a Firefly Algorithm.
 */

#include "FireflyAlgorithm.h"


// **********************************************************************************************************
// ************************************************ PUBLIC **************************************************
// **********************************************************************************************************
// ---------------------- CONSTRUCTORS ----------------------
/**
 * @brief The Firefly Algorithm constructor.
 *
 * @note This is the only constructor for the Firefly Algorithm, no default constructor exists.
 *
 * @param dimensions The number of elements per individual vector in the population.
 * @param populationSize The size of the population.
 * @param maxIterations The maximum number of iterations.
 * @param alpha The alpha scaling factor.
 * @param beta The minimum beta scaling factor.
 * @param gamma The gamma scaling factor.
 */
FireflyAlgorithm::FireflyAlgorithm(int dimensions, int populationSize, int maxIterations, double alpha, double beta, double gamma)
{
    // Check the number of dimensions.
    if(dimensions < 1)
    {
        cout << "\n***************************************************************\n";
        cout << "******            Dimension Size Is Too Small            ******\n";
        cout << "****** Dimension size must be greater than or equal to 1 ******\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Check the population size.
    if(populationSize < 4)
    {
        cout << "\n***************************************************************\n";
        cout << "*****            Population Size Is Too Small            ******\n";
        cout << "***** Population size must be greater than or equal to 4 ******\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Check the number of iterations.
    if(maxIterations < 1)
    {
        cout << "\n***************************************************************\n";
        cout << "**             Number of Iterations Is Too Small           ****\n";
        cout << "** Number of Iterations must be greater than or equal to 1 ****\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Check the alpha scaling factor.
    if(alpha < 0 || alpha > 1)
    {
        cout << "\n***************************************************************\n";
        cout << "************ Alpha Scaling Factor Is Out Of Range *************\n";
        cout << "************     Alpha Must Be In Range [0,1]     *************\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Set Firefly Algorithm Configuration Parameters.
    faConfig.dimensions = dimensions;
    faConfig.popSize = populationSize;
    faConfig.iterations = maxIterations;
    faConfig.alpha = alpha;
    faConfig.betaMin = beta;
    faConfig.gamma = gamma;
}

// ------------------------- METHODS ------------------------
/**
 * @brief Runs the Firefly Algorithm with set parameters.
 *
 * @param functionID The ID that references which Benchmark Function to use.
 * @param minBound, maxBound The minimum and maximum bounds of the individuals in the firefly algorithm.
 *
 * @return Returns the best global fitness.
 */
double FireflyAlgorithm::runFireflyAlgorithm(int functionID, double minBound, double maxBound)
{
    // Create a Mersenne Twister pseudo-random number generator.
    mt19937 randGenerator(time(NULL));

    // Record the start time.
    auto startTime = chrono::high_resolution_clock::now();

    // Create population with given parameters and evaluate it.
    FA_Population population;
    population.functionID = functionID;
    population.bounds.push_back(minBound);
    population.bounds.push_back(maxBound);
    generateFAPopulation(population, randGenerator);

    // Iterate through the Firefly Algorithm.
    for(int iteration = 0; iteration < faConfig.iterations; iteration++)
    {
        // Iterate the fireflies.
        iterateFireflies(population, faConfig.alpha, faConfig.betaMin, faConfig.gamma, randGenerator);

        // Save the global best fitness.
        population.bestGlobFit.push_back(population.fitness[0]);
    }

    // Record the end time.
    auto endTime = chrono::high_resolution_clock::now();

    // Calculate elapsed time in milliseconds it took to run the Firefly Algorithm.
    auto elapsedTime = endTime - startTime;
    double elapsedTimeMS = chrono::duration_cast<chrono::milliseconds>(elapsedTime).count();

    // Save elapsed time to the population.
    population.executionTime = elapsedTimeMS;

    // Add the population to the population list.
    popList.push_back(population);

    // Return the best fitness of the Firefly Algorithm.
    int bestFitListSize = population.bestGlobFit.size();
    return population.bestGlobFit[bestFitListSize-1];
}

/**
 * @brief Analyzes the results of the Firefly Algorithm.
 */
void FireflyAlgorithm::analyzeFAResults()
{
    if(popList.empty())
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze  ***********\n";
        cout << "**********        Please Run The Firefly Algorithm First         ***********";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new analysis object.
    faAnalysis = FA_Analysis();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].bestGlobFit.size();

        // Save the function ID.
        faAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].bestGlobFit);
        faAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].bestGlobFit);
        faAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].bestGlobFit[fitnessSize-1]);
        range.push_back(popList[i].bestGlobFit[0]);
        faAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        faAnalysis.medianFunctionFitness.push_back(popList[i].bestGlobFit[fitnessSize / 2]);

        // Save the execution time of data.
        faAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        faAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Prints the Results of the Firefly Algorithm.
 */
void FireflyAlgorithm::printFAResults()
{
    cout << "****************************************************************************\n";
    cout << "****** Printing Results of Firefly Algorithm on Current Population *********\n";
    cout << "----------------------------------------------------------------------------\n";

    // If the popList is empty, then Firefly Algorithm has not been run yet.
    if(popList.empty())
    {
        cout << "********** NO RESULTS FOR THIS POPULATION\n";
        cout << "********** PLEASE RUN THE FIREFLY ALGORITHM\n";
        cout << "----------------------------------------------------------------------------\n\n";
        return;
    }

    for(int i = 0; i < popList.size(); i++)
    {
        cout << "Function ID: " << popList[i].functionID << endl;
        cout << "Generation\t\tBest Fitness of Generation\n";
        cout.precision(12);
        for (int gen = 0; gen < popList[i].bestGlobFit.size(); gen++)
            cout << gen << "\t\t\t\t" << popList[i].bestGlobFit[gen] << endl;

        cout << "\nElapsed Time (ms) : " << popList[i].executionTime << endl;
        cout << "----------------------------------------------------------------------------\n\n";
    }
}

/**
 * @brief Prints the Analysis of the Firefly Algorithm Results.
 */
void FireflyAlgorithm::printFAAnalysis()
{
    if(faAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Print ********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    cout << "\n\n********************************************************\n";
    cout << "************** Printing Analysis Results ***************\n";
    cout << "--------------------------------------------------------\n";

    cout << "Function ID\t\tAverage Fitness\t\t\tStandard Deviation\t\t\tRange(min)\t\t\tRange(max)\t\t\t\tMedian\t\t\t\tTime(ms)\t\t\tFunction Calls\n";
    cout.precision(12);
    for(int row = 0; row < faAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << faAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(faAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << faAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(faAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << faAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(faAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << faAnalysis.ranges[row][0] << "\t\t\t";
        if(faAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << faAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(faAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << faAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << faAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << faAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Saves all Firefly Algorithm Results to file.
 */
void FireflyAlgorithm::saveFAResults()
{
    // If the popList is empty, exit the function.
    if(popList.empty())
    {
        cout << "\n******************************************************\n";
        cout << "**** THERE IS NO FIREFLY ALGORITHM DATA TO SAVE ******";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "FireflyAlgorithm-Results.csv";

    // Initialize the number of rows (iterations/generations).
    int rows = faConfig.iterations;

    // Create the file to where the matrix is saved.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    string header = "Iteration,";
    for(int pIndex = 0; pIndex < popList.size(); pIndex++)
    {
        header += "f" + to_string(popList[pIndex].functionID);
        if(pIndex == popList.size()-1)
            header += "\n";
        else
            header += ",";
    }
    outputFile << header;

    // Save the data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the generation.
        line += to_string(row) + ",";

        // Save the best fitness from generation <row> of each population.
        for(int pIndex = 0; pIndex < popList.size(); pIndex++)
        {
            line += to_string(popList[pIndex].bestGlobFit[row]);

            if(pIndex == popList.size()-1)
                line += "\n";
            else
                line += ",";
        }

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Save the averages if they exist.
    if(faAnalysis.functionIDs.size() > 0)
    {
        line = "Average,";
        for(int i = 0; i < faAnalysis.avgFunctionFitness.size(); i++)
        {
            line += to_string(faAnalysis.avgFunctionFitness[i]);

            if(i == faAnalysis.avgFunctionFitness.size()-1)
                line += "\n";
            else
                line += ",";
        }
        outputFile << line;
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Saves the Analysis of the Firefly Algorithm to file.
 */
void FireflyAlgorithm::saveFAAnalysis()
{
    if(faAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Save *********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = faAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create filename.
    string filename = "FireflyAlgorithm-Analysis.csv";


    // Create the file where to save the analysis.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    outputFile << faAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(faAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(faAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(faAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(faAnalysis.ranges[row][0]) + ",";
        line += to_string(faAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(faAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(faAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(faAnalysis.functionCalls[row]) + "\n";

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Saves the ending population solutions to file.
 */
void FireflyAlgorithm::saveEndingPopulation()
{
    // If the popList is empty, exit the function.
    if(popList.empty())
    {
        cout << "\n******************************************************\n";
        cout << "******** THERE IS NO POPULATION DATA TO SAVE *********";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "FA-EndingPop-Function";
    string currentFilename;

    // Create a output file object.
    ofstream outputFile;

    // Save ending populations of all functions to file.
    string line = "";
    for(int i = 0; i < popList.size(); i++)
    {
        // Set the current filename.
        currentFilename = filename + to_string(popList[i].functionID) + ".csv";

        // Open/Create the file with the current filename.
        outputFile.open (currentFilename);

        // Save the population.
        for(int sol = 0; sol < popList[i].pop.size(); sol++)
        {
            // For all the dimensions in each solution vector.
            for(int dim = 0; dim < popList[i].pop[sol].size(); dim++)
            {
                line += to_string(popList[i].pop[sol][dim]);
                if(dim == popList[i].pop[sol].size()-1)
                    line += "\n";
                else
                    line += ",";
            }

            // Save the line to file.
            outputFile << line;

            // Reset the line.
            line = "";
        }

        // Close the file.
        outputFile.close();
    }
}


// **********************************************************************************************************
// ************************************************ PRIVATE *************************************************
// **********************************************************************************************************
// ------------------------- METHODS ------------------------
/**
 * @brief Generates the initial population.
 *
 * @note Makes function call to SwarmUtilities.h --> createMatrixMT().
 *
 * @param population The FA_Population structure that holds the population.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void FireflyAlgorithm::generateFAPopulation(FA_Population &population, mt19937 &randGenerator)
{
    double minBound = population.bounds[0];
    double maxBound = population.bounds[1];

    // -------- Create Initial Population --------
    population.pop = createMatrixMT(faConfig.popSize, faConfig.dimensions, minBound, maxBound, randGenerator);

    // -------- Calculate the fitness --------
    population.fitness.resize(faConfig.popSize);
    evaluatePopulation(population.functionID, population.pop, population.fitness, population.functionCounter);
}

/**
 * @brief Calculates fitness of all solutions in population.
 *
 * @note Makes function call to SwarmUtilities.h --> calculateFitnessOfVector().
 *
 * @param functionID The ID of the benchmark function to use.
 * @param pop The population matrix.
 * @param fitness The fitness vector for each solution from the population.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void FireflyAlgorithm::evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter)
{
    for(int i = 0; i < pop.size(); i++)
    {
        fitness[i] = calculateFitnessOfVector(pop[i], functionID);
        functionCounter++;
    }
}

/**
 * @brief Calculate the fitness of an individual solution of the population.
 *
 * @note Makes function call to SwarmUtilities.h --> calculateFitnessOfVector().
 *
 * @param functionID The ID of the benchmark function to use.
 * @param indiv The individual of the population.
 * @param fitness The fitness variable for the individual.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void FireflyAlgorithm::evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter)
{
    fitness = calculateFitnessOfVector(indiv, functionID);
    functionCounter++;
}

/**
 * @brief Calculates the distance between two fireflies.
 *
 * @param firefly1 The first firefly.
 * @param firefly2 The second firefly.
 *
 * @return Returns the distance between the two fireflies.
 */
double FireflyAlgorithm::calculateDistanceBetweenFireflies(const vector<double> &firefly1, const vector<double> &firefly2)
{
    // Declare the distance variable.
    double r;

    // Initialize the summation variable to 0.
    double summedUp = 0;

    // Sum up the dimensions.
    for(int dim = 0; dim < firefly1.size(); dim++)
        summedUp += pow((firefly1[dim] - firefly2[dim]), 2);

    // Get the distance using the summed up values.
    r = sqrt(summedUp);

    // Return the distance.
    return r;
}

/**
 * @brief Calculates the light intensity of two fireflies.
 *
 * @note lIntensity1/lIntensity2 can be null, initialized, or uninitialized since they'll be changed in the function.
 *
 * @param ffFit1 Fitness of the first firefly.
 * @param ffFit2 Fitness of the second firefly.
 * @param r The distance between the two fireflies.
 * @param gamma The gamma scaling factor.
 * @param lIntensity1 The light intensity of the first firefly.
 * @param lIntensity2 The light intensity of the second firefly.
 */
void FireflyAlgorithm::calculateLightIntensity(const double &ffFit1, const double &ffFit2, const double &r, const double &gamma, double &lIntensity1, double &lIntensity2)
{
    // Calculate light intensity of the first firefly.
    lIntensity1 = ffFit1 * exp(-(gamma) * pow(r, 2));

    // Calculate light intensity of the second firefly.
    lIntensity2 = ffFit2 * exp(-(gamma) * pow(r, 2));
}

/**
 * @brief Calculates the attractiveness between two fireflies.
 *
 * @param betaMin The beta scaling factor.
 * @param gamma The gamma scaling factor.
 * @param r The distance between the two fireflies.
 *
 * @return Returns the attractiveness value between two fireflies.
 */
double FireflyAlgorithm::calculateAttractiveness(const double &betaMin, const double &gamma, const double &r)
{
    // Calculate the attractiveness.
    double beta = betaMin * exp(-(gamma) * pow(r, 2));

    // Return the attractiveness.
    return beta;
}

/**
 * @brief Moves the less brighter firefly towards the brighter firefly.
 *
 * @param brightFF The brighter firefly.
 * @param lessBrightFF The less brighter firefly.
 * @param alpha The alpha scaling factor.
 * @param betaMin The beta scaling factor.
 * @param gamma The gamma scaling factor.
 * @param r The distance between the two fireflies.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void FireflyAlgorithm::moveFirefly(const vector<double> &brightFF, vector<double> &lessBrightFF, const double &alpha, const double &betaMin, const double &gamma, const double &r, mt19937 &randGenerator)
{
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);
    double randNum;

    // Get the attractiveness between the two fireflies.
    double beta = calculateAttractiveness(betaMin, gamma, r);

    // Move the less brighter firefly towards the brighter firefly.
    for(int dim = 0; dim < lessBrightFF.size(); dim++)
    {
        randNum = dis(randGenerator);   // Generate random number that will be used to determine the Gaussian random number below.
        lessBrightFF[dim] = lessBrightFF[dim] + beta * (brightFF[dim] - lessBrightFF[dim]) + alpha * (randNum - 0.5);
    }
}

/**
 * @brief Moves all fireflies that are less brighter than the current firefly, towards the current firefly.
 *
 * @param currFFIndex Index of the current firefly in the population.
 * @param population The population of fireflies.
 * @param alpha The alpha scaling factor.
 * @param betaMin The beta scaling factor.
 * @param gamma The gamma scaling factor.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void FireflyAlgorithm::moveLessBrightFireflies(const int &currFFIndex, FA_Population &population, const double &alpha, const double &betaMin, const double &gamma, mt19937 &randGenerator)
{
    // Declare the distance and light intensity variables.
    // The current firefly will always assume lightIntensity1.
    double r, lightIntensity1, lightIntensity2;

    // Loop through the population of fireflies.
    for(int ffIndex = 0; ffIndex < population.pop.size(); ffIndex++)
    {
        // Calculate distance between the current firefly and the firefly at index ffIndex.
        r = calculateDistanceBetweenFireflies(population.pop[currFFIndex], population.pop[ffIndex]);

        // Calculate light intensity of the two fireflies.
        calculateLightIntensity(population.fitness[currFFIndex], population.fitness[ffIndex], r, gamma, lightIntensity1, lightIntensity2);

        // If the current firefly is brighter than the firefly at index ffIndex.
        if(lightIntensity2 < lightIntensity1)
        {
            // Move the firefly at index ffIndex towards the current firefly.
            moveFirefly(population.pop[currFFIndex], population.pop[ffIndex], alpha, betaMin, gamma, r, randGenerator);

            // Update the fitness of the firefly at index ffIndex.
            evaluateIndividual(population.functionID, population.pop[ffIndex], population.fitness[ffIndex], population.functionCounter);
        }
    }
}

/**
 * @brief Iterates the firefly population to the next generation.
 *
 * @param population The population of fireflies.
 * @param alpha The alpha scaling factor.
 * @param betaMin The beta scaling factor.
 * @param gamma The gamma scaling factor.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void FireflyAlgorithm::iterateFireflies(FA_Population &population, const double &alpha, const double &betaMin, const double &gamma, mt19937 &randGenerator)
{
    // Iterate the fireflies to the next generation.
    for(int currFFIndex = 0; currFFIndex < population.pop.size(); currFFIndex++)
        moveLessBrightFireflies(currFFIndex, population, alpha, betaMin, gamma, randGenerator);

    // Sort the new population of fireflies.
    quicksort(population.fitness, population.pop, 0, population.fitness.size()-1);
}