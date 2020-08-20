//
// Created by AlTimofeyev on 5/16/2019.
//

#include "HarmonySearch.h"



// **********************************************************************************************************
// ************************************************ PUBLIC **************************************************
// **********************************************************************************************************
// ---------------------- CONSTRUCTORS ----------------------
/**
 * @brief The Harmony Search constructor.
 *
 * @note This is the only constructor for the Harmony Search, no default constructor exists.
 *
 * @param dimensions The number of elements per individual vector in the population.
 * @param populationSize The size of the population.
 * @param maxIterations The maximum number of iterations.
 * @param HMCR The Harmony Memory Consideration Rate.
 * @param PAR The Pitch Adjustment Rate.
 * @param bandwidth The bandwidth range.
 */
HarmonySearch::HarmonySearch(int dimensions, int populationSize, int maxIterations, double HMCR, double PAR, double bandwidth)
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

    // Set Harmony Search Configuration Parameters.
    hsConfig.dimensions = dimensions;
    hsConfig.popSize = populationSize;
    hsConfig.iterations = maxIterations;
    hsConfig.HMCR = HMCR;
    hsConfig.PAR = PAR;
    hsConfig.bandwidth = bandwidth;
}

// ------------------------- METHODS ------------------------
/**
 * @brief Runs the Harmony Search with set parameters.
 *
 * @param functionID The ID that references which Benchmark Function to use.
 * @param minBound, maxBound The minimum and maximum bounds of the individuals in Harmony Search.
 *
 * @return Returns the best global fitness.
 */
double HarmonySearch::runHarmonySearch(int functionID, double minBound, double maxBound)
{
    // Create a Mersenne Twister pseudo-random number generator.
    mt19937 randGenerator(time(NULL));

    // Record the start time.
    auto startTime = chrono::high_resolution_clock::now();

    // Create population with given parameters and evaluate it.
    HS_Population population;
    population.functionID = functionID;
    population.bounds.push_back(minBound);
    population.bounds.push_back(maxBound);
    generateHSPopulation(population, randGenerator);

    // Sort the population.
    quicksort(population.fitness, population.pop, 0, hsConfig.popSize-1);

    // Begin Harmony Search.
    int iteration = 0;
    while(iteration < hsConfig.iterations)
    {
        // Iterate the population.
        iterateHarmony(population, hsConfig.HMCR, hsConfig.PAR, hsConfig.bandwidth, randGenerator);

        // Save the global best fitness.
        population.bestGlobFit.push_back(population.fitness[0]);

        iteration++;    // Increment the iteration.
    }

    // Record the end time.
    auto endTime = chrono::high_resolution_clock::now();

    // Calculate elapsed time in milliseconds it took to run the Harmony Search.
    auto elapsedTime = endTime - startTime;
    double elapsedTimeMS = chrono::duration_cast<chrono::milliseconds>(elapsedTime).count();

    // Save elapsed time to the population.
    population.executionTime = elapsedTimeMS;

    // Add the population to the population list.
    popList.push_back(population);

    // Return the best fitness of the Harmony Search.
    int bestFitListSize = population.bestGlobFit.size();
    return population.bestGlobFit[bestFitListSize-1];
}

/**
 * @brief Analyzes the results of the Harmony Search.
 */
void HarmonySearch::analyzeHSResults()
{
    if(popList.empty())
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze  ***********\n";
        cout << "**********         Please Run The Harmony Search First           ***********";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new analysis object.
    hsAnalysis = HS_Analysis();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].bestGlobFit.size();

        // Save the function ID.
        hsAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].bestGlobFit);
        hsAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].bestGlobFit);
        hsAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].bestGlobFit[fitnessSize-1]);
        range.push_back(popList[i].bestGlobFit[0]);
        hsAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        hsAnalysis.medianFunctionFitness.push_back(popList[i].bestGlobFit[fitnessSize / 2]);

        // Save the execution time of data.
        hsAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        hsAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Analyses the worst results of the Harmony Search.
 */
void HarmonySearch::analyzeHSWorstResults()
{
    if(popList.empty())
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze  ***********\n";
        cout << "**********         Please Run The Harmony Search First           ***********";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new worst analysis object.
    hsWorstAnalysis = HS_Analysis_Worst();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].worstFitness.size();

        // Save the function ID.
        hsWorstAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].worstFitness);
        hsWorstAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].worstFitness);
        hsWorstAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].worstFitness[fitnessSize-1]);
        range.push_back(popList[i].worstFitness[0]);
        hsWorstAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        hsWorstAnalysis.medianFunctionFitness.push_back(popList[i].worstFitness[fitnessSize / 2]);

        // Save the execution time of data.
        hsWorstAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        hsWorstAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Prints the Results of the Harmony Search.
 */
void HarmonySearch::printHSResults()
{
    cout << "****************************************************************************\n";
    cout << "******* Printing Results of Harmony Search on Current Population ***********\n";
    cout << "----------------------------------------------------------------------------\n";

    // If the popList is empty, then Harmony Search has not been run yet.
    if(popList.empty())
    {
        cout << "********** NO RESULTS FOR THIS POPULATION\n";
        cout << "********** PLEASE RUN THE HARMONY SEARCH\n";
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
 * @brief Prints the Analysis of the Harmony Search Results.
 */
void HarmonySearch::printHSAnalysis()
{
    if(hsAnalysis.functionIDs.empty())
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
    for(int row = 0; row < hsAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << hsAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(hsAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << hsAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(hsAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << hsAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(hsAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << hsAnalysis.ranges[row][0] << "\t\t\t";
        if(hsAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << hsAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(hsAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << hsAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << hsAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << hsAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Prints the Analysis of the worst Harmony Search Results.
 */
void HarmonySearch::printHSWorstAnalysis()
{
    if(hsWorstAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "***************** There Is No Worst Analysis Data To Print *****************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    cout << "\n\n********************************************************\n";
    cout << "*********** Printing Worst Analysis Results ************\n";
    cout << "--------------------------------------------------------\n";

    cout << "Function ID\t\tAverage Fitness\t\t\tStandard Deviation\t\t\tRange(min)\t\t\tRange(max)\t\t\t\tMedian\t\t\t\tTime(ms)\t\t\tFunction Calls\n";
    cout.precision(12);
    for(int row = 0; row < hsWorstAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << hsWorstAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(hsWorstAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << hsWorstAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(hsWorstAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << hsWorstAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(hsWorstAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << hsWorstAnalysis.ranges[row][0] << "\t\t\t";
        if(hsWorstAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << hsWorstAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(hsWorstAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << hsWorstAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << hsWorstAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << hsWorstAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Saves all Harmony Search Results to file.
 */
void HarmonySearch::saveHSResults()
{
    // If the popList is empty, exit the function.
    if(popList.empty())
    {
        cout << "\n******************************************************\n";
        cout << "***** THERE IS NO HARMONY SEARCH DATA TO SAVE ********";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "HarmonySearch-Results.csv";

    // Initialize the number of rows (iterations/generations).
    int rows = hsConfig.iterations;

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
    if(hsAnalysis.functionIDs.size() > 0)
    {
        line = "Average,";
        for(int i = 0; i < hsAnalysis.avgFunctionFitness.size(); i++)
        {
            line += to_string(hsAnalysis.avgFunctionFitness[i]);

            if(i == hsAnalysis.avgFunctionFitness.size()-1)
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
 * @brief Saves the Analysis of the Harmony Search to file.
 */
void HarmonySearch::saveHSAnalysis()
{
    if(hsAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Save *********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = hsAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create filename.
    string filename = "HarmonySearch-Analysis.csv";


    // Create the file where to save the analysis.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    outputFile << hsAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(hsAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(hsAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(hsAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(hsAnalysis.ranges[row][0]) + ",";
        line += to_string(hsAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(hsAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(hsAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(hsAnalysis.functionCalls[row]) + "\n";

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Saves the Analysis of the worst Harmony Search Results to file.
 */
void HarmonySearch::saveHSWorstAnalysis()
{
    if(hsWorstAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "***************** There Is No Worst Analysis Data To Save ******************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = hsWorstAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create filename.
    string filename = "HarmonySearch-Analysis-Worst.csv";


    // Create the file where to save the analysis.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    outputFile << hsWorstAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(hsWorstAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(hsWorstAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(hsWorstAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(hsWorstAnalysis.ranges[row][0]) + ",";
        line += to_string(hsWorstAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(hsWorstAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(hsWorstAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(hsWorstAnalysis.functionCalls[row]) + "\n";

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
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
 * @param population The HS_Population structure that holds the population.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void HarmonySearch::generateHSPopulation(HS_Population &population, mt19937 &randGenerator)
{
    double minBound = population.bounds[0];
    double maxBound = population.bounds[1];

    // -------- Create Initial Population --------
    population.pop = createMatrixMT(hsConfig.popSize, hsConfig.dimensions, minBound, maxBound, randGenerator);

    // -------- Calculate the fitness --------
    population.fitness.resize(hsConfig.popSize);
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
void HarmonySearch::evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter)
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
void HarmonySearch::evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter)
{
    fitness = calculateFitnessOfVector(indiv, functionID);
    functionCounter++;
}

/**
 * @brief Choose a random harmonic (dimension) from the population.
 *
 * @param pop The matrix population.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 *
 * @return Returns a random harmonic (dimension) from the population.
 */
double HarmonySearch::chooseRandomHarmonic(const vector<vector<double>> &pop, mt19937 &randGenerator)
{
    // Create distributions for the Mersenne Twister pseudo-random number generator.
    uniform_int_distribution<int> solDis(0, pop.size()-1);      // For solution.
    uniform_int_distribution<int> dimDis(0, pop[0].size()-1);   // For harmonic (dimension).

    // Choose a random solution and harmonic index from the population.
    int randSolution = solDis(randGenerator);
    int randDimension = dimDis(randGenerator);

    // Return the harmonic.
    return pop[randSolution][randDimension];
}

/**
 * @brief Adjust the pitch of a random solution's harmonic (dimension).
 *
 * @param pop The matrix population.
 * @param dim The harmonic (dimension) where to adjust the pitch.
 * @param PAR The Pitch Adjustment Rate.
 * @param bandwidth The bandwidth range for the adjustment.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 *
 * @return Returns the harmonic with the adjusted pitch.
 */
double HarmonySearch::adjustHarmonicPitch(const vector<vector<double>> &pop, const int &dim, const double &PAR, const double &bandwidth, mt19937 &randGenerator)
{
    // Create distributions for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> pitchDis(0.0, 1.0);       // For random pitch.
    uniform_real_distribution<double> eDis(-1.0, 1.0);          // For epsilon.
    uniform_int_distribution<int> solDis(0, pop.size()-1);      // For solution.

    double randPitch = pitchDis(randGenerator); // Generate a random pitch
    double epsilon = eDis(randGenerator);       // Generate a random epsilon value.
    int randSolution = solDis(randGenerator);   // Choose a random solution from the population.

    double adjHarmonic; // Declare an adjusted harmonic.

    // Adjust the harmonic pitch.
    if(randPitch < PAR)
        adjHarmonic = pop[randSolution][dim] + bandwidth * epsilon;
    else
        adjHarmonic = pop[randSolution][dim] - bandwidth * epsilon;

    // Return the altered harmonic.
    return adjHarmonic;
}

/**
 * @brief Generate a new harmonic (dimension) within the bounds.
 *
 * @param minBound, maxBound The minimum and maximum bounds of the individual (harmonic).
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 *
 * @return Returns a newly generated harmonic within the specified bounds.
 */
double HarmonySearch::generateNewRandHarmonic(const double &minBound, const double &maxBound, mt19937 &randGenerator)
{
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);

    double randNum = dis(randGenerator);                            // Generate a random number.
    double newHarmonic = minBound + randNum*(maxBound - minBound);  // Generate a new harmonic within max/min bounds.

    // Return the new harmonic.
    return newHarmonic;
}

/**
 * @brief Generates a new random harmony from existing population.
 *
 * @param population The HS_Population structure that holds the population.
 * @param HMCR The Harmony Memory Consideration Rate.
 * @param PAR The Pitch Adjustment Rate.
 * @param bandwidth The bandwidth range.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 *
 * @return Returns a new random harmony (solution).
 */
vector<double> HarmonySearch::generateNewRandHarmony(HS_Population &population, const double &HMCR, const double &PAR, const double &bandwidth, mt19937 &randGenerator)
{
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Declare a new harmony and its fitness variables.
    vector<double> newHarmony;

    // Declare a temp harmonic and a random number variable.
    double harmonic;
    double randNum;

    //  Generate the new harmony.
    for(int dim = 0; dim < population.pop[0].size(); dim++)
    {
        randNum = dis(randGenerator);   // Generate a random value.

        // If less than HMCR, choose a random harmonic from the population.
        if(randNum < HMCR)
            harmonic = chooseRandomHarmonic(population.pop, randGenerator);

            // Else if it's less than PAR, pick a random solution with dimension 'dim' and adjust its pitch.
        else if(randNum < PAR)
            harmonic = adjustHarmonicPitch(population.pop, dim, PAR, bandwidth, randGenerator);

            // Else generate a new harmonic within bounds via randomization.
        else
            harmonic = generateNewRandHarmonic(population.bounds[0], population.bounds[1], randGenerator);

        // Add the harmonic to the new harmony.
        newHarmony.push_back(harmonic);
    }

    // Return the new harmony.
    return newHarmony;
}

/**
 * @brief Iterate the harmonic population.
 *
 * @param population The HS_Population structure that holds the population.
 * @param HMCR The Harmony Memory Consideration Rate.
 * @param PAR The Pitch Adjustment Rate.
 * @param bandwidth The bandwidth range.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void HarmonySearch::iterateHarmony(HS_Population &population, const double &HMCR, const double &PAR, const double &bandwidth, mt19937 &randGenerator)
{
    // Generate a new harmony from existing population.
    vector<double> newHarmony = generateNewRandHarmony(population, HMCR, PAR, bandwidth, randGenerator);

    // Evaluate the new harmony.
    double newHarmonyFitness;
    evaluateIndividual(population.functionID, newHarmony, newHarmonyFitness, population.functionCounter);

    int fitnessSize = population.fitness.size();    // Get the size of the fitness list.

    // Determine if the worst harmony in population can be replaced with new one.
    if(newHarmonyFitness < population.fitness[fitnessSize-1])
    {
        // Save worst harmony and its fitness.
        population.worstSol.push_back(population.pop[fitnessSize-1]);
        population.worstFitness.push_back(population.fitness[fitnessSize-1]);

        // Replace worst harmony/fitness with the new harmony/fitness.
        population.pop[fitnessSize-1] = newHarmony;
        population.fitness[fitnessSize-1] = newHarmonyFitness;

        // Sort the updated harmony population.
        quicksort(population.fitness, population.pop, 0, fitnessSize-1);
    }

    // Else save the new harmony and its fitness to the worst list.
    else
    {
        population.worstSol.push_back(newHarmony);
        population.worstFitness.push_back(newHarmonyFitness);
    }
}