/**
 * @file ParticleSwarm.cpp
 * @author  Al Timofeyev
 * @date    May 10, 2019
 * @brief   This is an implementation of a Particle Swarm.
 */

#include "ParticleSwarm.h"


// **********************************************************************************************************
// ************************************************ PUBLIC **************************************************
// **********************************************************************************************************
// ---------------------- CONSTRUCTORS ----------------------
/**
 * @brief The Particle Swarm constructor.
 *
 * @note This is the only constructor for the Particle Swarm, no default constructor exists.
 *
 * @param dimensions The number of elements per individual vector in the population.
 * @param populationSize The size of the population.
 * @param maxIterations The maximum number of iterations.
 * @param kDampeningFactor The dampening factor, to dampen the velocities of particles.
 * @param c1 The scaling factor to bring velocity closer to personal best.
 * @param c2 The scaling factor to bring velocity closer to global best.
 */
ParticleSwarm::ParticleSwarm(int dimensions, int populationSize, int maxIterations, double kDampeningFactor, double c1, double c2)
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

    // Check the k dampening factor.
    if(kDampeningFactor < 0.8 || kDampeningFactor > 1.2)
    {
        cout << "\n***************************************************************\n";
        cout << "*******       K Dampening Factor Is Out-Of-Range       ********\n";
        cout << "******* K Dampening Factor Must Be In Range [0.8, 1.2] ********\n";
        cout << "*******      Your K Dampening Factor Is Set To: " << kDampeningFactor << "     ********\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Check the c1 and c2 scaling factors.
    if(c1 < 0 || c1 > 2 || c2 < 0 || c2 > 2)
    {
        cout << "\n***************************************************************\n";
        cout << "****       c1 and c2 Scaling Factor Is Out-Of-Range       *****\n";
        cout << "**** c1 and c2 Scaling Factor Must Be In Range [0.0, 2.0] *****\n";
        cout << "****        Your c1 Scaling Factor Is Set To: " << c1 << "         *****\n";
        cout << "****        Your c2 Scaling Factor Is Set To: " << c2 << "         *****\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // Set the Particle Swarm Configuration Parameters.
    psConfig.dimensions = dimensions;
    psConfig.popSize = populationSize;
    psConfig.iterations = maxIterations;
    psConfig.k = kDampeningFactor;
    psConfig.c1 = c1;
    psConfig.c2 = c2;
}

// ------------------------- METHODS ------------------------
/**
 * @brief Runs the Particle Swarm with set parameters.
 *
 * @param functionID The ID that references which Benchmark Function to use.
 * @param minBound, maxBound The minimum and maximum bounds of the individuals in the particle swarm.
 *
 * @return The best global fitness of the Particle Swarm.
 */
double ParticleSwarm::runParticleSwarm(int functionID, double minBound, double maxBound)
{
    // Create a Mersenne Twister pseudo-random number generator.
    mt19937 randGenerator(time(NULL));

    // Record the start time of Particle Swarm.
    auto startTime = chrono::high_resolution_clock::now();

    // Create population with given parameters and evaluate and sort it.
    PS_Population population;
    population.functionID = functionID;
    population.bounds.push_back(minBound);
    population.bounds.push_back(maxBound);
    generatePSPopulation(population, randGenerator);

    // Iterate through the Particle Swarm.
    for(int iteration = 0; iteration < psConfig.iterations; iteration++)
    {
        iteratePopulation(population, randGenerator);

        // Save the global best fitness.
        population.gBestFitnessList.push_back(population.gBestFitness);
    }

    // Record the end time of Particle Swarm.
    auto endTime = chrono::high_resolution_clock::now();

    // Calculate elapsed time in milliseconds it took to run the Particle Swarm.
    auto elapsedTime = endTime - startTime;
    double elapsedTimeMS = chrono::duration_cast<chrono::milliseconds>(elapsedTime).count();

    // Save elapsed time to the population.
    population.executionTime = elapsedTimeMS;

    // Add the population to the population list.
    popList.push_back(population);

    // Return the best fitness of the Particle Swarm.
    int bestFitListSize = population.gBestFitnessList.size();
    return population.gBestFitnessList[bestFitListSize-1];
}

/**
 * @brief Analyzes the results of the Particle Swarm.
 */
void ParticleSwarm::analyzePSResults()
{
    if(popList.empty())
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze  ***********\n";
        cout << "**********   Please Run The Particle Swarm Optimization First    ***********";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new analysis object.
    psAnalysis = PS_Analysis();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].gBestFitnessList.size();

        // Save the function ID.
        psAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].gBestFitnessList);
        psAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].gBestFitnessList);
        psAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].gBestFitnessList[fitnessSize-1]);
        range.push_back(popList[i].gBestFitnessList[0]);
        psAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        psAnalysis.medianFunctionFitness.push_back(popList[i].gBestFitnessList[fitnessSize / 2]);

        // Save the execution time of data.
        psAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        psAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Prints the Results of the Particle Swarm.
 */
void ParticleSwarm::printPSResults()
{
    cout << "****************************************************************************\n";
    cout << "******** Printing Results of Particle Swarm on Current Population **********\n";
    cout << "----------------------------------------------------------------------------\n";

    // If the popList is empty, then Particle Swarm has not been run yet.
    if(popList.size() == 0)
    {
        cout << "********** NO RESULTS FOR THIS POPULATION\n";
        cout << "********** PLEASE RUN THE PARTICLE SWARM\n";
        cout << "----------------------------------------------------------------------------\n\n";
        return;
    }

    for(int i = 0; i < popList.size(); i++)
    {
        cout << "Function ID: " << popList[i].functionID << endl;
        cout << "Generation\t\tBest Fitness of Generation\n";
        cout.precision(12);
        for (int gen = 0; gen < popList[i].gBestFitnessList.size(); gen++)
            cout << gen << "\t\t\t\t" << popList[i].gBestFitnessList[gen] << endl;

        cout << "\nElapsed Time (ms) : " << popList[i].executionTime << endl;
        cout << "----------------------------------------------------------------------------\n\n";
    }
}

/**
 * @brief Prints the Analysis of the Particle Swarm Results.
 */
void ParticleSwarm::printPSAnalysis()
{
    if(psAnalysis.functionIDs.empty())
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
    for(int row = 0; row < psAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << psAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(psAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << psAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(psAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << psAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(psAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << psAnalysis.ranges[row][0] << "\t\t\t";
        if(psAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << psAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(psAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << psAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << psAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << psAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Saves all Particle Swarm Results to file.
 */
void ParticleSwarm::savePSResults()
{
    // If the popList is empty, exit the function.
    if(popList.size() == 0)
    {
        cout << "\n******************************************************\n";
        cout << "****** THERE IS NO PARTICLE SWARM DATA TO SAVE *******";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "ParticleSwarm-Results.csv";

    // Initialize the number of rows (iterations/generations).
    int rows = psConfig.iterations;

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
            line += to_string(popList[pIndex].gBestFitnessList[row]);

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
    if(psAnalysis.functionIDs.size() > 0)
    {
        line = "Average,";
        for(int i = 0; i < psAnalysis.avgFunctionFitness.size(); i++)
        {
            line += to_string(psAnalysis.avgFunctionFitness[i]);

            if(i == psAnalysis.avgFunctionFitness.size()-1)
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
 * @brief Saves the Analysis of the Particle Swarm to file.
 */
void ParticleSwarm::savePSAnalysis()
{
    if(psAnalysis.functionIDs.empty())
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Save *********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = psAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create filename.
    string filename = "ParticleSwarm-Analysis.csv";


    // Create the file where to save the analysis.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    outputFile << psAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(psAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(psAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(psAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(psAnalysis.ranges[row][0]) + ",";
        line += to_string(psAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(psAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(psAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(psAnalysis.functionCalls[row]) + "\n";

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
void ParticleSwarm::saveEndingPopulation()
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
    string filename = "PSO-EndingPop-Function";
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
 * @brief Generates the initial population for Particle Swarm.
 *
 * @note Makes function call to SwarmUtilities.h --> createMatrixMT().
 *
 * @param population The PS_Population structure that holds the particle swarm population.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void ParticleSwarm::generatePSPopulation(PS_Population &population, mt19937 &randGenerator)
{
    double minBound = population.bounds[0];
    double maxBound = population.bounds[1];

    // -------- Create Initial Population. --------
    population.pop = createMatrixMT(psConfig.popSize, psConfig.dimensions, minBound, maxBound, randGenerator);

    // -------- Calculate initial particle velocity. --------
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 0.5*(maxBound - minBound));

    // Set the velocity size to be the same as population size.
    population.velocity.resize(psConfig.popSize);

    // Fill velocity vectors.
    for(int pSize = 0; pSize < psConfig.popSize; pSize++)
    {
        vector<double> pVelocity;
        for(int dSize = 0; dSize < psConfig.dimensions; dSize++)
            pVelocity.push_back(dis(randGenerator));

        // Add particle velocity to velocity matrix.
        population.velocity[pSize] = pVelocity;
    }

    // -------- Calculate the fitness of all particles. --------
    population.fitness.resize(psConfig.popSize);
    evaluatePopulation(population.functionID, population.pop, population.fitness, population.functionCounter);

    // -------- Set personal best for each particle. --------
    population.pBestInd = population.pop;
    population.pBestFitness = population.fitness;

    // -------- Set global best from all the particles. --------
    int gBestIndex = 0;
    for(int pSize = 0; pSize < psConfig.popSize; pSize++)
    {
        if(population.fitness[pSize] < population.fitness[gBestIndex])
            gBestIndex = pSize;
    }
    population.gBestIndividual = population.pop[gBestIndex];
    population.gBestFitness = population.fitness[gBestIndex];
}

/**
 * @brief Calculates the velocity of a particle in the Particle Swarm.
 *
 * @note If c1 > c2 then velocity is changed towards the personal best.
 * @note If c1 < c2 then velocity is changed towards the global best.
 *
 * @param pVelocity The velocity vector of the particle.
 * @param particle The particle (individual) of the population.
 * @param pBest The personal best solution of the particle.
 * @param gBest The global best particle in the population.
 * @param k The velocity dampening factor.
 * @param c1 The scaling factor to bring velocity closer to personal best.
 * @param c2 The scaling factor to bring velocity closer to global best.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void ParticleSwarm::calculateParticleVelocity(vector<double> &pVelocity, vector<double> particle, vector<double> pBest, vector<double> gBest, const double &k, const double &c1, const double &c2, mt19937 &randGenerator)
{
    // Create a distribution for the Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);

    double randNum1, randNum2;

    // Update each element of the particle velocity vector.
    for(int i = 0; i < pVelocity.size(); i++)
    {
        randNum1 = dis(randGenerator);
        randNum2 = dis(randGenerator);
        pVelocity[i] = k * (pVelocity[i] + c1 * randNum1 * (pBest[i] - particle[i]) + c2 * randNum2 * (gBest[i] - particle[i]));
    }
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
void ParticleSwarm::evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter)
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
void ParticleSwarm::evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter)
{
    fitness = calculateFitnessOfVector(indiv, functionID);
    functionCounter++;
}

/**
 * @brief Updates a single particle in the population.
 *
 * @param population The PS_Population structure that holds the particle swarm population.
 * @param particleIndex The index of the particle in the population to update.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void ParticleSwarm::updateParticle(PS_Population &population, int particleIndex, mt19937 &randGenerator)
{
    // Update the velocity of the particle.
    calculateParticleVelocity(population.velocity[particleIndex], population.pop[particleIndex], population.pBestInd[particleIndex], population.gBestIndividual, psConfig.k, psConfig.c1, psConfig.c2, randGenerator);

    // Update the particle.
    for(int i = 0; i < population.pop[particleIndex].size(); i++)
        population.pop[particleIndex][i] = population.pop[particleIndex][i] + population.velocity[particleIndex][i];

    // Calculate the fitness of the particle.
    evaluateIndividual(population.functionID, population.pop[particleIndex], population.fitness[particleIndex], population.functionCounter);

    // Update the particle's personal best if its fitness has improved.
    if(population.fitness[particleIndex] < population.pBestFitness[particleIndex])
    {
        population.pBestInd[particleIndex] = population.pop[particleIndex];
        population.pBestFitness[particleIndex] = population.fitness[particleIndex];
    }
}

/**
 * @brief Iterates the Particle Swarm population to the next generation.
 *
 * @param population The PS_Population structure that holds the particle swarm population.
 * @param randGenerator The Mersenne Twister pseudo-random number generator.
 */
void ParticleSwarm::iteratePopulation(PS_Population &population, mt19937 &randGenerator)
{
    // Update each individual particle in the population.
    for(int p = 0; p < population.pop.size(); p++)
    {
        // Update particle at index p.
        updateParticle(population, p, randGenerator);

        // Check if global best has improved.
        if(population.pBestFitness[p] < population.gBestFitness)
        {
            // Update the global best.
            population.gBestIndividual = population.pBestInd[p];
            population.gBestFitness = population.pBestFitness[p];
        }
    }
}