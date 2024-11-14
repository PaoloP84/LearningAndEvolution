#include "myevoalgo.h"
#include <cmath>

using namespace std;

struct FitnessAndId
{
    /**
     * \brief Fitness
     */
    double fitness;

    /**
     * \brief Id
     */
    int id;
};

/**
 * \brief Lesser-than operator overloading for FitnessAndId
 */
bool operator<(FitnessAndId first, FitnessAndId second)
{
	return (first.fitness < second.fitness);
}

//------------------------------------- MYEVOALGO CLASS ------------------------------------------------------
//---------------------------------- Inherited from EVOGA ----------------------------------------------------

MyEvoAlgo::MyEvoAlgo()
	: farsa::Evoga()
	, m_nevalsteps(1000)
	, m_input()
	, m_inputSize(1)
	, m_numBooleanOperators(1)
	, m_numLayers(1)
	, m_pop()
	, m_boolOps(1)
	, m_output()
	, m_outputSize(1)
	, m_nodes()
	, m_mutationRate(0.2)
	, m_test(false)
	, m_bestGenBuf()
	, m_statBuf()
	, m_logicOpBuf()
	, m_fullLogicOpBuf()
	, m_evalstatBuf()
	, m_saveBestGenInfo(false)
	, m_bestGenInterval(1)
	, m_fitNoise(0.0)
	, m_annealing(false)
	, m_numAnnealingGenerations(0)
	, m_randomOrder(false)
	, m_logicOpArray()
	, m_logicOpInfo(2)
	, m_annCnt(0)
	, m_rng(1)
{
}

MyEvoAlgo::~MyEvoAlgo()
{
}

void MyEvoAlgo::configure(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::Evoga::configure(params, prefix);
	m_nevalsteps = farsa::ConfigurationHelper::getInt(params, prefix + "nevalsteps", m_nevalsteps);
	m_inputSize = farsa::ConfigurationHelper::getInt(params, prefix + "inputSize", m_inputSize);
	if (m_inputSize <= 0)
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "inputSize", QString::number(m_inputSize), "The input size must be >= 1");
	}
	m_numBooleanOperators = farsa::ConfigurationHelper::getInt(params, prefix + "numBooleanOperators", m_numBooleanOperators);
	if (m_numBooleanOperators <= 0)
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "numBooleanOperators", QString::number(m_numBooleanOperators), "The number of boolean operators must be >=1 ");
	}
	m_numLayers = farsa::ConfigurationHelper::getInt(params, prefix + "numLayers", m_numLayers);
	if (m_numLayers <= 0)
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "numLayers", QString::number(m_numLayers), "The number of layers must be >=1 ");
	}
	m_boolOps = farsa::ConfigurationHelper::getInt(params, prefix + "boolOps", m_boolOps);
	if (m_boolOps > 5)
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "boolOps", QString::number(m_boolOps), "The number of logic operators cannot exceed the number of functions! ");
	}
	m_outputSize = farsa::ConfigurationHelper::getInt(params, prefix + "outputSize", m_outputSize);
	if (m_outputSize <= 0)
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "outputSize", QString::number(m_outputSize), "The output size must be >= 1");
	}
	m_mutationRate = farsa::ConfigurationHelper::getDouble(params, prefix + "mutationRate", m_mutationRate);
	if (m_mutationRate > 1.0)
	{
		m_mutationRate /= 100.0;
	}
	m_saveBestGenInfo = farsa::ConfigurationHelper::getBool(params, prefix + "saveBestGenInfo", m_saveBestGenInfo);
	if (m_saveBestGenInfo)
	{
		m_bestGenInterval = farsa::ConfigurationHelper::getInt(params, prefix + "bestGenInterval", m_bestGenInterval);
	}
	m_fitNoise = farsa::ConfigurationHelper::getDouble(params, prefix + "fitNoise", m_fitNoise);
	m_annealing = farsa::ConfigurationHelper::getBool(params, prefix + "annealing", m_annealing);
	if (m_annealing)
	{
		m_numAnnealingGenerations = farsa::ConfigurationHelper::getInt(params, prefix + "numAnnealingGenerations", m_numAnnealingGenerations);
	}
	m_randomOrder = farsa::ConfigurationHelper::getBool(params, prefix + "randomOrder", m_randomOrder);
}

void MyEvoAlgo::describe(QString type)
{
	// Calling parent function
	farsa::Evoga::describe(type);

	Descriptor d = addTypeDescription(type, "A modified version of the steadyState algorithm that uses floating-point genes");
	d.describeInt("inputSize").limits(1,INT_MAX).def(1).help("The number of bits of the input bit-string");
	d.describeInt("numBooleanOperators").limits(1,INT_MAX).def(1).help("The number of boolean operators");
	d.describeInt("numLayers").limits(1,INT_MAX).def(1).help("The number of layers");
	d.describeInt("boolOps").limits(1,INT_MAX).def(1).help("The number of logic operators used to solve the task");
	d.describeInt("outputSize").limits(1,INT_MAX).def(1).help("The number of bits of the output bit-string");
	d.describeReal("mutationRate").def(0.2).help("The mutation rate");
	d.describeInt("resetBufNumGen").limits(1,INT_MAX).def(1).help("The number of generations after which the buffers should be reset");
	d.describeBool("saveBestGenInfo").def(false).help("Flags whether or not the information about the best individual must be saved");
	d.describeInt("bestGenInterval").limits(1,INT_MAX).def(1).help("Interval for saving information about the current best individual");
	d.describeReal("fitNoise").def(0.0).help("The noise applied to the individual's fitness");
	d.describeBool("annealing").def(false).help("Flags whether or not the annealing must be run");
	d.describeInt("numAnnealingGenerations").limits(0,INT_MAX).def(0).help("Number of annealing generations");
	d.describeBool("randomOrder").def(false).help("Flags whether or not randomisation must be performed in order to avoid the replacement of the same worst individuals");
}

void MyEvoAlgo::postConfigureInitialization()
{
	// Calling parent function
	farsa::Evoga::postConfigureInitialization();
	//! Resize the input bit-string
	m_input.resize(m_inputSize);
	//! Resize the output bit-string
	m_output.resize(m_outputSize);
	//! Resize the circuit node array
	m_nodes.resize(m_input.size() + m_numBooleanOperators);
	int numPatterns = (int)pow(2.0,(double)m_inputSize);
	// Parent's stats and selected individuals' stats
	m_logicOpInfo = 2; // Dummy, it has already been done in the initialisation
	if (m_annealing)
	{
		m_logicOpInfo = 3; // Annealing population's stats
	}
	m_logicOpArray.resize(m_logicOpInfo);

    //allocating memory for the mutation vector
    glen = (m_numBooleanOperators * 3) + m_outputSize;

	// Resize the population
	m_pop.resize(popSize * 2); // To take into account the offspring

	for (int i = 0; i < m_pop.size(); i++)
	{
		m_pop[i].resize(glen);
		for (int j = 0; j < glen; j++)
		{
			m_pop[i][j] = j; // Initialisation
		}
	}

    farsa::Logger::info("MyEvoAlgo Configured - Number of genes: " + QString::number(glen));
}

void MyEvoAlgo::evolveAllReplicas()
{
	m_bestGenBuf.clear();
	m_statBuf.clear();
	m_logicOpBuf.clear();
	m_fullLogicOpBuf.clear();
	m_evalstatBuf.clear();

	// It doesn't matter the option chosen in evolutionType,
	// this class will always call my evolutionary Algorithm.
	// If the <test> flag is set to true, it performs a test
	// of the individual and all its one-mutation neighbours
	if (!m_test)
		evolve();
	else
		hamming();
}

farsa::RandomGenerator* MyEvoAlgo::getRng()
{
	return &m_rng;
}

void MyEvoAlgo::setSeed(const int s)
{
	m_rng.setSeed(s);
	currentSeed = s;
}

void MyEvoAlgo::performMutation(int from, int to, int mut)
{
	// If <mut> is 0, no mutations are performed (i.e., simple copy)
    if (mut == 0)
	{
		for (int i = 0; i < glen; i++) 
		{
			m_pop[to][i] = m_pop[from][i];
		}
    }
    else
	{
		// Mutation of the internal logic
		int layer = 0;
        int k = 0;
        while (k < m_numBooleanOperators)
		{
            for (int i = 0; i < (m_numBooleanOperators / m_numLayers); i++)
			{
                for (int j = 0; j < 3; j++)
				{
                    if (j == 0)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry1 = getRng()->getInt(0,(k + m_input.size()) - 1);
							m_pop[to][3 * (i + k)+j] = entry1;
                        }
                        else 
						{
							m_pop[to][3*(i+k)+j] = m_pop[from][3*(i+k)+j];
						}
                    }
                    if (j == 1)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry2 = getRng()->getInt(0,(k + m_input.size())-1);
							m_pop[to][3*(i+k)+j] = entry2;
                        }
                        else 
						{
							m_pop[to][3*(i+k)+j] = m_pop[from][3*(i+k)+j];
						}
                    }
                    if (j == 2)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry3 = getRng()->getInt(0,m_boolOps - 1);
							m_pop[to][3*(i+k)+j] = entry3;
                        }
                        else 
						{
							m_pop[to][3*(i+k)+j] = m_pop[from][3*(i+k)+j];
						}
                    }
                }
            }
			k = k + (m_numBooleanOperators / m_numLayers);
			layer++;
        }
		// Mutation of the output genes
        for (int i = m_numBooleanOperators * 3; i < glen; i++)
		{
			if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
			{
				m_pop[to][i] = getRng()->getInt(m_input.size(), (m_numBooleanOperators + m_input.size() - 1));
            }
			else
			{
				m_pop[to][i] = m_pop[from][i];
			}
        }
    }
}

void MyEvoAlgo::immediateMutation(const QVector<int> from, QVector<int>& to, int mut)
{
	// If <mut> is 0, no mutations are performed (i.e., simple copy)
    if (mut == 0)
	{
		for (int i = 0; i < glen; i++) 
		{
			to[i] = from[i];
		}
    }
    else
	{
		// Mutation of the internal logic
		int layer = 0;
        int k = 0;
        while (k < m_numBooleanOperators)
		{
            for (int i = 0; i < (m_numBooleanOperators / m_numLayers); i++)
			{
                for (int j = 0; j < 3; j++)
				{
                    if (j == 0)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry1 = getRng()->getInt(0,(k + m_input.size()) - 1);
							to[3 * (i + k)+j] = entry1;
                        }
                        else 
						{
							to[3*(i+k)+j] = from[3*(i+k)+j];
						}
                    }
                    if (j == 1)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry1 = getRng()->getInt(0,(k + m_input.size()) - 1);
							to[3 * (i + k)+j] = entry1;
                        }
                        else 
						{
							to[3*(i+k)+j] = from[3*(i+k)+j];
						}
                    }
                    if (j == 2)
					{
                        double number = getRng()->getDouble(0.0,1.0);
						if (number < m_mutationRate)
						{
							int entry3 = getRng()->getInt(0,m_boolOps - 1);
							to[3*(i+k)+j] = entry3;
                        }
                        else 
						{
							to[3*(i+k)+j] = from[3*(i+k)+j];
						}
                    }
                }
            }
			k = k + (m_numBooleanOperators / m_numLayers);
			layer++;
        }
		// Mutation of the output genes
        for (int i = m_numBooleanOperators * 3; i < glen; i++)
		{
			if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
			{
				to[i] = getRng()->getInt(m_input.size(), (m_numBooleanOperators + m_input.size() - 1));
            }
			else
			{
				to[i] = from[i];
			}
        }
    }
}

void MyEvoAlgo::neutralMutation(const QVector<int> from, QVector<int>& to)
{
	// Extract a gene to be mutated (not an output gene)
	bool found = false;
	int geneIdx = -1;
	while (!found)
	{
		geneIdx = getRng()->getInt(0, glen);
		// Security check (perhaps it is useless!!!)
		found = (geneIdx < glen);
	}
	// Store the old value of the gene (to avoid the generation of the same individual)
	const int oldVal = from[geneIdx];
	// Modify the gene
	int currVal;
	bool modified = false;
	if (geneIdx < (glen - m_outputSize))
	{
		// Internal gene
		while (!modified)
		{
			if ((geneIdx % 3) != 2)
			{
				// Operator's input
				int layerCounter = 1;
				int lastValidEntry = m_input.size() - 1;
				bool layerFound = false;
				while (!layerFound)
				{
					if (geneIdx >= layerCounter * (m_numBooleanOperators / m_numLayers) * 3)
					{
						layerCounter++;
						lastValidEntry += (m_numBooleanOperators / m_numLayers);
					}
					else
					{
						layerFound = true;
					}
				}
				currVal = getRng()->getInt(0, lastValidEntry);
			}
			else
			{
				// Operator's type
				currVal = getRng()->getInt(0, (m_boolOps - 1));
			}
			modified = (currVal != oldVal);
		}
	}
	else
	{
		// Output gene
		while (!modified)
		{
			currVal = getRng()->getInt(m_input.size(), (m_numBooleanOperators + m_input.size() - 1));
			modified = (currVal != oldVal);
		}
	}
	for (int g = 0; g < glen; g++)
	{
		if (g == geneIdx)
		{
			to[g] = currVal;
		}
		else
		{
			to[g] = from[g];
		}
	}
}

double MyEvoAlgo::evaluate(const QVector<int> ind)
{
    double fitness = 0.0;
	QVector<int> currOutput(m_outputSize);
    int k = 0;
	int j;
	// Clear the output
	for (j = 0; j < m_output.size(); j++)
	{
		m_output[j] = 0;
	}
	// First of all, copy the input in the circuit node array
	for (j = 0; j < m_input.size(); j++)
	{
		m_nodes[j] = m_input[j];
	}
    while(k < m_numBooleanOperators)
	{
        for (int i = 0; i < (m_numBooleanOperators / m_numLayers); i++)
		{
            int output;
			// Retrieve the operands and the boolean operator
            int op1 = ind[3 * (i + k)];
            int op2 = ind[3 * (i + k) + 1];
            int optype = ind[3 * (i + k) + 2];

            // Compute the output of the boolean operator, given the operands
            switch(optype)
			{
                case 0:
                    output = AND(m_nodes[op1],m_nodes[op2]);
                    break;
                case 1:
                    output = OR(m_nodes[op1],m_nodes[op2]);
                    break;
				case 2:
                    output = NAND(m_nodes[op1],m_nodes[op2]);
                    break;
				case 3:
                    output = NOR(m_nodes[op1],m_nodes[op2]);
                    break;
				case 4:
					output = XOR(m_nodes[op1],m_nodes[op2]);
                    break;
				default:
					farsa::Logger::error("Individual not evaluated! Reason: invalid value for optype " + QString::number(optype) + "!!!");
					exit(-1);
					break;
            }
			// Copy the output value in the circuit node array
			m_nodes[j] = output;
			j++;
        }
        k += (m_numBooleanOperators / m_numLayers);
    }
	// Now compute the output of the evolved circuit
	for (int i = 0; i < currOutput.size(); i++)
	{
		currOutput[i] = m_nodes[ind[ind.size() - m_outputSize + i]];
		m_output[i] = currOutput[i];
	}
	// Expected output
	int expOut = 0;
	int count = 0;
	for (int i = 0; i < m_input.size(); i++)
	{
		if (m_input[i] == 1)
		{
			count++;
		}
	}
	expOut = ((count % 2) == 0) ? 1 : 0;
	// Compute the fitness
	fitness = (currOutput[0] == expOut) ? 1.0 : 0.0;
	// Clear values 
	for (int i = 0; i < m_nodes.size(); i++)
	{
		m_nodes[i] = 0;
	}
    return fitness;
}

void MyEvoAlgo::loadGen(FILE *fp, int ind)
{
	int v;

	fscanf(fp, "DYNAMICAL NN\n");
	for (int g = 0; g < glen; g++) 
	{
		fscanf(fp, "%d\n", &v);
		m_pop[ind][g] = v;
	}
	fscanf(fp, "END\n");
}

int MyEvoAlgo::loadAllGen(int gen, char* filew)
{
	char filename[512];
	char message[512];
	char flag[512];
	FILE* fp;
	if (gen >= 0)
	{
		sprintf(filename, "G%dS%d.gen", gen, seed);
	}
	else
	{
		sprintf(filename, "%s", filew);
	}

	fp = fopen(filename, "r");
	if (fp != NULL) 
	{
		resetPop();
		bool cond = true;
		while (cond) {
			flag[0] = '\0';
			fscanf(fp, "%s : %s\n", flag, message);
			if (strcmp(flag, "**NET") == 0) 
			{
				QVector<int> ind(glen);
				m_pop.append(ind);
				loadGen(fp, m_pop.size() - 1);
			} 
			else 
			{
				cond = false;
			}
		}
		farsa::Logger::info(QString("Loaded ind: %1").arg(m_pop.size()));
		fclose(fp);
	} 
	else 
	{
		farsa::Logger::error(QString("File %1 could not be opened").arg(filename));
	}

	loadedIndividuals = m_pop.size();

	return m_pop.size();
}

double MyEvoAlgo::annealing(const QVector<int> original, double fitness, double actualFitness, int steps, QVector<int>& novel, bool& solved)
{
	double novelFit;
	double currFit;
	// Compute all the bit strings of length <m_input_size>
    	const QVector<QVector<int> > bitStrings = allBitArrays(m_input.size());
	const int numInputs = (int)pow(2.0,(double)m_input.size());
	QVector<int> orig = original;
	QVector<int> evalInd(glen);
	int annGn = 0;
	double fit = 0.0;
	double storedFit = actualFitness;
        int cstep;
	// Reset temporary individual
	invalidFill(evalInd);
	cstep = steps;
	while (annGn < m_numAnnealingGenerations)
	{
		// We must generate new noise at each iteration to the original fitness
		currFit = storedFit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));//+ getRng()->getDouble(-m_fitNoise, m_fitNoise);
		// Mutate the individual
		immediateMutation(orig, evalInd, 1);
		// Evaluate
		double tmpFit = 0.0;
		double noisyFit;
		for (int i = 0; i < numInputs; i++)
		{
			fit = 0.0;
			// Set the input
			for (int j = 0; j < m_input.size(); j++)
			{
				m_input[j] = bitStrings[i][j];
			}
			// evaluate the parent
			fit = evaluate(evalInd);
			tmpFit += fit;
		}
		// Normalise fitness
		tmpFit /= numInputs;
		cstep += numInputs;
		if ((tmpFit == 1.0) && !solved)
		{
			solved = true;
			// Save the number of evaluation steps
			char evalfname[1024];
			sprintf(evalfname, "evalS%d.txt", currentSeed);
			FILE* evalfp = fopen(evalfname, "w");
			if (evalfp != NULL)
			{
				fprintf(evalfp, "%d", cstep);
				fclose(evalfp);
			}
		}
		// Add noise
		noisyFit = tmpFit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));//+ getRng()->getDouble(-m_fitNoise, m_fitNoise);
		/*if (noisyFit < 0.0)
			noisyFit = 0.0;
		if (noisyFit > 1.0)
			noisyFit = 1.0;*/
		if (noisyFit >= currFit)
		{
			immediateMutation(evalInd, orig, 0);
			storedFit = tmpFit;
			currFit = noisyFit;
		}
		annGn++;
	}
	novel = orig;
	novelFit = storedFit;
	return novelFit;
}

int** MyEvoAlgo::generateBitstrings(int n){

    vector<string> input;
    int num_strings = (int)pow(2.0,(double)n);
    input.push_back("0");
    input.push_back("1");
    int i,j = 0;
    for(i = 2; i < (1<<n);i = i<<1){
        for(j = i-1; j>=0 ; j--)
            input.push_back(input[j]);

        for(j=0; j < i;j++)
            input[j] = "0" + input[j];

        for(j=i; j<2*i; j++)
            input[j] = "1" + input[j];
    }

    int** bit_input = new int*[num_strings];
    for(int i = 0; i< (int)pow(2.0,(double)n); i++)
        bit_input[i] = new int[n];


    for(i = 0; i< num_strings;i++){
        for(int j = 0; j < n; j++){
            string tmp ;
            tmp = input[i][j];
            bit_input[i][j] = atoi(tmp.c_str());
        }
    }
    return bit_input;

}

void MyEvoAlgo::evolve()
{
    farsa::Logger::warning("Evolving SteadyState from MyEvoAlgo class");

    int rp; //replication
    int gn; //generation
	int cstep; // evaluation step
    int id; //individuals
    double fit;
    int startGeneration = 0;
    char statfile[128];
    char genFile[128];
    char filename[64];
	bool solved;
	bool finish;
	double bestfit;
	int bestid;
	double cbestfit;
	int cbestid;
	double tmpfit;

    farsa::Logger::info("EVOLUTION: steady state with custom evaluation function");
    farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Compute all the bit strings of length <m_input_size>
    // Store the number of inputs (given n bits --> 2^n input bit-strings)
	const int numInputs = (int)pow(2.0,(double)m_input.size());
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(m_input.size());
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(m_input.size());
		for (int j = 0; j < m_input.size(); j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}
	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    for(rp = 0; rp < nreplications; rp++) 
	{	
        startGeneration = 0;
		setSeed(getStartingSeed() + rp);
        farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
        resetGenerationCounter();
		// Initialise the population (i.e. only the parents, offspring are mutated copies)
        initPop();
        // Set fbest to a very low value
        this->fbest = -99999.0;

        emit startingReplication( rp );

        QTime evotimer;
        evotimer.start();
		// Reset tfitness and ntfitness
        for (int i = 0; i < m_pop.size(); i++) 
		{
            tfitness[i] = 0.0;
            ntfitness[i] = 0.0;
        }
        //code to recovery a previous evolution: Experimental
        sprintf(statfile, "statS%d.fit", getStartingSeed() + rp);
        //now check if the file exists
        farsa::DataChunk statTest(QString("stattest"), Qt::blue, 2000, false);
        if (statTest.loadRawData(QString(statfile),0))
		{
            startGeneration = statTest.getIndex();
            sprintf(genFile,"G%dS%d.gen", startGeneration, getStartingSeed() + rp);
            farsa::Logger::info("Recovering from startGeneration: " + QString::number(startGeneration));
            farsa::Logger::info(QString("Loading file: ") + genFile);
            loadAllGen(-1, genFile);
            cgen = startGeneration;
            emit recoveredInterruptedEvolution( QString(statfile) );
        } //end evolution recovery code

		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_logicOpBuf = "";
		m_fullLogicOpBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpLogicOpBuf = "";
		QString tmpFullLogicOpBuf = "";
		QString tmpEvalStatBuf = "";

		// generations
		gn = startGeneration;
		cstep = 0;
		solved = false;
		finish = false;
		bestfit = -9999.0;
		bestid = -1;
		while (cstep < m_nevalsteps)// && !solved)
        //for(gn = startGeneration; gn < nogenerations; gn++) 
		{
			evotimer.restart();
			vector<int> identity;
            farsa::Logger::info(" Generation " + QString::number(gn + 1));
            exp->initGeneration(gn);
            if ( commitStep() ) 
			{ 
				return;
			}
			// initialise logic operators' statistics
			for (int i = 0; i < m_logicOpInfo; i++)
			{
				m_logicOpArray[i] = 0.0;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				identity.push_back(id);
				performMutation(id, popSize + id, 1);
			}
			cbestfit = -9999.0;
			cbestid = -1;
			for(id = 0; id < popSize; id++) 
			{
				// Test the individual and its child
				// Parents are tested only at the first generation
				if (gn == startGeneration)
				{
					tfitness[id] = 0.0;
					ntfitness[id] = 0.0;
					for (int i = 0; i < numInputs; i++)
					{
						fit = 0.0;
						// Set the input
						for (int j = 0; j < m_input.size(); j++)
						{
							m_input[j] = bitStrings[i][j];
						}
						// evaluate the parent
						fit = evaluate(m_pop[id]);
						// update its fitness
						tfitness[id] += fit;
						ntfitness[id]++;
						if (isStopped()) 
						{ // stop evolution
							return;
						}
						cstep++;
					}
					tmpfit = tfitness[id] / ntfitness[id];
					if ((tmpfit == 1.0) && !solved)
					{
						solved = true;
						// Save the number of evaluation steps
						char evalfname[1024];
						sprintf(evalfname, "evalS%d.txt", currentSeed);
						FILE* evalfp = fopen(evalfname, "w");
						if (evalfp != NULL)
						{
							fprintf(evalfp, "%d", cstep);
							fclose(evalfp);
						}
					}
					// Check whether or not individual is better than current best
					if (tmpfit > cbestfit)
					{
						cbestfit = tmpfit;
						cbestid = id;
					}
				}
				tfitness[popSize + id] = 0.0;
				ntfitness[popSize + id] = 0.0;
				for (int i = 0; i < numInputs; i++)
				{
					fit = 0.0;
					// Set the input
					for (int j = 0; j < m_input.size(); j++)
					{
						m_input[j] = bitStrings[i][j];
					}
					// evaluate the child
					fit = evaluate(m_pop[popSize + id]);
					// update its fitness
					tfitness[popSize + id] += fit;
					ntfitness[popSize + id]++;
					if (isStopped()) 
					{ // stop evolution
						return;
					}
					cstep++;
				}
				tmpfit = tfitness[popSize + id] / ntfitness[popSize + id];
				if ((tmpfit == 1.0) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				if (tmpfit > cbestfit)
				{
					cbestfit = tmpfit;
					cbestid = popSize + id;
				}
			}
            exp->endGeneration(gn);
            if ( commitStep() ) 
			{ 
				return; 
			}
			
			// ========= Selection part ========== //
            // Finally, we look for the worst individuals (parents) and substitute them with best children.
            // What we do is: we order both the parents and the children in descending order, then we take the
            // popSize best individuals. This is not the same as what is done in the sequential version of
            // the algorithm, but should be similar. We overwrite the worst parents with the best children
			if (m_randomOrder)
				random_shuffle(identity.begin(),identity.end());
            QVector<FitnessAndId> parents(popSize);
            QVector<FitnessAndId> children(popSize);
            for (int i = 0; i < popSize; i++)
			{
				parents[i].fitness = tfitness[identity[i]] / ntfitness[identity[i]];
                parents[i].id = identity[i];
                /*parents[i].fitness = tfitness[i] / ntfitness[i];
                parents[i].id = i;*/
				parents[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));//+= getRng()->getDouble(-m_fitNoise, m_fitNoise);
				/*if (parents[i].fitness < 0.0)
					parents[i].fitness = 0.0;
				if (parents[i].fitness > 1.0)
					parents[i].fitness = 1.0;*/
            }

            for (int i = 0; i < popSize; i++)
			{
				children[i].fitness = tfitness[popSize + identity[i]] / ntfitness[popSize + identity[i]];
                children[i].id = popSize + identity[i];
                /*children[i].fitness = tfitness[popSize + i] / ntfitness[popSize + i];
                children[i].id = popSize + i;*/
				children[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));//+= getRng()->getDouble(-m_fitNoise, m_fitNoise);
				/*if (children[i].fitness < 0.0)
					children[i].fitness = 0.0;
				if (children[i].fitness > 1.0)
					children[i].fitness = 1.0;*/
            }
			// Compute the average number of logic operators for what concerns the initial population
			QVector<int> parentLogicOp(popSize);
			QVector<int> parentSingleLogicOp(popSize);
			float parentLogicOpCount = 0.0;
			float parentSingleLogicOpCount = 0.0;
			for (int i = 0; i < popSize; i++)
			{
				QVector<int> tree;
				QVector<int> singleOpTree;
				pathLogicOperator(m_pop[i][(glen - 1)], m_pop[i], tree, 0);
				countLogicOperator(tree, singleOpTree);
				parentLogicOp[i] = tree.size();
				parentLogicOpCount += parentLogicOp[i];
				parentSingleLogicOp[i] = singleOpTree.size();
				parentSingleLogicOpCount += parentSingleLogicOp[i];
			}
			parentLogicOpCount /= popSize;
			parentSingleLogicOpCount /= popSize;
			m_logicOpArray[0] = parentLogicOpCount;
			// Sorting both parents and children. They are sorted in ascending order but we need the best
            // individuals (those with the highest fitness) first
            qSort(parents);   //Error de operadores aqui#include <QTime>
            qSort(children);
			// Store the list of selected and replaced individuals
			QVector<int> selected(popSize);
			QVector<double> selectedFit(popSize);
			QVector<int> replaced(popSize);
            int p = popSize - 1;
            int c = popSize - 1;
            for (int i = 0; i < popSize; i++)
			{
                if (parents[p].fitness > children[c].fitness) 
				{
                    selected[i] = parents[p].id;
					selectedFit[i] = parents[p].fitness;
					replaced[i] = parents[p].id;
					p--;
            	} 
				else 
				{
					selected[i] = children[c].id;
					selectedFit[i] = children[c].fitness;
					replaced[i] = parents[popSize - 1 - c].id;
                    c--;
                }
            }
			// Check whether or not to perform annealing
			if (m_annealing)
			{
				// Annealing
				// Parents
				for (int i = 0; i < popSize; i++)
				{
					// Get the individual index
					const int index = selected[i];
					m_annCnt = index;
					QVector<int> original(glen);
					immediateMutation(m_pop[index], original, 0);
					double fitness = selectedFit[i];
					double actualFitness = tfitness[index] / ntfitness[index];
					QVector<int> novel(glen);
					double novelFit = annealing(original, fitness, actualFitness, cstep, novel, solved);
					// Update the individual
					immediateMutation(novel, m_pop[index], 0);
					// Update structure data
					selectedFit[i] = novelFit;
					tfitness[selected[i]] = selectedFit[i];
					ntfitness[selected[i]] = 1.0;
					// Check whether or not the individual returned by the annealing routine is better than current best
					if (novelFit > cbestfit)
					{
						cbestfit = novelFit;
						cbestid = index;
					}
					cstep += (numInputs * m_numAnnealingGenerations);
				}
			}
			if (cbestfit > bestfit)
			{
				bestfit = cbestfit;
				bestid = cbestid;
			}
			// Swap individuals
            for (int i = 0; i < popSize; i++)
			{
                performMutation(selected[i], replaced[i], 0);
                tfitness[replaced[i]] = tfitness[selected[i]];
                ntfitness[replaced[i]] = ntfitness[selected[i]];
            }
			// Compute and save statistics
            computeFStat2();
            saveFStat();
			// Count the average number of logic operators (i.e., how complex the circuit is)
			QVector<int> logicOp(popSize);
			QVector<int> singleLogicOp(popSize);
			float logicOpCount = 0.0;
			float singleLogicOpCount = 0.0;
			for (int i = 0; i < popSize; i++)
			{
				QVector<int> tree;
				QVector<int> singleOpTree;
				pathLogicOperator(m_pop[i][(glen - 1)], m_pop[i], tree, 0);
				countLogicOperator(tree, singleOpTree);
				logicOp[i] = tree.size();
				logicOpCount += logicOp[i];
				singleLogicOp[i] = singleOpTree.size();
				singleLogicOpCount += singleLogicOp[i];
			}
			logicOpCount /= popSize;
			singleLogicOpCount /= popSize;
			m_logicOpArray[1] = logicOpCount;
			saveLogicOperatorCounter(logicOpCount, singleLogicOpCount);
			saveFullLogicOpInfo();
            emit endGeneration( cgen, fmax, faverage, fmin );
            if (commitStep()) 
			{
                return; // stop the evolution process
            }

            cgen++;

			// Save eval stats
			saveEvalStats(cstep, bestfit);

			farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmax).arg(bestfit).arg(cstep));
            fflush(stdout);

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
				finish = true;

			if (finish)//if (solved || finish)
			{
				// Save the number of evaluation steps
				if (!solved)
				{
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}

				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				tmpLogicOpBuf += m_logicOpBuf;
				{
				   const QString logicOpFilename = QString("logicOperatorS%1.txt").arg(currentSeed);
				   QFile logicOpFile(logicOpFilename);
				   if (!logicOpFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
				       farsa::Logger::warning("Cannot save stats into " + logicOpFilename);
				   } else {
					   QTextStream s(&logicOpFile);
					   s << tmpLogicOpBuf;
				   }
				   m_logicOpBuf.clear();
				}

				tmpFullLogicOpBuf += m_fullLogicOpBuf;
				{
				   const QString fullLogicOpFilename = QString("fullLogicOperatorS%1.txt").arg(currentSeed);
				   QFile fullLogicOpFile(fullLogicOpFilename);
				   if (!fullLogicOpFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
				       farsa::Logger::warning("Cannot save stats into " + fullLogicOpFilename);
				   } else {
					   QTextStream s(&fullLogicOpFile);
					   s << tmpFullLogicOpBuf;
				   }
				   m_fullLogicOpBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					m_evalstatBuf.clear();
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        }
        saveAllGen();
    }
}

void MyEvoAlgo::hamming()
{
	int rp; //replication
	int startGeneration;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;
    	double h;
	int nh;
	int ng;
	QVector<int> cind;
	QVector<int> oind;
	int i;
	int j;
	int g;
	farsa::Logger::info("Compute Hamming distance");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		//code to recovery a previous evolution: Experimental
		sprintf(statfile, "statS%d.fit", getStartingSeed() + rp);
		//now check if the file exists
		farsa::DataChunk statTest(QString("stattest"), Qt::blue, 2000, false);
		if (statTest.loadRawData(QString(statfile),0))
		{
		    startGeneration = statTest.getIndex();
		    sprintf(genFile,"G%dS%d.gen", startGeneration, getStartingSeed() + rp);
		    farsa::Logger::info("Recovering from startGeneration: " + QString::number(startGeneration));
		    farsa::Logger::info(QString("Loading file: ") + genFile);
		    loadAllGen(-1, genFile);
		    cgen = startGeneration;
		    emit recoveredInterruptedEvolution( QString(statfile) );
		} //end evolution recovery code

		// Hamming distance between all pairs of individuals in the population
		h = 0.0;
		nh = 0;
		i = 0;
		while (i < (popSize - 1))
		{
			cind = m_pop[i];
			j = (i + 1);
			while (j < popSize)
			{
				oind = m_pop[j];
				ng = 0;
				for (g = 0; g < glen; g++)
				{
					if (cind[g] != oind[g])
						ng++;
				}
				h += ((double)ng / (double)glen);
				nh++;
				j++;
			}
			i++;
		}
		if (nh > 0)
			h /= ((double)nh);

		// Save data about robot preferences/types
		sprintf(filename, "hammingS%d.txt", currentSeed);
		fp = fopen(filename, "w");
		if (fp != NULL)
		{
			fprintf(fp, "%lf\n", h);
			fclose(fp);
		}
	}		
}

QVector<int> MyEvoAlgo::oneBitSum(const QVector<int> bitArray, int n)
{
	QVector<int> retBitArray(n);
	int i = n - 2;
	int remainder = 0;
	int tmpSum = (bitArray[n - 1] + 1) % 2;
	// Check whether or not there is a remainder
	// (remember that the sum of two bits could
	// be either 0 or 1)
	if (tmpSum == 0)
	{
        remainder = 1;
	}
	retBitArray[n - 1] = tmpSum;
	// Compute the remaining bits of the bit-array
	// by taking into account the remainder
	while (i >= 0)
    {
        retBitArray[i] = bitArray[i];
        if (remainder > 0)
        {
            retBitArray[i] = (retBitArray[i] + 1) % 2;
            if (retBitArray[i] == 0)
            {
                remainder = 1;
            }
            else
            {
                remainder = 0;
            }
        }
		i--;
    }
	return retBitArray;
}

QVector<QVector<int> > MyEvoAlgo::allBitArrays(int n)
{
	// Compute the number of combinations
	// (there are 2^n combinations of bit-array
	// of length n)
	const int numCombs = (int)pow(2.0,(double)n);
	int j;
	// The bit matrix
	QVector<QVector<int> > bitArray(numCombs);
	for (int i = 0; i < numCombs; i++)
	{
		bitArray[i].resize(n);
	}
	// Compute the combinations
	// The first array contains all zeroed values
	for (int i = 0; i < n; i++)
	{
		bitArray[0][i] = 0;
	}
	// All the remaining arrays are obtained by summing
	// up the previous bit sequence and 1
	j = 1;
	while (j < numCombs)
	{
		bitArray[j] = oneBitSum(bitArray[j-1], n);
		j++;
	}
	return bitArray;
}

void MyEvoAlgo::initPop()
{
	// Valid output ranges
	const int minOutRange = m_input.size();
	int maxOutRange = m_input.size() + m_numBooleanOperators - 1;
	for (int i = 0; i < m_pop.size(); i++)
	{
		int j = 0;
		int layer = 0;
		while (j < m_numBooleanOperators)
		{
			for (int k = 0; k < (m_numBooleanOperators / m_numLayers); k++)
			{
				for (int l = 0; l < 3; l++)
				{
					if (l == 0)
					{
						m_pop[i][3 * (j + k) + l] = getRng()->getInt(0, (j + minOutRange - 1));
					}
					if (l == 1)
					{
						m_pop[i][3 * (j + k) + l] = getRng()->getInt(0, (j + minOutRange - 1));
					}
					if (l == 2)
					{
						m_pop[i][3 * (j + k) + l] = getRng()->getInt(0, m_boolOps - 1);
					}
				}
			}
			j = j + (m_numBooleanOperators / m_numLayers);
			layer++;
		}
		for (j = m_numBooleanOperators * 3; j < glen; j++)
		{
			m_pop[i][j] = getRng()->getInt(minOutRange, maxOutRange);
		}
	}
}

void MyEvoAlgo::resetPop()
{
	for (int i = 0; i < m_pop.size(); i++)
	{
		m_pop[i].resize(0);
	}
	m_pop.resize(0);
}

QVector<int> MyEvoAlgo::circularShift(QVector<int> bitArray)
{
	const int len = bitArray.size();
	QVector<int> retBitArray(len);
	int i = 1;
	retBitArray[0] = bitArray[len - 1];
	while (i < len)
	{
		retBitArray[i] = bitArray[i - 1];
		i++;
	}
	return retBitArray;
}

QVector<int> MyEvoAlgo::bitSum(QVector<int> bitArray)
{
	const int len = bitArray.size() / 2;
	QVector<int> retBitArray(len + 1);
	int i = len - 1;
	int remainder = 0;
	int tmpSum;
	while (i >= 0)
	{
		tmpSum = (bitArray[i] + bitArray[len + i] + remainder) % 2;
		if (tmpSum == 0)
		{
			if ((bitArray[i] != 0) || (bitArray[len + i] != 0))
			{
				remainder = 1;
			}
		}
		else
		{
			if ((bitArray[i] == 1) && (bitArray[len + i] == 1))
			{
				remainder = 1;
			}
			else
			{
				remainder = 0;
			}
		}
		retBitArray[i] = tmpSum;
		i--;
	}
	retBitArray[0] = remainder;
	return retBitArray;
}

QVector<int> MyEvoAlgo::bitMultiplier(QVector<int> bitArray)
{
	const int len = bitArray.size();
	QVector<int> retBitArray(len);
	const int halfLen = len / 2;
	int m1 = 0;
	int m2 = 0;
	int i = 0;
	while (i < halfLen)
	{
		m1 += bitArray[halfLen - 1 - i] * ((int)pow(2, (double)i));
		m2 += bitArray[len - 1 - i] * ((int)pow(2, (double)i));
		i++;
	}
	int p = m1 * m2;
	i = len - 1;
	while (i >= 0)
	{
		retBitArray[i] = (p % 2);
		p /= 2;
		i--;
	}
	return retBitArray;
}

QVector<int> MyEvoAlgo::evenParity(QVector<int> bitArray)
{
	QVector<int> parity(1);
	int count = 0;
	parity[0] = 0;
	for (int i = 0; i < bitArray.size(); i++)
	{
		if (bitArray[i] == 1)
		{
			count++;
		}
	}
	if ((count % 2) == 0)
	{
		parity[0] = 1;
	}
	return parity;
}

QVector<int> MyEvoAlgo::oddParity(QVector<int> bitArray)
{
	QVector<int> parity(1);
	int count = 0;
	parity[0] = 0;
	for (int i = 0; i < bitArray.size(); i++)
	{
		if (bitArray[i] == 1)
		{
			count++;
		}
	}
	if ((count % 2) != 0)
	{
		parity[0] = 1;
	}
	return parity;
}

int MyEvoAlgo::booleanAnd( int a, int b ) 
{
	return (((a == 1) && (b == 1)) ? 1 : 0);
}

int MyEvoAlgo::booleanOr( int a, int b ) 
{
	return (((a == 1) || (b == 1)) ? 1 : 0);
}

int MyEvoAlgo::booleanXor( int a, int b ) 
{
	return ((a != b) ? 1 : 0);
}

int MyEvoAlgo::booleanNand( int a, int b ) 
{
	return (1 - booleanAnd(a,b));
}

int MyEvoAlgo::booleanNor( int a, int b )
{
	return (1 - booleanOr(a,b));
}

void MyEvoAlgo::saveBestInd()
{
	//to do
	//selecting best fathers
	int i;
	int bi;
	double bn;
	double ffit;
	bn = -999999.0;// it was 0 possible source of bug in case of negative fitness
	bi = -1;

	//finding the best simply the best, one individual
	for(i = 0; i < this->popSize; i++) 
	{
		ffit = this->tfitness[i] / this->ntfitness[i];
		if(ffit > bn) 
		{
			bn = ffit;
			bi = i;
		}
	}

	//now saving
	m_bestGenBuf += QString("**NET : s%1_%2.wts\n").arg(cgen).arg(bi);
	m_bestGenBuf += "DYNAMICAL NN\n";
	for (int j = 0; j < glen; j++) {
		m_bestGenBuf += QString("%1\n").arg(m_pop[bi][j]);
	}
	m_bestGenBuf += "END\n";
}

void MyEvoAlgo::saveGen(FILE *fp, int ind)
{
    int j;
    fprintf(fp, "DYNAMICAL NN\n");
    for (j = 0; j < glen; j++)
	{
        fprintf(fp, "%d\n", m_pop[ind][j]);
	}
    fprintf(fp, "END\n");
}

void MyEvoAlgo::saveAllGen()
{
    FILE *fp;
    char filename[64];
    int i;

    sprintf(filename, "G%dS%d.gen", cgen, currentSeed);
    if ((fp = fopen(filename, "w+")) == NULL) 
	{
        farsa::Logger::error(QString("Cannot open file %1").arg(filename));
    } 
	else {
        //we save
        for(i = 0; i < popSize; i++) 
		{
            fprintf(fp, "**NET : %d_%d_%d.wts\n", cgen, 0, i);
            saveGen(fp, i);
        }
        fclose( fp );
    }
}

void MyEvoAlgo::saveFStat()
{
	m_statBuf += QString("%1 %2 %3\n").arg(fmax).arg(faverage).arg(fmin);
}

void MyEvoAlgo::saveLogicOperatorCounter(float logicOpCount, float singleLogicOpCount)
{
	m_logicOpBuf += QString("%1 %2 \n").arg(logicOpCount).arg(singleLogicOpCount);
}

void MyEvoAlgo::invalidFill(QVector<int>& ind)
{
	for (int i = 0; i < ind.size(); i++)
	{
		ind[i] = -1;
	}
}

void MyEvoAlgo::saveFullLogicOpInfo()
{
	if (m_annealing)
	{
		m_fullLogicOpBuf += QString("%1 %2 %3\n").arg(m_logicOpArray[0]).arg(m_logicOpArray[1]).arg(m_logicOpArray[2]);
	}
	else
	{
		m_fullLogicOpBuf += QString("%1 %2\n").arg(m_logicOpArray[0]).arg(m_logicOpArray[1]);
	}
}

void MyEvoAlgo::saveEvalStats(const int steps, const double fit)
{
	m_evalstatBuf += QString("%1 %2\n").arg(steps).arg(fit);
}

/*
 * Code for printing the binary tree of the evolved logic circuit
 */

int MyEvoAlgo::maxHeight(BinaryTree *p) 
{
	if (!p) 
		return 0;
	int leftHeight = maxHeight(p->left);
	int rightHeight = maxHeight(p->right);
	return ((leftHeight > rightHeight) ? (leftHeight + 1) : (rightHeight + 1));
}

string MyEvoAlgo::intToString(int val) 
{
	ostringstream ss;
	ss << val;
	return ss.str();
}

void MyEvoAlgo::printBranches(int branchLen, int nodeSpaceLen, int startLen, int nodesInThisLevel, const deque<BinaryTree*>& nodesQueue, ostream& out) 
{
	deque<BinaryTree*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel / 2; i++) 
	{
		out << ((i == 0) ? setw(startLen - 1) : setw(nodeSpaceLen - 2)) << "" << ((*iter++) ? "/" : " ");
		out << setw(2 * branchLen + 2) << "" << ((*iter++) ? "\\" : " ");
	}
	out << endl;
}

void MyEvoAlgo::printNodes(int branchLen, int nodeSpaceLen, int startLen, int nodesInThisLevel, const deque<BinaryTree*>& nodesQueue, ostream& out) 
{
	deque<BinaryTree*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel; i++, iter++) 
	{
		out << ((i == 0) ? setw(startLen) : setw(nodeSpaceLen)) << "" << ((*iter && (*iter)->left) ? setfill('_') : setfill(' '));
		out << setw(branchLen + 2) << ((*iter) ? intToString((*iter)->data) : "");
		out << ((*iter && (*iter)->right) ? setfill('_') : setfill(' ')) << setw(branchLen) << "" << setfill(' ');
	}
	out << endl;
}

// Print the leaves only (just for the bottom row)
void MyEvoAlgo::printLeaves(int indentSpace, int level, int nodesInThisLevel, const deque<BinaryTree*>& nodesQueue, ostream& out) 
{
	deque<BinaryTree*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel; i++, iter++) 
	{
		out << ((i == 0) ? setw(indentSpace + 2) : setw(2 * level + 2)) << ((*iter) ? intToString((*iter)->data) : "");
	 }
	out << endl;
}

// Pretty formatting of a binary tree to the output stream
// @ param
// level  Control how wide you want the tree to sparse (eg, level 1 has the minimum space between nodes, while level 2 has a larger space between nodes)
// indentSpace  Change this to add some indent space to the left (eg, indentSpace of 0 means the lowest level of the left node will stick to the left margin)
void MyEvoAlgo::printPretty(BinaryTree *root, int level, int indentSpace, ostream& out) 
{
	int h = maxHeight(root);
	int nodesInThisLevel = 1;

	int branchLen = 2 * ((int)pow(2.0, h) - 1) - (3 - level) * (int)pow(2.0, h-1);  // eq of the length of branch for each node of each level
	int nodeSpaceLen = 2 + (level + 1) * (int)pow(2.0, h);  // distance between left neighbor node's right arm and right neighbor node's left arm
	int startLen = branchLen + (3 - level) + indentSpace;  // starting space to the first node to print of each level (for the left most node of each level only)

	deque<BinaryTree*> nodesQueue;
	nodesQueue.push_back(root);
	for (int r = 1; r < h; r++) 
	{
		printBranches(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);
		branchLen = branchLen/2 - 1;
		nodeSpaceLen = nodeSpaceLen/2 + 1;
		startLen = branchLen + (3 - level) + indentSpace;
		printNodes(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);

		for (int i = 0; i < nodesInThisLevel; i++) 
		{
			BinaryTree *currNode = nodesQueue.front();
			nodesQueue.pop_front();
			if (currNode) 
			{
				nodesQueue.push_back(currNode->left);
				nodesQueue.push_back(currNode->right);
			} 
			else 
			{
				nodesQueue.push_back(NULL);
				nodesQueue.push_back(NULL);
			}
		}
		nodesInThisLevel *= 2;
	}
	printBranches(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);
	printLeaves(indentSpace, level, nodesInThisLevel, nodesQueue, out);
}

void MyEvoAlgo::treeLogicOperator(int node, QVector<int> genotype, BinaryTree* root)
{
	const int bit_length = m_input.size();
    if ((genotype[(node - bit_length) * 3] <= (bit_length - 1)) && (genotype[(node - bit_length) * 3 + 1] <= bit_length - 1))
	{
        root->left = new BinaryTree(genotype[(node - bit_length) * 3]);
        root->right = new BinaryTree(genotype[(node - bit_length) * 3 + 1]);
        return;
    }
	else
    {
        if (genotype[(node - bit_length) * 3] > bit_length - 1)
		{
            root->left = new BinaryTree(genotype[(node - bit_length) * 3]);
            treeLogicOperator(genotype[(node - bit_length) * 3], genotype, root->left);
            if (genotype[(node - bit_length) * 3 + 1] <= bit_length - 1)
			{
                root->right = new BinaryTree(genotype[(node - bit_length) * 3 + 1]);
            }
        }
        if (genotype[(node - bit_length) * 3 + 1] > bit_length - 1)
		{
            root->right = new BinaryTree(genotype[(node - bit_length) * 3 + 1]);
            treeLogicOperator(genotype[(node - bit_length) * 3 + 1], genotype, root->right);
            if (genotype[(node - bit_length) * 3] <= bit_length - 1)
			{
                root->left = new BinaryTree(genotype[(node - bit_length) * 3]);
            }
        }
    }
}

void MyEvoAlgo::pathLogicOperator(int node, QVector<int> genotype, QVector<int>& tree, int type)
{
	const int bit_length = m_input.size();
	tree.append(node);
    if(type == 1)
	{
        tree.append(genotype[(node - bit_length) * 3]);
		tree.append(genotype[(node - bit_length) * 3 + 1]);
		tree.append(genotype[(node - bit_length) * 3 + 2]);
    } 
    if (genotype[(node - bit_length) * 3] <= (bit_length - 1) && genotype[(node - bit_length) * 3 + 1] <= (bit_length - 1))
	{
        return;
    }
	else
    {
        if (genotype[(node - bit_length) * 3] > (bit_length - 1))
		{
            pathLogicOperator(genotype[(node - bit_length) * 3], genotype, tree, type);
        }
        if (genotype[(node - bit_length) * 3 + 1] > (bit_length - 1))
		{
            pathLogicOperator(genotype[(node-bit_length)* 3 + 1], genotype, tree, type);
        }
    }
}

void MyEvoAlgo::typeLogicOperator(int node, QVector<int> genotype, QVector<int>& boolean_function)
{
	const int bit_length = m_input.size();
    boolean_function.append(genotype[(node - bit_length) * 3 + 2]);
    if (genotype[(node - bit_length) * 3] <= (bit_length - 1) && genotype[(node - bit_length) * 3 + 1] <= (bit_length - 1))
	{
        return;
    }
	else
    {
        if (genotype[(node - bit_length) * 3] > (bit_length - 1))
		{
			typeLogicOperator(genotype[(node - bit_length) * 3], genotype, boolean_function);
        }
        if (genotype[(node - bit_length) * 3 + 1] > (bit_length - 1))
		{
            typeLogicOperator(genotype[(node - bit_length) * 3 + 1], genotype, boolean_function);
		}
    }
}

void MyEvoAlgo::countLogicOperator(QVector<int> tree, QVector<int>& singleOpTree)
{
	const int bit_length = m_input.size();
    for (int i = 0; i < tree.size(); i++)
	{
        int count = 0;
        int n = singleOpTree.size();
        if(tree[i] < bit_length)
		{
            count++;
		}
		for (int j = 0; j < n; j++)
		{
            if (tree[i] == singleOpTree[j])
			{
                count++;
            }
		}
        if (count == 0)
		{
            singleOpTree.push_back(tree[i]);
        }
	}
}
