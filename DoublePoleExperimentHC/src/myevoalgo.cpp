#include "myevoalgo.h"
#include "evonet.h"
#include <cmath>
#include <QTextStream>
#include <QFile>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

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
	, m_pop()
	, m_mutationRate(0.2)
	, m_bestGenBuf()
	, m_statBuf()
	, m_gstatBuf()
	, m_evalstatBuf()
	, m_saveBestGenInfo(false)
	, m_bestGenInterval(1)
	, m_fitNoise(0.0)
	, m_cartPole()
	, m_taskType(0)
	, m_taskEpisodes(1)
	, m_generalize(false)
	, m_rng(1)
	, m_calcCorr(false)
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
	m_taskType = farsa::ConfigurationHelper::getInt(params, prefix + "taskType", m_taskType);
	if ((m_taskType < 0) || (m_taskType >= NTASKS))
	{
		farsa::ConfigurationHelper::throwUserConfigError(prefix + "taskType", QString::number(m_taskType), "Task type can be one of the following values: 0, 1, 2 or 3!!!");
	}
	m_generalize = farsa::ConfigurationHelper::getBool(params, prefix + "generalize", m_generalize);
	m_calcCorr = farsa::ConfigurationHelper::getBool(params, prefix + "calcCorr", m_calcCorr);
}

void MyEvoAlgo::describe(QString type)
{
	// Calling parent function
	farsa::Evoga::describe(type);

	Descriptor d = addTypeDescription(type, "A modified version of the steadyState algorithm that uses floating-point genes");
    d.describeReal("mutationRate").def(0.2).help("The mutation rate");
	d.describeBool("saveBestGenInfo").def(false).help("Flags whether or not the information about the best individual must be saved");
	d.describeInt("bestGenInterval").limits(1,INT_MAX).def(1).help("Interval for saving information about the current best individual");
	d.describeReal("fitNoise").def(0.0).help("The noise applied to the individual's fitness");
	d.describeInt("taskType").limits(0,NTASKS).def(0).help("Task type");
	d.describeBool("generalize").def(false).help("Flags whether or not generalization test must be run");
}

void MyEvoAlgo::postConfigureInitialization()
{
	// Calling parent function
	farsa::Evoga::postConfigureInitialization();
	// Set episodes
	m_taskEpisodes = EVAL_EPS;
	// Initialize cart pole
	bool fixedState = true;
	bool markov = true;
	if ((m_taskType == 1) || (m_taskType == 3))
		fixedState = false;
	if ((m_taskType == 2) || (m_taskType == 3))
		markov = false;
	m_cartPole.nmarkov_long = false;
	m_cartPole.generalization_test = false;
    	m_cartPole.testCart = false;
    	m_cartPole.initializeCartPole(markov, fixedState);
	// Resize population
	m_pop.resize(popSize * 2);
	for (int i = 0; i < (popSize * 2); i++)
	{
		m_pop[i].resize(glen);
	}
}

void MyEvoAlgo::evolveAllReplicas()
{
	m_bestGenBuf.clear();
	m_statBuf.clear();
	m_gstatBuf.clear();
	m_evalstatBuf.clear();
	
	stopEvolution = false;
	// It doesn't matter the option chosen in evolutionType,
	// this class will always call my evolutionary Algorithm.
	evolve();
}

void MyEvoAlgo::setSeed(const int s)
{
	m_rng.setSeed(s);
	currentSeed = s;
}

farsa::RandomGenerator* MyEvoAlgo::getRng()
{
	return &m_rng;
}

void MyEvoAlgo::initPop()
{
	for (int i = 0; i < m_pop.size(); i++)
	{
		for (int j = 0; j < glen; j++)
		{
			m_pop[i][j] = getRng()->getInt(0,255);
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

void MyEvoAlgo::performMutation(int from, int to, int mut)
{
	bool ok;
	int num;
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
		// Mutation
        	for (int i = 0; i < glen; i++)
		{
			if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
			{
				ok = false;
				num = -1;
				while (!ok)
				{
					num = getRng()->getInt(0,256);
					ok = ((num > 0) && (num < 256));
				}
				m_pop[to][i] = num;
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
	bool ok;
	int num;
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
		// Mutation
        	for (int i = 0; i < glen; i++)
		{
			if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
			{
				ok = false;
				num = -1;
				while (!ok)
				{
					num = getRng()->getInt(0,256);
					ok = ((num > 0) && (num < 256));
				}
				to[i] = num;
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
	bool ok;
	int num;
	
	currVal = -1;
	while (!modified)
	{
		ok = false;
		num = -1;
		while (!ok)
		{
			num = getRng()->getInt(0,256);
			ok = ((num > 0) && (num < 256));
		}
		currVal = num;
		modified = (currVal != oldVal);
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

double MyEvoAlgo::evaluate(const QVector<int> ind, int& steps, int& completed)
{
	double fitness = 0.0;
	int i;
	const double len = m_cartPole.getTaskLength();
	double fit;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	// Reset <completed>
	completed = 0;
	// Evaluate network
	for (i = 0; i < m_taskEpisodes; i++)
	{
		evonet->resetNet();
		fit = m_cartPole.evalNet(evonet, i);
		fitness += (fit / len);
		steps += int(fit);
		if (m_cartPole.generalization_test)
		{
			// Store how many generalization episodes have been solved
			if ((fit / len) == 1.0)
				completed++;
		}
	}
	fitness /= m_taskEpisodes;
	//printf("Fitness %.3f\n", fitness);
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
		sprintf(filename, "G%dS%d.gen", gen, currentSeed);
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

void MyEvoAlgo::saveBestgInd(QVector<int> ind)
{
	FILE *fp;
    char filename[64];
	int j;

	sprintf(filename, "bestgS%d.gen", currentSeed);
    if ((fp = fopen(filename, "w+")) == NULL) 
	{
        farsa::Logger::error(QString("Cannot open file %1").arg(filename));
    } 
	else {
		//now saving
		fprintf(fp, "**NET : %d_%d_%d.wts\n", cgen, 0, 0);
		fprintf(fp, "DYNAMICAL NN\n");
		for (j = 0; j < glen; j++)
		{
			fprintf(fp, "%d\n", ind[j]);
		}
		fprintf(fp, "END\n");
	}
}

void MyEvoAlgo::saveGStat(double fit, double gfit, int completed)
{
	m_gstatBuf += QString("%1 %2 %3\n").arg(fit).arg(gfit).arg(completed);
}

void MyEvoAlgo::saveEvalStats(int steps, double fit, double gfit, int completed, double fixfit)
{
	if ((m_taskType != 1) && (m_taskType != 3))
		m_evalstatBuf += QString("%1 %2 %3 %4\n").arg(steps).arg(fit).arg(gfit).arg(completed);
	else
		m_evalstatBuf += QString("%1 %2 %3 %4 %5\n").arg(steps).arg(fit).arg(fixfit).arg(gfit).arg(completed);
}

void MyEvoAlgo::invalidFill(QVector<int>& ind)
{
	for (int i = 0; i < ind.size(); i++)
	{
		ind[i] = -1;
	}
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolve()
{
	farsa::Logger::warning("Evolving Stochastic Hill-Climber from MyEvoAlgo class");

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
	int steps;
	double bestfit;
	int bestid;
	int completed; // Episodes solved during generalization
	int bestcompleted;
	// The following three variables refer to random episodes
	double bestfixedfit;
	int bestfixedid;
	double fixedfit;
	double tmpfit;
        double pfit;
        double cfit;
	
	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	initPop();
		// Set fbest to a very low value
		this->fbest = -99999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
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
		m_gstatBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpgStatBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		solved = false;
		finish = false;
		bestfit = -9999.0;
		bestid = -1;
		bestcompleted = -9999;
		bestfixedfit = -9999.0;
		while (cstep < m_nevalsteps)// && !solved)
		//while (!solved && !finish)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				performMutation(id, popSize + id, 1);
			}
			for(id = 0; id < popSize; id++) 
			{
				// Test the individual and its child
				// Parents are tested only at the first generation
				if ((gn == startGeneration) || (m_fitNoise > 0.0) || ((m_taskType == 1) || (m_taskType == 3)))
				{
					m_cartPole.generalization_test = false;
					m_taskEpisodes = EVAL_EPS;
					tfitness[id] = 0.0;
					ntfitness[id] = 0.0;
					completed = 0;
					// evaluate the parent
					steps = 0;
					fit = evaluate(m_pop[id], steps, completed);
					// update its fitness
					tfitness[id] += fit;
					ntfitness[id]++;
					if (isStopped()) 
					{ // stop evolution
						return;
					}
					cstep += steps;
					if ((m_taskType == 0) || (m_taskType == 2))
					{
						if ((fit == 1.0) && !solved)
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
					}
					// Check whether or not individual is better than current best
					if (fit > bestfit)
					{
						bestfit = fit;
						bestid = id;
					}
					pfit = fit;
				}
				m_cartPole.generalization_test = false;
				m_taskEpisodes = EVAL_EPS;
				tfitness[popSize + id] = 0.0;
				ntfitness[popSize + id] = 0.0;
				completed = 0;
				// evaluate the child
				steps = 0;
				fit = evaluate(m_pop[popSize + id], steps, completed);
				// update its fitness
				tfitness[popSize + id] += fit;
				ntfitness[popSize + id]++;
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((m_taskType == 0) || (m_taskType == 2))
				{
					if ((fit == 1.0) && !solved)
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
				}
				// Check whether or not individual is better than current best
				if (fit > bestfit)
				{
					bestfit = fit;
					bestid = (popSize + id);
				}
				cfit = fit;
				// Apply noise to both parent fitness and offspring fitness
		                pfit = tfitness[id] / ntfitness[id]; // Useless during first generation
		                pfit *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
				cfit = tfitness[popSize + id] / ntfitness[popSize + id]; // Useless
		                cfit *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
		                if (cfit >= pfit)
				{
					performMutation(popSize + id, id, 0);
	       				tfitness[id] = tfitness[popSize + id];
					ntfitness[id] = ntfitness[popSize + id];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 
			}
			
			// If task type corresponds to random initial states condition, evaluate current best on fixed initial states condition
			if ((m_taskType == 1) || (m_taskType == 3))
			{
				// Evaluate current bestid on the fixed episodes
		                m_cartPole.generalization_test = false;
				// Set fixed state flag to true
		                m_cartPole.fixed_state = true;
		                // Set the number of test episodes
				m_taskEpisodes = EVAL_EPS;
				// Evaluate
				steps = 0;
				fixedfit = evaluate(m_pop[bestid], steps, completed);
				cstep += steps;
				// Check whether or not the individual generalizes better than current bestg
				if (fixedfit > bestfixedfit)
				{
					bestfixedfit = fixedfit;
					bestfixedid = bestid; // Useless
				}
				// Reset fixed state flag to false
				m_cartPole.fixed_state = false;
				// Check if the fitness is 1.0
				if ((fixedfit == 1.0) && !solved)
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
			}
			else
			{
				// The best fixed fitness matches the best fitness
				bestfixedfit = bestfit;
			}
			// Compute and save statistics
		    	computeFStat2();
		    	saveFStat();

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		cgen++;

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Best gfitness = %5 - Generalization episodes solved = %6 - Fixed episodes fit = %7 - Evaluation steps = %8").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmax).arg(bestfit).arg(-9999.0).arg(bestcompleted).arg(bestfixedfit).arg(cstep));
            		fflush(stdout);

			saveEvalStats(cstep, bestfit, -9999.0, bestcompleted, bestfixedfit);

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
				finish = true;

			//if (((gn + 1) % m_resetBufNumGen) == 0)
			if (finish)//if (solved || finish)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
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
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}

				if (m_calcCorr)
				{
					// Store population's fitness in order to compute the correlation between fitness and other metrics (i.e., population size, mutation rate, etc.)
					char corrfname[1024];
					sprintf(corrfname, "corrPopFitS%d.txt", currentSeed);
					FILE* corrfp = fopen(corrfname, "w");
					if (corrfp != NULL)
					{
						for (id = 0; id < popSize; id++)
						{
							if (ntfitness[id] > 0.0)
								fprintf(corrfp, "%lf\n", tfitness[id] / ntfitness[id]);
							else
								fprintf(corrfp, "0.0\n");
						}
						fclose(corrfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(corrfname));
					}
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

MyEvoAlgo::CartPole::CartPole()
	: MUP(0.000002)
	, MUC(0.0005)
	, GRAVITY(-9.8)
	, MASSCART(1.0)
	, MASSPOLE_1(0.1)
	, LENGTH_1(0.5)
	, FORCE_MAG(10.0)
	, TAU(0.01)
	, one_degree(0.0174532)
	, six_degrees(0.1047192)
	, twelve_degrees(0.2094384)
	, fifteen_degrees(0.2617993)
	, thirty_six_degrees(0.628329)
	, fifty_degrees(0.87266)
	, fixed_state(true)
	, rng(1)
{
}

MyEvoAlgo::CartPole::~CartPole()
{
}

void MyEvoAlgo::CartPole::initializeCartPole(bool velocity, bool fixedState){
  maxFitness = 100000;

  MARKOV=velocity;
  fixed_state = fixedState;
  stampOuput = "";

  LENGTH_2 = 0.05;
  MASSPOLE_2 = 0.01;

  // MyEvoAlgo::CartPole::reset() which is called here
}

void MyEvoAlgo::CartPole::setCartPole(const int s)
{
	// Set the seed
	setSeed(s);
	// Initialize generalization states
	initGStates();
}

farsa::RandomGenerator* MyEvoAlgo::CartPole::getRandomGenerator()
{
	return &rng;
}

void MyEvoAlgo::CartPole::setSeed(const int s)
{
	rng.setSeed(s);
}

double MyEvoAlgo::CartPole::evalNet(farsa::Evonet *evonet, int episode)
{
	int steps=0;
	double input[NUM_INPUTS];
	double output;

	int nmarkovmax;

	double nmarkov_fitness;

	double jiggletotal; //total jiggle in last 100
	int count;  //step counter
	int ninputs; // number of inputs

	if (nmarkov_long) nmarkovmax=100000;
	else if (generalization_test) nmarkovmax=1000;
	else nmarkovmax=1000;

	init(episode);

	if (MARKOV)
	{ // Velocities provided as inputs
		//printf("Episode %d - Initial state: %lf %lf %lf %lf %lf %lf\n", episode, state[0], state[1], state[2], state[3], state[4], state[5]);
		//while (steps++ < maxFitness)
		while (steps < nmarkovmax)
		{
			input[0] = state[0] / 4.8;
			input[1] = state[1] / 2;
			input[2] = state[2] / 0.52;
			input[3] = state[3] / 2;
			input[4] = state[4] / 0.52;
			input[5] = state[5] / 2;
			input[6] = 0.5;
			
			for (int j = 0; j < NUM_INPUTS; j++)
				evonet->setInput(j, input[j]);
			//Activate the net
			evonet->updateNet();
			// Get the output
			output = evonet->getOutput(0);
			//printf("steps %d - output %lf\n", steps, output);
			//Sleep(1000);
			// Update steps
			steps++;
			// Perform the action
			performAction(output,steps);
			// Check whether or not the pole(s) is(are) fallen
			if (outsideBounds())	// if failure
			{
				//printf("Out of bounds! cart state %.3f limit 2.4 - pole angle %.3f limit %.3f\n", state[0], state[2], thirty_six_degrees);
				//Sleep(1000);
				break;			// stop it now
			}
		}
		//printf("Completed %d steps!!!\n", steps);
		//Sleep(1000);
		return ((double) steps);
	}
	else
	{  //NON MARKOV CASE
		while (steps < nmarkovmax)
		{
			input[0] = state[0] / 4.8;
			input[1] = state[2] / 0.52;
			input[2] = state[4] / 0.52;
			input[3] = 0.5;
			ninputs = 4;
			
			for (int j = 0; j < ninputs; j++)
				evonet->setInput(j, input[j]);
			//Activate the net
			evonet->updateNet();
			// Get the output
			output=evonet->getOutput(0);
			// Update steps
			steps++;
			// Perform the action
			performAction(output,steps);
			// Check whether or not the pole(s) is(are) fallen
			if (outsideBounds())	// if failure
				break;			// stop it now

			if (nmarkov_long && (outsideBounds()))	// if failure
				break;			// stop it now
		}

		//If we are generalizing we just need to balance it a while
		/*if (generalization_test)
		{
			return (double) steps;//balanced_sum;
		}*/

		//Sum last 100
		/*if ((steps>100)&&(!nmarkov_long)) 
		{
			jiggletotal=0;
			//Adjust for array bounds and count
			for (count=steps-99-2;count<=steps-2;count++)
				jiggletotal+=jigglestep[count];
		}*/

		/*if (!nmarkov_long) 
		{
			if (balanced_sum > 100)
			{
				nmarkov_fitness=((0.1*(((double) balanced_sum)/1000.0))+(0.9*(0.75/(jiggletotal))));
			}
			else
			{
				nmarkov_fitness=((double) balanced_sum);
			}
			return nmarkov_fitness;
		}
		else*/
		{
			return (double) steps;
		}
	}
}

void MyEvoAlgo::CartPole::init(int episode)
{
	int j;

	if (!MARKOV) {
		//Clear all fitness records
		cartpos_sum=0.0;
		cartv_sum=0.0;
		polepos_sum=0.0;
		polev_sum=0.0;
	}

	balanced_sum=0; //Always count # balanced

	last_hundred=false;
	
	if (!generalization_test)
	{
		// Evaluation phase
		if (fixed_state)
		{
			// Fixed initial states
			if (episode < 0)
			{
				state[0] = state[1] = state[3] = state[4] = state[5] = 0;
				state[2] = 0.07; // Approximately 4°
				//farsa::Logger::info(QString("state: %1 %2 %3 %4 0.5").arg(state[0]).arg(state[1]).arg(state[2]).arg(state[3]));
			}
			else
			{
				if (episode == 0)
				{
					state[0] = -1.944;
					state[1] = state[2] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 1)
				{
					state[0] = 1.944;
					state[1] = state[2] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 2)
				{
					state[1] = -1.215;
					state[0] = state[2] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 3)
				{
					state[1] = 1.215;
					state[0] = state[2] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 4)
				{
					state[2] = -0.10472;
					state[0] = state[1] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 5)
				{
					state[2] = 0.10472;
					state[0] = state[1] = state[3] = state[4] = state[5] = 0;
				}
				if (episode == 6)
				{
					state[3] = -0.135088;
					state[0] = state[1] = state[2] = state[4] = state[5] = 0;
				}
				if (episode == 7)
				{
					state[3] = 0.135088;
					state[0] = state[1] = state[2] = state[4] = state[5] = 0;
				}
			}
			//Sleep(100000);
		}
		else
		{
			// Random initial states
			state[0] = getRandomGenerator()->getDouble(-1.944,1.944);
			state[1] = getRandomGenerator()->getDouble(-1.215,1.215);
			state[2] = getRandomGenerator()->getDouble(-0.10472,0.10472);
			state[3] = getRandomGenerator()->getDouble(-0.135088,0.135088);
			state[4] = getRandomGenerator()->getDouble(-0.10472,0.10472);
			state[5] = getRandomGenerator()->getDouble(-0.135088,0.135088);
		}
	}
	else
	{
		// Consistency check about the value of variable episode
		if ((episode < 0) || (episode >= TEST_EPS))
		{
			farsa::Logger::error(QString("Invalid episode %1! Test episodes are %2! Check your code...").arg(episode).arg(TEST_EPS));
			exit(-1);
		}
		// Generalization
		for (j = 0; j < NUM_INPUTS - 1; j++)
			state[j] = gStates[episode][j];
	
	}
}

double MyEvoAlgo::CartPole::getTaskLength()
{
	if (nmarkov_long)
		return 100000.0;
	else
		return 1000.0;
}

void MyEvoAlgo::CartPole::performAction(double output, int stepnum)
{

	int i;
	double  dydx[6];

	const bool RK4 = true; //Set to Runge-Kutta 4th order integration method
	const double EULER_TAU = TAU / 4;

    	if(testCart)
        	stampOuput += QString("%1 %2 %3 \n").arg(state[0]).arg(state[2]).arg(state[4]);
  
	/*--- Apply action to the simulated cart-pole ---*/
	if(RK4)
	{
		for(i=0;i<2;++i)
		{
			dydx[0] = state[1];
			dydx[2] = state[3];
			dydx[4] = state[5];
			step(output,state,dydx);
			rk4(output,state,dydx,state);
		}
	}
	else
	{
		for(i=0;i<8;++i)
		{
			dydx[0] = state[1];
			dydx[2] = state[3];
			dydx[4] = state[5];
			step(output,state,dydx);
			state[0] += EULER_TAU * dydx[0];
			state[1] += EULER_TAU * dydx[1];
			state[2] += EULER_TAU * dydx[2];
			state[3] += EULER_TAU * dydx[3];
			state[4] += EULER_TAU * dydx[4];
			state[5] += EULER_TAU * dydx[5];
		}
	}

	//Record this state
	cartpos_sum += fabs(state[0]);
	cartv_sum += fabs(state[1]);
	polepos_sum += fabs(state[2]);
	polev_sum += fabs(state[3]);
	if (stepnum <= 1000)
		jigglestep[stepnum-1]=fabs(state[0])+fabs(state[1])+fabs(state[2])+fabs(state[3]);

	if (!(outsideBounds()))
		++balanced_sum;
}

void MyEvoAlgo::CartPole::step(double action, double *st, double *derivs)
{
    	double force,costheta_1,costheta_2,sintheta_1,sintheta_2,
          gsintheta_1,gsintheta_2,temp_1,temp_2,ml_1,ml_2,fi_1,fi_2,mi_1,mi_2;
	
	force =  (action - 0.5) * FORCE_MAG * 2;
	costheta_1 = cos(st[2]);
	sintheta_1 = sin(st[2]);
	gsintheta_1 = GRAVITY * sintheta_1;
	costheta_2 = cos(st[4]);
	sintheta_2 = sin(st[4]);
	gsintheta_2 = GRAVITY * sintheta_2;

	ml_1 = LENGTH_1 * MASSPOLE_1;
	ml_2 = LENGTH_2 * MASSPOLE_2;
	temp_1 = MUP * st[3] / ml_1;
	temp_2 = MUP * st[5] / ml_2;
	fi_1 = (ml_1 * st[3] * st[3] * sintheta_1) +
			(0.75 * MASSPOLE_1 * costheta_1 * (temp_1 + gsintheta_1));
	fi_2 = (ml_2 * st[5] * st[5] * sintheta_2) +
			(0.75 * MASSPOLE_2 * costheta_2 * (temp_2 + gsintheta_2));
	mi_1 = MASSPOLE_1 * (1 - (0.75 * costheta_1 * costheta_1));
	mi_2 = MASSPOLE_2 * (1 - (0.75 * costheta_2 * costheta_2));

	derivs[1] = (force + fi_1 + fi_2)
					/ (mi_1 + mi_2 + MASSCART);

	derivs[3] = -0.75 * (derivs[1] * costheta_1 + gsintheta_1 + temp_1)
					/ LENGTH_1;
	derivs[5] = -0.75 * (derivs[1] * costheta_2 + gsintheta_2 + temp_2)
					/ LENGTH_2;
}

void MyEvoAlgo::CartPole::rk4(double f, double y[], double dydx[], double yout[])
{

    int i;

    double hh,h6,dym[6],dyt[6],yt[6];


    hh=TAU*0.5;
    h6=TAU/6.0;
    for (i=0;i<=5;i++) yt[i]=y[i]+hh*dydx[i];
    step(f,yt,dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i=0;i<=5;i++) yt[i]=y[i]+hh*dyt[i];
    step(f,yt,dym);
    dym[0] = yt[1];
    dym[2] = yt[3];
    dym[4] = yt[5];
    for (i=0;i<=5;i++) {
        yt[i]=y[i]+TAU*dym[i];
        dym[i] += dyt[i];
    }
    step(f,yt,dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i=0;i<=5;i++)
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

bool MyEvoAlgo::CartPole::outsideBounds()
{
	const double failureAngle = thirty_six_degrees;

	return
		(state[0] < -2.4             ||
		state[0] > 2.4               ||
		state[2] < -failureAngle     ||
		state[2] > failureAngle      ||
		state[4] < -failureAngle     ||
		state[4] > failureAngle);
}

void MyEvoAlgo::CartPole::initGStates()
{
	int i;
	for (i = 0; i < TEST_EPS; i++)
	{
		gStates[i][0] = getRandomGenerator()->getDouble(-1.944,1.944);
		gStates[i][1] = getRandomGenerator()->getDouble(-1.215,1.215);
		gStates[i][2] = getRandomGenerator()->getDouble(-0.10472,0.10472);
		gStates[i][3] = getRandomGenerator()->getDouble(-0.135088,0.135088);
		gStates[i][4] = getRandomGenerator()->getDouble(-0.10472,0.10472);
		gStates[i][5] = getRandomGenerator()->getDouble(-0.135088,0.135088);
	}
}
