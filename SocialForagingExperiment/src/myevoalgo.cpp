#include "myevoalgo.h"
#include "evonet.h"
#include <cmath>
#include <QTextStream>
#include <QFile>

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
	, m_annealing(false)
	, m_numAnnealingGenerations(0)
	, m_randomOrder(false)
	, m_algo(0)
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
	m_algo = farsa::ConfigurationHelper::getInt(params, prefix + "algo", m_algo);
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
	d.describeBool("annealing").def(false).help("Flags whether or not the annealing must be run");
	d.describeInt("numAnnealingGenerations").limits(0,INT_MAX).def(0).help("Number of annealing generations");
	d.describeBool("randomOrder").def(false).help("Flags whether or not randomisation must be performed in order to avoid the replacement of the same worst individuals");
}

void MyEvoAlgo::postConfigureInitialization()
{
	// Calling parent function
	farsa::Evoga::postConfigureInitialization();
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
	// this class will always call my evolutionary algorithm.
	if (m_algo == 0)
		evolve();
	else if (m_algo == 1)
		evolveHC();
	else
		hamming();
}
void MyEvoAlgo::setSeed(int s)
{
	currentSeed = s;
	m_rng.setSeed(s);
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

double MyEvoAlgo::annealing(const QVector<int> original, double fitness, double actualFitness, QVector<int>& novel, int& steps, bool& solved)
{
	double novelFit;
	double currFit;
	// Compute all the bit strings of length <m_input_size>
    	QVector<int> orig = original;
	QVector<int> evalInd(glen);
	int annGn = 0;
	double fit = 0.0;
	double storedFit = actualFitness;
	int nsteps;
	int dummy;
	// Reset temporary individual
	invalidFill(evalInd);
	while (annGn < m_numAnnealingGenerations)
	{
		// We must generate new noise at each iteration to the original fitness
		currFit = storedFit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
		// Mutate the individual
		immediateMutation(orig, evalInd, 1);
		// Evaluate
		nsteps = exp->getNSteps();
		exp->setNetParameters(evalInd);
		exp->doAllTrialsForIndividual(0);
		fit = exp->getFitness();
		steps += nsteps;
		// Add noise
		double noisyFit = fit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
		if (noisyFit >= currFit)
		{
			immediateMutation(evalInd, orig, 0);
			storedFit = fit;
			currFit = noisyFit;
		}
		annGn++;
	}
	novel = orig;
	novelFit = storedFit;
	return novelFit;
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

void MyEvoAlgo::saveEvalStats(int steps, double fit)
{
	m_evalstatBuf += QString("%1 %2\n").arg(steps).arg(fit);
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
	farsa::Logger::info("EVOLUTION: steady state with custom evaluation function");
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
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		// generations
		gn = startGeneration;
		cstep = 0;
		solved = false;
		finish = false;
		bestfit = -9999.0;
		bestid = -1;
		while (cstep < m_nevalsteps)// && !solved)
		//while (!solved && !finish)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	vector<int> identity;
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				identity.push_back(id);
				performMutation(id, popSize + id, 1);
			}
			for(id = 0; id < popSize; id++) 
			{
				// Test the individual and its child
				// Parents are tested only at the first generation
				if ((gn == startGeneration) || (m_fitNoise > 0.0))
				{
					tfitness[id] = 0.0;
					ntfitness[id] = 0.0;
					// evaluate the parent
					steps = exp->getNSteps();
					exp->setNetParameters(m_pop[id]);
					exp->doAllTrialsForIndividual(id);
					fit = exp->getFitness();
					// update its fitness
					tfitness[id] = fit;
					ntfitness[id] = 1;
					if (isStopped()) 
					{ // stop evolution
						return;
					}
					cstep += steps;
					// Check whether or not individual is better than current best
					if (fit > bestfit)
					{
						bestfit = fit;
						bestid = id;
					}
				}
				tfitness[popSize + id] = 0.0;
				ntfitness[popSize + id] = 0.0;
				// evaluate the child
				steps = exp->getNSteps();
				exp->setNetParameters(m_pop[popSize + id]);
				exp->doAllTrialsForIndividual(popSize + id);
				fit = exp->getFitness();
				// update its fitness
				tfitness[popSize + id] = fit;
				ntfitness[popSize + id] = 1;
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				// Check whether or not individual is better than current best
				if (fit > bestfit)
				{
					bestfit = fit;
					bestid = (popSize + id);
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
                		parents[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
            		}

            		for (int i = 0; i < popSize; i++)
			{
				children[i].fitness = tfitness[popSize + identity[i]] / ntfitness[popSize + identity[i]];
                		children[i].id = popSize + identity[i];
                		children[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
            		}
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
					QVector<int> original(glen);
					immediateMutation(m_pop[index], original, 0);
					double fitness = selectedFit[i];
					double actualFitness = tfitness[index] / ntfitness[index];
					QVector<int> novel(glen);
					double novelFit = annealing(original, fitness, actualFitness, novel, cstep, solved);
					// Update the individual
					immediateMutation(novel, m_pop[index], 0);
					// Update structure data
					selectedFit[i] = novelFit;
					tfitness[selected[i]] = selectedFit[i];
					// Check whether or not the individual returned by the annealing routine is better than current best
					if (novelFit > bestfit)
					{
						bestfit = novelFit;
						bestid = index;
					}
				}
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

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		cgen++;

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation: max %3 avg %4 min %5 - Best fitness = %6 - Evaluation steps = %7").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmax).arg(faverage).arg(fmin).arg(bestfit).arg(cstep));
            		fflush(stdout);

			saveEvalStats(cstep, bestfit);

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
				finish = true;

			if (finish)
			{
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


/*
 * Main function of the Hill-Climber Algorithm
 */
void MyEvoAlgo::evolveHC()
{
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
	double pfit;
	double cfit;
	farsa::Logger::info("EVOLUTION: Hill-Climber");
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
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		// generations
		gn = startGeneration;
		cstep = 0;
		solved = false;
		finish = false;
		bestfit = -9999.0;
		bestid = -1;
		while (cstep < m_nevalsteps)// && !solved)
		//while (!solved && !finish)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	vector<int> identity;
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				identity.push_back(id);
				performMutation(id, popSize + id, 1);
			}
			for(id = 0; id < popSize; id++) 
			{
				// Test the individual and its child
				// Parents are tested only at the first generation
				if ((gn == startGeneration) || (m_fitNoise > 0.0))
				{
					tfitness[id] = 0.0;
					ntfitness[id] = 0.0;
					// evaluate the parent
					steps = exp->getNSteps();
					exp->setNetParameters(m_pop[id]);
					exp->doAllTrialsForIndividual(id);
					fit = exp->getFitness();
					// update its fitness
					tfitness[id] = fit;
					ntfitness[id] = 1;
					if (isStopped()) 
					{ // stop evolution
						return;
					}
					cstep += steps;
					// Check whether or not individual is better than current best
					if (fit > bestfit)
					{
						bestfit = fit;
						bestid = id;
					}
				}
				tfitness[popSize + id] = 0.0;
				ntfitness[popSize + id] = 0.0;
				// evaluate the child
				steps = exp->getNSteps();
				exp->setNetParameters(m_pop[popSize + id]);
				exp->doAllTrialsForIndividual(popSize + id);
				fit = exp->getFitness();
				// update its fitness
				tfitness[popSize + id] = fit;
				ntfitness[popSize + id] = 1;
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				// Check whether or not individual is better than current best
				if (fit > bestfit)
				{
					bestfit = fit;
					bestid = (popSize + id);
				}
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
			
			// Compute and save statistics
		    	computeFStat2();
		    	saveFStat();

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		cgen++;

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation: max %3 avg %4 min %5 - Best fitness = %6 - Evaluation steps = %7").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmax).arg(faverage).arg(fmin).arg(bestfit).arg(cstep));
            		fflush(stdout);

			saveEvalStats(cstep, bestfit);

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
				finish = true;

			//if (((gn + 1) % m_resetBufNumGen) == 0)
			if (finish)//if (solved || finish)
			{
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

