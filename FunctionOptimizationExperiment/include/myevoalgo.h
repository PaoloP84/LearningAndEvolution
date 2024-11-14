/********************************************************************************
 *  FARSA - Total99                                                             *
 *  Copyright (C) 2012-2013 Gianluca Massera <emmegian@yahoo.it>                *
 *                                                                              *
 *  This program is free software; you can redistribute it and/or modify        *
 *  it under the terms of the GNU General Public License as published by        *
 *  the Free Software Foundation; either version 2 of the License, or           *
 *  (at your option) any later version.                                         *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful,             *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 *  GNU General Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License           *
 *  along with this program; if not, write to the Free Software                 *
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
 ********************************************************************************/

#ifndef MYEVOALGO_H
#define MYEVOALGO_H

#include "farsaplugin.h"

#include "evoga.h"
#include "logger.h"
#include "evodataviewer.h"
#include "randomgenerator.h"
#include "evorobotexperiment.h"
#include "evonet.h"
#include <parametersettable.h>
#include <configurationhelper.h>
#include <configurationparameters.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <deque>
#include <iomanip>
#include <sstream>
#include <string>

#include <QFile>
#include <QTime>
#include <QVector>
#include <QtAlgorithms>
#include <QThreadPool>
#include <QtConcurrentMap>

/**
 * \brief My evolutionary algorithm. It inherits from Evoga
 *        and runs a modified version of the steadyState
 *        algorithm (with floating-point genes).
 */
class FARSA_PLUGIN_API MyEvoAlgo : public farsa::Evoga
{
	Q_OBJECT
	FARSA_REGISTER_CLASS(Evoga)

public:
     /**
	 * \brief Constructor
	 */
	 MyEvoAlgo();


     /**
	 * \brief Destructor
	 */     
    ~MyEvoAlgo();
	
	/**
	 * \brief Configures the object using a ConfigurationParameters object
	 *
	 * \param params the configuration parameters object with parameters to
	 *               use
	 * \param prefix the prefix to use to access the object configuration
	 *               parameters. This is guaranteed to end with the
	 *               separator character when called by the factory, so you
	 *               don't need to add one
	 */
	virtual void configure(farsa::ConfigurationParameters& params, QString prefix);

	//! Save
	virtual void save(farsa::ConfigurationParameters& params, QString prefix)
	{
		// Calling parent function
		farsa::Evoga::save(params, prefix);
		farsa::Logger::warning("NOT IMPLEMENTED (MyEvoAlgo::save)");
		abort();
	}

	//! Describe
	static void describe(QString type);
	
	void postConfigureInitialization();
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolve();

	void evolveHC();

	//! Evaluation function
	double evaluate(QVector<double> x);

	/**
	 * \brief It invokes either evolve() or test() method independently of
	 *        the chosen type. This depends on the <test> flag.
	 */
	virtual void evolveAllReplicas();

	void setSeed(const int seed);

private:
	/**
	 * \brief Initialises the population
	 */
	void initPop();
	/**
	 * \brief Resets the population (new size = 0)
	 */
	void resetPop();
	/**
	 * \brief Mutates the gene
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 * \param mut the mutation rate
	 */
	void performMutation(int from, int to, int mut);
	/**
	 * \brief Mutates the gene (does not consider the population)
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 * \param mut the mutation rate
	 */
	void immediateMutation(const QVector<int> from, QVector<int>& to, int mut);
	/**
	 * \brief Mutates one gene
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 */
	void neutralMutation(const QVector<int> from, QVector<int>& to);
	/**
	 * \brief Compute fitness statistics (min, max, avg)
	 *
	 */
	void computeFStat();
	/**
	 * \brief Loads a genotype
	 *
	 * \param fp the file descriptor
	 * \param ind the individual to be loaded
	 */
	void loadGen(FILE *fp, int ind);
	/**
	 * \brief Loads all the genotypes
	 *
	 * \param gen the generation
	 * \param filew the filename
	 * \return the number of loaded individuals
	 */
	int loadAllGen(int gen, char* filew);
	/**
	 * \brief Saves the best individual
	 *
	 */
	void saveBestInd();
	/**
	 * \brief Saves a genotype
	 *
	 * \param fp The file pointer
	 * \param ind The individual to be saved
	 */
	void saveGen(FILE *fp, int ind);
	/**
	 * \brief Saves all the genotypes
	 *
	 */
	void saveAllGen();
	/**
	 * \brief Saves fitness statistics (min, max, avg)
	 *
	 */
	void saveFStat();
	/**
	 * \brief Saves evaluation stats
	 *
	 */
	void saveEvalStats(int steps, double fit);
	/**
	 * \brief Fills a QVector with -1 values (invalid)
	 *
	 * \param ind The QVector
	 */
	void invalidFill(QVector<int>& ind);
	/**
	 * \brief Annealing routine
	 *
	 * \param original The original individual (before the annealing routine)
	 * \param novel The new individual (after the annealing routine)
	 * \return The fitness of the novel individual
	 */
	double annealing(const QVector<int> original, double fitness, double actualFitness, QVector<int>& novel, int& steps, bool& solved);

	//! Rastrigin function
	double rastrigin(QVector<double> x);
	//! Rosenbrock function
	double rosenbrock(QVector<double> x);
	//! Sphere function
	double sphere(QVector<double> x);

	//! Convert geno to pheno
	void genoToPheno(const QVector<int> geno, QVector<double>& pheno);

	farsa::RandomGenerator* getRng();

	void hamming();

	//! Population
	QVector< QVector<int> > m_pop;
	//! Number of evaluation steps
	int m_nevalsteps;
	//! Mutation rate
	double m_mutationRate;
	
	//! Buffers for saving best genotypes and best fitness
	QString m_bestGenBuf;
	QString m_statBuf;
	QString m_gstatBuf;
	QString m_evalstatBuf;
	
	//! Flags whether or not the information about the best individual must be saved
	bool m_saveBestGenInfo;
	//! Interval for saving information about the best individual
	int m_bestGenInterval;
	//! Noise on the individual's fitness (i.e., stochasticity)
	double m_fitNoise;
	//! Flags whether or not the annealing must be run
	bool m_annealing;
	//! Number of annealing generations
	int m_numAnnealingGenerations;
	//! Flags whether or not randomisation must be performed in order to avoid the replacement of the same worst individuals
	bool m_randomOrder;
	//! Algorithm (0: SSS or SSSHC; 1: SSSHC)
	int m_algo;
	//! Function id (0: Rastrigin; 1: Rosenbrock; 2: Sphere)
	int m_funct;
	//! Size of vector (complexity of the problem)
	int m_n;
	//! Random generator
	farsa::RandomGenerator m_rng;
};

#endif
