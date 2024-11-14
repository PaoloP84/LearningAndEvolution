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

#define NTASKS 4
#define EVAL_EPS 8
#define TEST_EPS 1000

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
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluate(const QVector<int> ind, int& steps, int& completed);
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolve();

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
	 * \brief Saves the best generalizing individual
	 *
	 */
	void saveBestgInd(QVector<int> ind);
	/**
	 * \brief Saves generalization stats
	 *
	 */
	void saveGStat(double fit, double gfit, int completed);
	/**
	 * \brief Saves evaluation stats
	 *
	 */
	void saveEvalStats(int steps, double fit, double gfit, int completed, double fixfit);
	/**
	 * \brief Fills a QVector with -1 values (invalid)
	 *
	 * \param ind The QVector
	 */
	void invalidFill(QVector<int>& ind);
	farsa::RandomGenerator* getRng();

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
	//! Task type:
	//!  - 0: double-pole with velocity (fixed initial state(s))
	//!  - 1: double-pole with velocity (random initial states)
	//!  - 2: double-pole without velocity (fixed initial state(s))
	//!  - 3: double-pole without velocity (random initial states)
	int m_taskType;
	//! Task episodes
	int m_taskEpisodes;
	//! Generalization flag
	bool m_generalize;
	//! Random generator
	farsa::RandomGenerator m_rng;
	//! Correlation computation
	bool m_calcCorr;
	
	//! The class implementing the Cart-Pole
	class CartPole {
	public:
		CartPole();
		~CartPole();
		void initializeCartPole(bool velocity, bool fixedState);
		void setCartPole(const int s);
		virtual double evalNet(farsa::Evonet* evonet, int thresh);
		/**
		 * \brief Returns the random number generator
		 *
		 * \return The random number generator
		 */
		farsa::RandomGenerator* getRandomGenerator();
		void setSeed(const int s);
		double getTaskLength();
		double maxFitness;
		bool MARKOV;
		bool fixed_state;
		
		//! Local random number generator
		farsa::RandomGenerator rng;
		
		bool last_hundred;
		bool nmarkov_long;  //Flag that we are looking at the champ
		bool generalization_test;  //Flag we are testing champ's generalization
		bool testCart;
		QString stampOuput;
		
		double state[6];
		
		double jigglestep[1000];
	
	protected:
		virtual void init(int type);
	
	private:
	
		void performAction(double output, int stepnum);
		void step(double action, double *state, double *derivs);
		void rk4(double f, double y[], double dydx[], double yout[]);
		bool outsideBounds();

		void initGStates();
		
		const static int NUM_INPUTS=7;
		//const static double MUP = 0.000002;
		//const static double MUC = 0.0005;
		//const static double GRAVITY= -9.8;
		//const static double MASSCART= 1.0;
		//const static double MASSPOLE_1= 0.1;
		
		//const static double LENGTH_1= 0.5;		  /* actually half the pole's length */
		
		//const static double FORCE_MAG= 10.0;
		//const static double TAU= 0.01;		  //seconds between state updates
		
		//const static double one_degree= 0.0174532;	/* 2pi/360 */
		//const static double six_degrees= 0.1047192;
		//const static double twelve_degrees= 0.2094384;
		//const static double fifteen_degrees= 0.2617993;
		//const static double thirty_six_degrees= 0.628329;
		//const static double fifty_degrees= 0.87266;

		const double MUP;
		const double MUC;
		const double GRAVITY;
		const double MASSCART;
		const double MASSPOLE_1;
		
		const double LENGTH_1;		  /* actually half the pole's length */
		
		const double FORCE_MAG;
		const double TAU;		  //seconds between state updates
		
		const double one_degree;	/* 2pi/360 */
		const double six_degrees;
		const double twelve_degrees;
		const double fifteen_degrees;
		const double thirty_six_degrees;
		const double fifty_degrees;
		
		double LENGTH_2;
		double MASSPOLE_2;
		
		//Queues used for Gruau's fitness which damps oscillations
		int balanced_sum;
		double cartpos_sum;
		double cartv_sum;
		double polepos_sum;
		double polev_sum;

		// Generalization states
		double gStates[TEST_EPS][NUM_INPUTS - 1];
	};
	
	// The cart pole used for evaluating individuals
	CartPole m_cartPole;
};

#endif
