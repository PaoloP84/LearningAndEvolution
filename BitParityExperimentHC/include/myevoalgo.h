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

	struct BinaryTree {
      BinaryTree *left, *right;
      int data;
      BinaryTree(int val) : left(NULL), right(NULL), data(val) { }
    };

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
	double evaluate(const QVector<int> ind);
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolve();
	/**
	 * \brief It invokes evolve method independently of the chosen type
	 */
	virtual void evolveAllReplicas();

	void setSeed(const int s);
private:
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
	 * \brief Returns the sum of bitArray + 1
	 *
	 * \param bitArray the bit-array
	 * \param n the number of bits
	 */
	QVector<int> oneBitSum(const QVector<int> bitArray, int n);
	/**
	 * \brief Computes all the bit arrays of length n
	 *
	 * \param n the number of bits
	 */
	QVector<QVector<int> > allBitArrays(int n);
	/**
	 * \brief Initialises the population
	 */
	void initPop();
	/**
	 * \brief Resets the population (new size = 0)
	 */
	void resetPop();
	/**
	 * \brief Computes the circular shift of the input bit-array
	 *
	 * \param bitArray the bit-array
	 * \return the bit-array after the circular shift
	 */
	static QVector<int> circularShift(QVector<int> bitArray);
	/**
	 * \brief Computes the bit sum of the input bit-array
	 *        It is worth noting that the first half of
	 *        the bit-array represents the first number,
	 *        whereas the second half is the second number.
	 *
	 * \param bitArray the bit-array
	 * \return the bit-array after the bit sum
	 */
	static QVector<int> bitSum(QVector<int> bitArray);
	/**
	 * \brief Computes the bit multiplication of the input bit-array
	 *        It is worth noting that the first half of
	 *        the bit-array represents the first number,
	 *        whereas the second half is the second number.
	 *
	 * \param bitArray the bit-array
	 * \return the bit-array after the bit multiplication
	 */
	static QVector<int> bitMultiplier(QVector<int> bitArray);
	/**
	 * \brief Computes the even parity of the input bit-array
	 *
	 * \param bitArray the bit-array
	 * \return the even parity of the bit-array
	 */
	static QVector<int> evenParity(QVector<int> bitArray);
	/**
	 * \brief Computes the odd parity of the input bit-array
	 *
	 * \param bitArray the bit-array
	 * \return the odd parity of the bit-array
	 */
	static QVector<int> oddParity(QVector<int> bitArray);
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
	 * \brief Saves real fitness statistics (min, max, avg)
	 *
	 * \param logicOpCount how many (average) logic operators belong to the functional circuit
	 *        (it takes into consideration all times an operator is used)
	 * \param singleLogicOpCount how many (average) logic operators belong to the functional circuit
	 *        (each operator is considered once regardless of how many times it appears)
	 */
	void saveLogicOperatorCounter(float logicOpCount, float singleLogicOpCount);
	/**
	 * \brief Fills a QVector with -1 values (invalid)
	 *
	 * \param ind The QVector
	 */
	void invalidFill(QVector<int>& ind);
	/**
	 * \brief Saves all information concerning logic operators
	 *
	 */
	void saveFullLogicOpInfo();

	farsa::RandomGenerator* getRng();

	void saveEvalStats(const int steps, const double fit);

	int** generateBitstrings(int n);

	//! Code for printing a binary tree
	int maxHeight(BinaryTree *p);

    std::string intToString(int val);

	void printBranches(int branchLen, int nodeSpaceLen, int startLen, int nodesInThisLevel, const std::deque<BinaryTree*>& nodesQueue, std::ostream& out);

    void printNodes(int branchLen, int nodeSpaceLen, int startLen, int nodesInThisLevel, const std::deque<BinaryTree*>& nodesQueue, std::ostream& out);

    void printLeaves(int indentSpace, int level, int nodesInThisLevel, const std::deque<BinaryTree*>& nodesQueue, std::ostream& out);

    void printPretty(BinaryTree *root, int level, int indentSpace, std::ostream& out);

	void treeLogicOperator(int node, QVector<int> genotype, BinaryTree* root);

	void pathLogicOperator(int node, QVector<int> genotype, QVector<int>& tree, int type);

	void typeLogicOperator(int node, QVector<int> genotype, QVector<int>& boolean_function);

	void countLogicOperator(QVector<int> tree, QVector<int>& singleOpTree);
	//! End code

	//! BOOLEAN OPERATORS
	
	/**
	 * \brief Boolean and
	 */
	int booleanAnd(int, int);
	/**
	 * \brief Boolean or
	 */
    int booleanOr(int, int);
	/**
	 * \brief Boolean xor
	 */
    int booleanXor(int, int);
	/**
	 * \brief Boolean nand
	 */
    int booleanNand(int, int);
	/**
	 * \brief Boolean nor
	 */
    int booleanNor(int, int);

	//! Number of evaluation steps
	int m_nevalsteps;
	//! The input bit-string
	QVector<int> m_input;
	//! The number of bits of the input
	int m_inputSize;
	//! The number of boolean operators in the circuit
	int m_numBooleanOperators;
	//! The number of layers
	int m_numLayers;
	//! Population
	QVector<QVector<int> > m_pop;
	//! How many operators have been defined (= # MACROS)
	int m_boolOps;
	//! The input bit-string
	QVector<int> m_output;
	//! The number of bits of the output
	int m_outputSize;
	//! QVector for storing the values of the circuit nodes
	QVector<int> m_nodes;
	//! Mutation rate
	double m_mutationRate;
	//! Flags whether or not the post-evaluation test must be performed
	bool m_test;

	//! Buffers for saving best genotypes and best fitness
	QString m_bestGenBuf;
	QString m_statBuf;
	QString m_logicOpBuf;
	QString m_fullLogicOpBuf;
	QString m_evalstatBuf;

	//! Number of generations after which the buffers should be reset
	int m_resetBufNumGen;
	//! Flags whether or not the information about the best individual must be saved
	bool m_saveBestGenInfo;
	//! Interval for saving information about the best individual
	int m_bestGenInterval;
	//! Noise on the individual's fitness (i.e., stochasticity)
	double m_fitNoise;
	//! Logic operator array
	QVector<double> m_logicOpArray;
	//! Logic operator array size (i.e., information to store and save)
	int m_logicOpInfo;
	//! Random number generator
	farsa::RandomGenerator m_rng;

	//! Macros for boolean operators
	#define AND(a,b) booleanAnd(a,b);
	#define OR(a,b) booleanOr(a,b);
	#define XOR(a,b) booleanXor(a,b);
	#define NAND(a,b) booleanNand(a,b);
	#define NOR(a,b) booleanNor(a,b);
};

#endif
