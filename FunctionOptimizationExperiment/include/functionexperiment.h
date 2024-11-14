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

#ifndef FUNCTIONEXPERIMENT_H
#define FUNCTIONEXPERIMENT_H

#include "farsaplugin.h"
#include "evorobotexperiment.h"
#include "wvector.h"
#include <QVector>

/**
 * \brief An experiment in which a khepera robot has to discriminate between an
 *        object and the arena walls
 *
 * The resources used by this experiment are the same as the ones of the parent
 * class (EvoRobotExperiment)
 *
 * The parameters for this experiment are:
 * 	- distanceThreshold: the distance from the object below which the robot
 * 	                     is rewarded. This is the distance of the nearest
 * 	                     points of the robot and the object
 * 	- playgroundWidth: the width of the part of the arena surrounded by
 * 	                   walls (the playground is the area where the robot can
 * 	                   move)
 * 	- playgroundHeight: the height of the part of the arena surrounded by
 * 	                    walls (the playground is the area where the robot
 * 	                    can move)
 * 	- minObjectDistanceFromWall: the minimum distance from the walls at
 * 	                             which the object is placed. The default is
 * 	                             0.05
 * 	- minInitialRobotDistanceFromObject: the minimum distance from the
 * 	                                     object at which the robot is
 * 	                                     initialized. The default is 0.1
 * 	- minInitialRobotDistanceFromWall: the minimum distance from the walls
 * 	                                   at which the robot is initialized.
 * 	                                   The default is 0.1
 */
class FARSA_PLUGIN_API FunctionExperiment : public farsa::EvoRobotExperiment
{
	Q_OBJECT
	FARSA_REGISTER_CLASS(EvoRobotExperiment)

public:
	/**
	 * \brief Constructor
	 */
    FunctionExperiment();

	/**
	 * \brief Destructor
	 */
    ~FunctionExperiment();

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

	/**
	 * \brief Saves the actual status of parameters into the
	 *        ConfigurationParameters object passed
	 *
	 * \param params the configuration parameters object on which to save
	 *               the actual parameters
	 * \param prefix the prefix to use to access the object configuration
	 *               parameters
	 */
	virtual void save(farsa::ConfigurationParameters& params, QString prefix);

	/**
	 * \brief Add to Factory::typeDescriptions() the descriptions of all
	 *        parameters and subgroups
	 *
	 * It's mandatory in all subclasses where configure and save methods
	 * have been re-implemented for dealing with new parameters and
	 * subgroups to also implement the describe method
	 * \param type is the name of the type regarding the description. The
	 *             type is used when a subclass reuse the description of its
	 *             parent calling the parent describe method passing the
	 *             type of the subclass. In this way, the result of the
	 *             method describe of the parent will be the addition of the
	 *             description of the parameters of the parent class into
	 *             the type of the subclass
	 */
	static void describe( QString type );

	/**
	 * \brief This function is called after all linked objects have been
	 *        configured
	 *
	 * See the description of the ConfigurationParameters class for more
	 * information. This method creates the arena and fills it with objects
	 */
	virtual void postConfigureInitialization();

	/**
	 * \brief Called at the begin of each generation
	 *
	 * \param generation the generation about to start
	 */
	virtual void initGeneration(int generation);

	/**
	 * \brief Called before the evaluation of a new individual
	 *
	 * \param individual the id of the individual about to be tested
	 */
	virtual void initIndividual(int individual);

	/**
	 * \brief Called at the beginning of each trial
	 *
	 * \param trial the trial about to start
	 */
	virtual void initTrial(int trial);

	/**
	 * \brief Called at the beginning of each step (before world advance)
	 *
	 * \param step the step about to start
	 */
	virtual void initStep(int step);

	/**
	 * \brief Celled after all sensors have been updated but before network
	 *        spreading
	 *
	 * This is useful, for example, to overwrite the inputs of the neural
	 * network (i.e.: to silence some neurons during the experiment without
	 * modifing sensors classes)
	 */
	virtual void afterSensorsUpdate();

	/**
	 * \brief Called just before updating motors and after updating the
	 *        neural network
	 *
	 * This is useful, for example, to overwrite the outputs of the neural
	 * network
	 */
	virtual void beforeMotorsUpdate();

	/**
	 * \brief Called just before the world advances, after the update of
	 *        motors
	 *
	 * This is useful, for example, to manually actuate motors overriding
	 * the robot controller commands
	 */
	virtual void beforeWorldAdvance();

	/**
	 * \brief Called at the end of each step
	 *
	 * \param step the step about to end
	 */
	virtual void endStep(int step);

	/**
	 * \brief Called at the end of each trial
	 *
	 * \param trial the trial about to end
	 */
	virtual void endTrial(int trial);

	/**
	 * \brief Called at the end of an individual life
	 *
	 * \param individual the individual that has been just tested
	 */
	virtual void endIndividual(int individual);

	/**
	 * \brief Called at the end of each generation
	 *
	 * \param generation the generation about to end
	 */
	virtual void endGeneration(int generation);
};

#endif
