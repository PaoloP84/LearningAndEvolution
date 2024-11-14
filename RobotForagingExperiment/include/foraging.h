#ifndef BASICFORAGING_H
#define BASICFORAGING_H

#include "farsaplugin.h"
#include "evorobotexperiment.h"

namespace farsa {

class FARSA_PLUGIN_API BasicForagingExperiment : public EvoRobotExperiment
{
	Q_OBJECT
	FARSA_REGISTER_CLASS(EvoRobotExperiment)

public:
	/**
	 * \brief Constructor
	 */
    BasicForagingExperiment();

	/**
	 * \brief Destructor
	 */
    ~BasicForagingExperiment();

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
    virtual void configure(ConfigurationParameters& params, QString prefix);
	
	/**
	 * \brief Saves the actual status of parameters into the
	 *        ConfigurationParameters object passed
	 *
	 * \param params the configuration parameters object on which to save
	 *               the actual parameters
	 * \param prefix the prefix to use to access the object configuration
	 *               parameters
	 */
    virtual void save(ConfigurationParameters& params, QString prefix);
	

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
     *        spreadingm_playgroundWidth
	 *
	 * This is useful, for example, to overwrite the inputs of the neural
	 * network (i.e.: to silence some neurons during the experiment without
	 * modifing sensors classes)
	 */
	virtual void afterSensorsUpdate();

	/**
     * \brief
 Called just before updating motors and after updating the
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
     * This is useful, for example, to manually actuate motors overriding    virtual void resourceChanged(QString name, ResourceChangeType changeType);

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
    virtual void endIndividual(int);

	/**
	 * \brief Called at the end of each generation
	 *
	 * \param generation the generation about to end
	 */
    virtual void endGeneration(int generation);

private:
	/**
         * \brief Called to notify changes in resources
         *
         * We use this function to build a list of the robots in the experiment.
         * \param resourceName the name of the resource that has changed.
         * \param chageType the type of change the resource has gone through
         *                  (whether it was created, modified or deleted)
         */
    virtual void resourceChanged(QString resourceName, ResourceChangeType changeType);

    //---------- ARENA SETUP AND OBJECTS ----------//

    /** \brief Creates the arena */
	void setupArena();

	//! Place food items after being eaten
	void placeFood(int f);

	//! Get random generator
	farsa::RandomGenerator* getRng();

	//! Set seed
	void setSeed(int s);

	QVector<farsa::RobotOnPlane*> m_robots;
	QVector<QString> m_robotNames;

    	/** \brief If true re-creates the world at the beginning of the next trial */
    	bool m_recreateWorld;

	/** \brief Robot position */
	wVector m_robotPos;

    	/** \brief Number of objects in the arena */
    	int m_foodsNumber;

    	/** \brief The foods placed inside the arena */
   	Cylinder2DWrapper** m_foods;

	//! Minimum distance between objects
	double m_minDistance;

	//! Maximum number of attempts to place elements in the environment
	const unsigned int m_maxNumAttempts;

	//! Random number generator
	farsa::RandomGenerator m_rng;
};

}

#endif

