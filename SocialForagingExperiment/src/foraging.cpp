#include "foraging.h"
#include "configurationhelper.h"
#include "randomgenerator.h"
#include "helperresources.h"

// This is needed because the isnan and isinf functions are not present windows
#ifdef WIN32
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#else
#define isnan(x) std::isnan(x)
#define isinf(x) std::isinf(x)
#endif

#define LOCAL_SEED 1234

namespace farsa
{

BasicForagingExperiment::BasicForagingExperiment() :
	EvoRobotExperiment(),
	m_robotPos(),
	m_foods(),
	m_foodsNumber(10),
	m_recreateWorld(false),
	m_minDistance(0.25),
	m_maxNumAttempts(1000),
	m_collidedObjectInfo(),
	m_rng(1)
{
}


BasicForagingExperiment::~BasicForagingExperiment()
{
}

void BasicForagingExperiment::configure(ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	EvoRobotExperiment::configure(params, prefix);

	// Loading our parameters. Checks on m_playgroundWidth and m_playgroundHeight will be done in setupArena() because
	m_foodsNumber = ConfigurationHelper::getInt(params, prefix + "foodsNumber", m_foodsNumber);
	m_minDistance = ConfigurationHelper::getDouble(params, prefix + "minDistance", m_minDistance);
}

void BasicForagingExperiment::save(ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	EvoRobotExperiment::save(params, prefix);

	Logger::error("NOT IMPLEMENTED (BasicForagingExperiment::save)");
	abort();
}

void BasicForagingExperiment::describe(QString type)
{
	// Calling parent function
	EvoRobotExperiment::describe(type);
	Descriptor d = addTypeDescription(type, "The experiment in which a khepera robot has to eat foods");

	d.describeInt("foodsNumber").def(15).limits(1, +MaxInteger).help("The number of foods in the arena","This is the number of foods that are in the arena and that can give fitness to the robot");
	d.describeReal("minDistance").def(0.25).limits(0.0, 100.0).help("The minimum distance between objects and between robot and objects");
}

void BasicForagingExperiment::postConfigureInitialization()
{
	// Calling parent function
	EvoRobotExperiment::postConfigureInitialization();

	farsa::ResourcesLocker locker(this);
	farsa::Arena* arena = getResource<farsa::Arena>("arena");
	farsa::Evonet* evonet = getResource<farsa::Evonet>("agent[0]:evonet");
	m_collidedObjectInfo.resize(m_robots.size());
	// Here we create the arena and the object
    	setupArena();
}

void BasicForagingExperiment::initGeneration(int /*generation*/)
{
}

void BasicForagingExperiment::initIndividual(int /*individual*/)
{
	// Resetting the world to avoid numerical problems
	m_recreateWorld = true;

	setSeed(LOCAL_SEED);
}

void BasicForagingExperiment::initTrial(int trial)
{
    	if (m_recreateWorld) {
        	recreateWorld();
        	recreateAllRobots();
        	recreateArena();

        	setupArena();

        	m_recreateWorld = false;
    	}

    	ResourcesLocker locker(this);
    	Arena* arena = getResource<Arena>("arena");

	// Now placing the foods
    	if (m_foodsNumber > 0) {
		// The first food item is placed here
        	const real arenaWidthHalfLimitForObject = arena->getWidth() / 2.0 - m_foods[0]->phyObject()->radius() - 0.1;
        	const real arenaHeightHalfLimitForObject = arena->getHeight() / 2.0 - m_foods[0]->phyObject()->radius() - 0.1;

		real ox = globalRNG->getDouble(-arenaWidthHalfLimitForObject, arenaWidthHalfLimitForObject);
		real oy = globalRNG->getDouble(-arenaHeightHalfLimitForObject, arenaHeightHalfLimitForObject);
		m_foods[0]->setPosition(ox, oy);
		// Place the remaining food items
        	for (int i = 1; i < m_foodsNumber; i++) 
		{
			placeFood(i);	
   		}
	}

    	// Now placing the robots
	for (int i = 0; i < m_robots.size(); i++)
	{
    		RobotOnPlane* robot = getResource<RobotOnPlane>(m_robotNames[i]);
	    	//Epuck* r=dynamic_cast<Epuck*>(robot);
		farsa::PhyEpuck* r = dynamic_cast<farsa::PhyEpuck*>(robot);
		r->setColor(QColor(0,0,255));

	   	unsigned int numAttempts = 0;
	    	bool goOn = true;
	    	do 
		{
			// We have 100 attempts to place the robot
			if (numAttempts > m_maxNumAttempts) 
			{
		    		throwUserRuntimeError("Unable to place the robot (more than " + QString::number(m_maxNumAttempts) + " attempts made)");
			}
			numAttempts++;
			// Robot is placed here!!!
			const real arenaWidthHalfLimitForRobot = arena->getWidth() / 2.0 - robot->robotRadius() - 0.1;
			const real arenaHeightHalfLimitForRobot = arena->getHeight() / 2.0 - robot->robotRadius() - 0.1;
			const real rx = globalRNG->getDouble(-arenaWidthHalfLimitForRobot, arenaWidthHalfLimitForRobot);
			const real ry = globalRNG->getDouble(-arenaHeightHalfLimitForRobot, arenaHeightHalfLimitForRobot);
			robot->setPosition(arena->getPlane(), rx, ry);
			goOn = false;
			// Check whether robot position overlaps with foods
			if (m_foodsNumber > 0) 
			{
		    		for (int f = 0; f < m_foodsNumber; f++) 
				{
		        		double dist; 
					const wVector robotPosition(robot->position().x, robot->position().y, 0.0);
		        		const wVector objectPosition(m_foods[f]->position().x, m_foods[f]->position().y, 0.0);
		        		dist = (robotPosition - objectPosition).norm() - robot->robotRadius() - m_foods[f]->phyObject()->radius();

		        		if (dist < m_minDistance)
		            			goOn = true;
		    		}
			}
			if (!goOn)
			{
				// Check overlap with robots
				for (int j = 0; j < i; j++)
				{
					double dist; 
					const wVector robotPosition(robot->position().x, robot->position().y, 0.0);
		        		const wVector otherPosition(m_robots[j]->position().x, m_robots[j]->position().y, 0.0);
		        		dist = (robotPosition - otherPosition).norm() - robot->robotRadius() * 2.0;

		        		if (dist < m_minDistance)
		            			goOn = true;
				}
			}
	    	} while (goOn);					

	    	// Changing also the robot orentation (random orientation)
	    	robot->setOrientation(arena->getPlane(), getRng()->getDouble(-PI_GRECO, PI_GRECO));
	}

    	// Resetting fitness and event counters for the current trial
    	trialFitnessValue = 0.0;

	// Reset neural network
	const int nagents = getNAgents();
	for (int i = 0; i < nagents; i++)
	{
		Evonet* evonet = getNeuralNetwork(i);
		evonet->resetNet();
	}
}

void BasicForagingExperiment::initStep(int /*step*/)
{
}

void BasicForagingExperiment::afterSensorsUpdate()
{
}

void BasicForagingExperiment::beforeMotorsUpdate()
{
}

void BasicForagingExperiment::beforeWorldAdvance()
{
}

void BasicForagingExperiment::endStep(int step)
{
	double dist;
	double nAng;

    	ResourcesLocker locker(this);
    	Arena* arena = getResource<Arena>("arena");

	int rob = 0;
	foreach(farsa::RobotOnPlane* robot, m_robots)
	{
		// Check collisions with either food objects or walls
		QVector< farsa::PhyObject2DWrapper * > collidedObjects = arena->getKinematicRobotCollisions(m_robotNames[rob]);
		int n = collidedObjects.size();
		if (n != 0)
		{
			for (int i = 0; i < n; i++) 
			{
				if ( collidedObjects[i]->type() == farsa::PhyObject2DWrapper::Wall ) 
				{
					// Collision with a wall, we make the robot rebound
					// Set random orientation
					robot->setOrientation(arena->getPlane(), getRng()->getDouble(-PI_GRECO, PI_GRECO));
				}
				else if ( collidedObjects[i]->type() == farsa::PhyObject2DWrapper::BigCylinder ) 
				{
					// Collision with a food item
					// Look for the object (i.e., food item) the robot collided with
					farsa::wVector objPos = collidedObjects[i]->position();
					bool found = false;
					int f = -1;
					// Collision with one or more food objects
					for (int j = 0; (j < m_foodsNumber) && !found; j++)
					{
						if (objPos == m_foods[j]->position())
						{
							f = j;
							found = true;
						}
					}
					if (!found)
					{
						farsa::Logger::error("Unable to find the eaten object, stop execution!!!!");
						exit(-1);
					}
					// The robot collided with a different object.
					// We update its information
					m_collidedObjectInfo[rob] = f;
				}
			}
		}
		else
		{
			// Reset collided object info
			m_collidedObjectInfo[rob] = -1;
		}
		rob++;
	}
	bool foundPair = false;
	QVector<int> pairIdx(2);
	int i = 0;
	while (i < m_robots.size() - 1 && !foundPair)
	{
		if (m_collidedObjectInfo[i] != -1)
		{
			for (int j = i + 1; j < m_robots.size() && !foundPair; j++)
			{
				if (m_collidedObjectInfo[i] == m_collidedObjectInfo[j])
				{
					foundPair = true;
					pairIdx[0] = i;
					pairIdx[1] = j;
					break;
				}
			}
		}
		i++;
	}

	if (foundPair)
	{
		trialFitnessValue += 1.0;
		int objId = m_collidedObjectInfo[pairIdx[0]];
		// Place the eaten food object in a different position in the arena
		placeFood(objId);
		// We reset the object with which the robots collided
		m_collidedObjectInfo[pairIdx[0]] = -1;
		m_collidedObjectInfo[pairIdx[1]] = -1;
	}
}

void BasicForagingExperiment::endTrial(int /*trial*/)
{
    	// Add fitness for the trial to the total
    	totalFitnessValue += trialFitnessValue;
}

void BasicForagingExperiment::endIndividual(int /*individual*/)
{
	totalFitnessValue = totalFitnessValue / farsa::real(getNTrials());
}

void BasicForagingExperiment::endGeneration(int /*generation*/)
{
}

void BasicForagingExperiment::resourceChanged(QString resourceName, ResourceChangeType changeType)
{
        // Calling parent function
        EvoRobotExperiment::resourceChanged(resourceName, changeType);

        // Here we are only interested in robots, so we build a regular expression to only check when
        // a robot changes
        QRegExp checkRobots("agent\\[(\\d)\\]:robot");
        if (checkRobots.indexIn(resourceName) != -1) {
                // Getting the index
                int index = checkRobots.cap(1).toUInt();

                // The code below should work as expected. If there is a crash, then there is probably a bug
                // in resource management
                if (changeType == Deleted) {
                        // Removing robot
                        m_robots[index] = NULL;
                } else {
                        // Adding the robot to our list
                        farsa::RobotOnPlane* const robot = getResource<farsa::RobotOnPlane>();

                        if (m_robots.size() == index) {
                                m_robots.append(robot);
                                m_robotNames.append(QString("agent["+QString::number(index)+"]:robot"));
                        } else {
                                m_robots[index] = robot;
                                m_robotNames[index] = QString("agent["+QString::number(index)+"]:robot");
                        }
                }
        }
}


void BasicForagingExperiment::setupArena()
{
    	ResourcesLocker locker(this);

    	// Getting the arena
    	Arena* arena = getResource<Arena>("arena");

	//arena->getPlane()->setColor(Qt::black);

    	// Creating walls all around the arena
    	const real m_wallThickness = 0.01f;
    	const real m_wallHeight = 0.05f;
    	const real halfHeight = arena->getHeight() / 2.0;
    	const real halfWidth = arena->getWidth() / 2.0;
    	const real topWallPos = halfHeight + m_wallThickness / 2.0;
    	const real bottomWallPos = -(halfHeight + m_wallThickness / 2.0);
    	const real rightWallPos = halfWidth + m_wallThickness / 2.0;
    	const real leftWallPos = -(halfWidth + m_wallThickness / 2.0);
    	arena->createWall(Qt::green, wVector(-halfWidth, topWallPos, 0.0), wVector(halfWidth, topWallPos, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::green, wVector(-halfWidth, bottomWallPos, 0.0), wVector(halfWidth, bottomWallPos, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::green, wVector(rightWallPos, -halfHeight, 0.0), wVector(rightWallPos, halfHeight, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::green, wVector(leftWallPos, -halfHeight, 0.0), wVector(leftWallPos, halfHeight, 0.0), m_wallThickness, m_wallHeight);

    	// Now creating the foods
    	if (m_foodsNumber > 0)
	{
        	m_foods = new Cylinder2DWrapper*[m_foodsNumber];

        	for(int i = 0; i < m_foodsNumber; i++) 
		{
			m_foods[i] = arena->createBigCylinder(Qt::red, m_wallHeight);
			m_foods[i]->setStatic(true);
        	}
    	}

    	// Finally checking if the arena is big enough. We only check against the robot radius because it is bigger than the object radius
    	RobotOnPlane* robot = getResource<RobotOnPlane>("agent[0]:robot");
    	const real minDimension = robot->robotRadius() * 10;
    	if (arena->getWidth() < minDimension) 
	{
        	throwUserRuntimeError("Cannot run the experiment because the arena is too small (width too small)");
    	} 
	else if (arena->getHeight() < minDimension) 
	{
        	throwUserRuntimeError("Cannot run the experiment because the arena is too small (height too small)");
    	}
}

void BasicForagingExperiment::placeFood(int f)
{
	ResourcesLocker locker(this);

    	// Getting the arena
    	Arena* arena = getResource<Arena>("arena");
	RobotOnPlane* robot = getResource<RobotOnPlane>("agent[0]:robot");
	const real arenaWidthHalfLimitForObject = arena->getWidth() / 2.0 - m_foods[0]->phyObject()->radius() - 0.1;
	const real arenaHeightHalfLimitForObject = arena->getHeight() / 2.0 - m_foods[0]->phyObject()->radius() - 0.1;
	bool placed = false;
	while (!placed)
	{
		real ox = getRng()->getDouble(-arenaWidthHalfLimitForObject, arenaWidthHalfLimitForObject);
		real oy = getRng()->getDouble(-arenaHeightHalfLimitForObject, arenaHeightHalfLimitForObject);
		m_foods[f]->setPosition(ox, oy);
		placed = true;
		for (int i = 0; i < m_foodsNumber && placed; i++)
		{
			if(i != f)
			{
				double dist; 
				const wVector foodPosition(m_foods[f]->position().x, m_foods[f]->position().y, 0.0);
		        	const wVector otherPosition(m_foods[i]->position().x, m_foods[i]->position().y, 0.0);
				dist = (foodPosition - otherPosition).norm() - m_foods[f]->phyObject()->radius() - m_foods[i]->phyObject()->radius();
				if (dist < m_minDistance)
					placed = false;
			}
		}
		if (placed)
		{
			double dist; 
			const wVector foodPosition(m_foods[f]->position().x, m_foods[f]->position().y, 0.0);
			const wVector otherPosition(robot->position().x, robot->position().y, 0.0);
			dist = (foodPosition - otherPosition).norm() - m_foods[f]->phyObject()->radius() - robot->robotRadius();
			if (dist < m_minDistance)
				placed = false;
		}
	}
}

farsa::RandomGenerator* BasicForagingExperiment::getRng()
{
	return &m_rng;
}

void BasicForagingExperiment::setSeed(int s)
{
	m_rng.setSeed(s);
}

}
