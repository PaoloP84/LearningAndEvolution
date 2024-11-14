/********************************************************************************
 *  FARSA - Total99                                                             *
 *  Copyright (C) 2005-2011 Gianluca Massera <emmegian@yahoo.it>                *
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

#include "fastevolvabilityexperiment.h"
#include "utilitiesexceptions.h"
#include "randomgenerator.h"
#include "robots.h"
#include "world.h"
#include "phybox.h"
#include "wheeledexperimenthelper.h"
#include "logger.h"
#include "configurationhelper.h"

// This is needed because the isnan and isinf functions are not present windows
#ifdef WIN32
	#include <float.h>
	#define isnan(x) _isnan(x)
	#define isinf(x) (!_finite(x))
#else
	#define isnan(x) std::isnan(x)
	#define isinf(x) std::isinf(x)
#endif

FastEvolvabilityExperiment::FastEvolvabilityExperiment()
	: farsa::EvoRobotExperiment()
{
	// No need to add resources, they are added by EvoRobotExperiment
}

FastEvolvabilityExperiment::~FastEvolvabilityExperiment()
{
}

void FastEvolvabilityExperiment::configure(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::EvoRobotExperiment::configure(params, prefix);
}

void FastEvolvabilityExperiment::save(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::EvoRobotExperiment::save(params, prefix);

    farsa::Logger::error("NOT IMPLEMENTED (FastEvolvabilityExperiment::save)");
	abort();
}

void FastEvolvabilityExperiment::describe(QString type)
{
	// Calling parent function
	farsa::EvoRobotExperiment::describe(type);
	Descriptor d = addTypeDescription(type, "The experiment in which the electronic circuit should compute the circular shift of the input bits");
}

void FastEvolvabilityExperiment::postConfigureInitialization()
{
	// Calling parent function
	farsa::EvoRobotExperiment::postConfigureInitialization();
}

void FastEvolvabilityExperiment::initGeneration(int /*generation*/)
{
}

void FastEvolvabilityExperiment::initIndividual(int /*individual*/)
{
}

void FastEvolvabilityExperiment::initTrial(int /*trial*/)
{
	// Resetting fitness for the current trial
	trialFitnessValue = 0;
}

void FastEvolvabilityExperiment::initStep(int /*step*/)
{
}

void FastEvolvabilityExperiment::afterSensorsUpdate()
{
}

void FastEvolvabilityExperiment::beforeMotorsUpdate()
{
}

void FastEvolvabilityExperiment::beforeWorldAdvance()
{
}

void FastEvolvabilityExperiment::endStep(int /*step*/)
{
	farsa::ResourcesLocker locker(this);
	trialFitnessValue += 1.0;
}

void FastEvolvabilityExperiment::endTrial(int /*trial*/)
{
	totalFitnessValue += trialFitnessValue / farsa::real(getNSteps());
}

void FastEvolvabilityExperiment::endIndividual(int /*individual*/)
{
	totalFitnessValue = totalFitnessValue / farsa::real(getNTrials());
}

void FastEvolvabilityExperiment::endGeneration(int /*generation*/)
{
}
