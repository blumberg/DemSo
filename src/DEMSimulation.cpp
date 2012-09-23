/*
 *   DemSo - 2D Discrete Element Method for Soil application
 *   Copyright (C) 2012  UNICAMP FEM/DMC
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* Functions for the DEMSimulation class defined in datatypes.hpp */

#include <iostream>
#include <string>
#include <cstdlib>
//#include "cutil_math.h"
#include "datatypes.hpp"
#include "parser.hpp"

using namespace std;

// FIXME: move this to a DEMParticleType.cpp
void DEMProperties::addParticleType (DEMParticleType newType)
{
	particleTypes.push_back(newType);
}

int DEMProperties::particleTypeIndexById (string id)
{
	for (register int i = 0;i < particleTypes.size();i++)
		if (particleTypes[i].id == id) return i;
	throw "Nonexistant particle id: " + id;
}
// --- end FIXME

void DEMSimulation::loadFromFile (const char *filename)
{
	DEMParser parser(filename);

	try {
		parameters = parser.loadParameters();
		environment = parser.loadEnvironment();
		properties = parser.loadProperties();
		particles = parser.loadParticles(&properties);
	}
	catch (string errorMsg) {
		cerr << "Parsing Error: " << errorMsg << endl;
		exit (1);
	}
}

void DEMSimulation::printConfiguration (void)
{
	cout << "-- Parameters" << endl;
	cout << "timeStep: " << parameters.timeStep << endl;
	cout << endl;
	cout << "-- Environment" << endl;
	cout << "dimensions: (" << environment.dimension.x;
	cout << ", " << environment.dimension.y;
	cout << ", " << environment.dimension.z << ")" << endl;
	cout << "gravity: (" << environment.gravity.x;
	cout << ", " << environment.gravity.y;
	cout << ", " << environment.gravity.z << ")" << endl;
	cout << "boundaryNormalStiffness: " << environment.boundaryNormalStiffness << endl;
	cout << "boundaryShearStiffness: " << environment.boundaryShearStiffness << endl;
	cout << "boundaryDamping: " << environment.boundaryDamping << endl;
	cout << endl;
	cout << "-- Properties" << endl;
	cout << "numParticleTypes: " << properties.particleTypes.size() << endl;
	cout << endl;
	for (register int i = 0;i < properties.particleTypes.size();i++)
	{
		cout << "\t---- Particle Type " << i+1 << endl;
		cout << "\tid: " << properties.particleTypes[i].id << endl;
		cout << "\tname: " << properties.particleTypes[i].name << endl;
		cout << "\tcolor: (" << properties.particleTypes[i].color.x;
		cout << ", " << properties.particleTypes[i].color.y;
		cout << ", " << properties.particleTypes[i].color.z << ")" << endl;
		cout << "\tmass: " << properties.particleTypes[i].mass << endl;
		cout << "\tradius: " << properties.particleTypes[i].radius << endl;
		cout << "\tnormalStiffness: " << properties.particleTypes[i].normalStiffness << endl;
		cout << "\tshearStiffness: " << properties.particleTypes[i].shearStiffness << endl;
		cout << "\tnormalDamping: " << properties.particleTypes[i].normalDamping << endl;
		cout << endl;
	}
}
