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
	cout << "FPS: " << parameters.framesPerSecond << endl;
	cout << endl;
	cout << "-- Environment" << endl;
	cout << "dimensions: (" << environment.dimension.x;
	cout << ", " << environment.dimension.y;
	cout << ", " << environment.dimension.z << ")" << endl;
	cout << "gravity: (" << environment.gravity.x;
	cout << ", " << environment.gravity.y;
	cout << ", " << environment.gravity.z << ")" << endl;
	cout << endl;
	cout << "-- Properties" << endl;
	cout << "numParticleTypes: " << properties.particleTypes.size() << endl;
	cout << endl;
	for (register int i = 0;i < properties.particleTypes.size();i++)
	{
		cout << "\t---- Particle Type " << i+1 << endl;
		cout << "\tid: " << properties.particleTypes[i].id << endl;
		cout << "\tname: " << properties.particleTypes[i].name << endl;
		cout << "\tmass: " << properties.particleTypes[i].mass << endl;
		cout << "\tradius: " << properties.particleTypes[i].radius << endl;
		cout << "\tnormalStiffness: " << properties.particleTypes[i].normalStiffness << endl;
		cout << "\tnormalDamping: " << properties.particleTypes[i].normalDamping << endl;
		cout << "\tboundaryDamping: " << properties.particleTypes[i].boundaryDamping << endl;
		cout << endl;
	}
	cout << "numParticles: " << particles.positions.size() << endl;
}
