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
	throw "Nonexistant particle id: '" + id + "'";
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
	cout << ", " << environment.dimension.y << ")" << endl;
	cout << "gravity: (" << environment.gravity.x;
	cout << ", " << environment.gravity.y << ")" << endl;
	cout << "boundaryNormalStiffness: " << environment.boundaryNormalStiffness << endl;
	cout << "boundaryShearStiffness: " << environment.boundaryShearStiffness << endl;
	cout << endl;
	cout << "-- Properties" << endl;
	cout << "numParticleTypes: " << properties.particleTypes.size() << endl;
	cout << endl;
	for (register int i = 0;i < properties.particleTypes.size();i++)
	{
		cout << "\t---- Particle Type " << i << endl;
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
		cout << "\tboundaryDamping: " << properties.particleTypes[i].boundaryDamping << endl;
		cout << "\tfrictionCoefficient: " << properties.particleTypes[i].frictionCoefficient << endl;
		cout << endl;
	}
	cout << "-- Particle Rectangles" << endl;
	cout << "numRectangles: " << particles.start.size () << endl;
	for (register int i = 0; i < particles.start.size(); i++)
	{
		cout << "\t---- Particle Rectangle " << i << endl;
		cout << "\tstart: (" << particles.start[i].x << ", " << particles.start[i].y << ")" << endl;
		cout << "\tend: (" << particles.end[i].x << ", " << particles.end[i].y << ")" << endl;
		cout << "\tnum: (" << particles.num[i].x << ", " << particles.num[i].y << ")" << endl;
		cout << "\ttypes:";
		for (register int j = 0; j < particles.types[i].size(); j++)
			cout << " [" << particles.types[i][j] << "] = "
				 << properties.particleTypes[particles.types[i][j]].id;
		cout << endl;
	}
	cout << endl;
	cout << "-- Particle Triangles" << endl;
	cout << "numTriangles: " << particles.t_pos.size() << endl;
	for (register int i = 0; i < particles.t_pos.size(); i++)
	{
		cout << "\t---- Particle Triangle " << i << endl;
		cout << "\tpos: (" << particles.t_pos[i].x << ", " << particles.t_pos[i].y << ")" << endl;
		cout << "\twidth: " << particles.width[i] << endl;
		cout << "\tnum: " << particles.t_num[i] << endl;
		cout << "\ttypes:";
		for (register int j = 0; j < particles.t_types[i].size(); j++)
			cout << " [" << particles.t_types[i][j] << "] = "
				 << properties.particleTypes[particles.t_types[i][j]].id;
		cout << endl;
	}
	cout << endl;
	cout << "-- Single Particles" << endl;
	cout << "numSingleParticles: " << particles.pos.size () << endl;
	for (register int i = 0; i < particles.pos.size(); i++)
	{
		cout << "\t---- Particle " << i << endl;
		cout << "\ttype: [" << particles.type[i] << "] = "
			 << properties.particleTypes[particles.type[i]].id << endl;
		cout << "\tpos: (" << particles.pos[i].x << ", " << particles.pos[i].y << ")" << endl;
		cout << "\tvel: (" << particles.vel[i].x << ", " << particles.vel[i].y << ")" << endl;
		cout << "\ttheta: " << particles.theta[i] << endl;
		cout << "\tomega: " << particles.omega[i] << endl;
		cout << endl;
	}
}
