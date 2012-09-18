/* Functions for the DEMParticles class defined in datatypes.hpp */

#include <iostream>
#include <string>
#include <cstdlib>
#include "cutil_math.h"
#include "datatypes.hpp"
#include "parser.hpp"

using namespace std;

void DEMParticles::addParticles (DEMParticles newparts)
{
	positions.insert	(positions.end(),		newparts.positions.begin(),		newparts.positions.end());
	velocities.insert	(velocities.end(),		newparts.velocities.begin(),	newparts.velocities.end());
	accelerations.insert(accelerations.end(),	newparts.accelerations.begin(), newparts.accelerations.end());
	typeIndexes.insert	(typeIndexes.end(), 	newparts.typeIndexes.begin(),	newparts.typeIndexes.end());
}

void DEMParticles::generateBlock (int typeIndex, float3 start, float3 end, float3 spacing, DEMProperties *properties)
{
	float3 cubeSize = end - start;
	float radius = properties->particleTypes[typeIndex].radius;

	long int numParticles = (long int) floor((cubeSize.x + spacing.x)/(2*radius + spacing.x))
									  *floor((cubeSize.y + spacing.y)/(2*radius + spacing.y))
/*									  *floor((cubeSize.z + spacing.z)/(2*radius + spacing.z))*/; 
										// Descomentar a linha acima para o caso 3D

	cout << "Bloco tipo: " << properties->particleTypes[typeIndex].id << endl;
	cout << "numParticles: " << numParticles << endl;	//DEBUG
	cout << "Cubesize: (" << cubeSize.x << ", " << cubeSize.y << ", " << cubeSize.z << ")" << endl;
	cout << "Spacing: (" << spacing.x << ", " << spacing.y << ", " << spacing.z << ")" << endl;

	positions.reserve(numParticles);
	velocities.reserve(numParticles);
	accelerations.reserve(numParticles);
	typeIndexes.reserve(numParticles);

	register int counter = 0;
	for(float y = start.y; y + 2*radius < end.y; y += 2*radius + spacing.y)
	{
		for(float x = start.x; x + 2*radius < end.x; x += 2*radius + spacing.x)
		{
			positions.push_back(make_float3(x+radius, y+radius, 0));// TODO: Caso 3D
			velocities.push_back(make_float3(0.0));					// TODO: Deveria ser possível dar vel. iniciais
			accelerations.push_back(make_float3(0.0));				// TODO: idem
			typeIndexes.push_back(typeIndex);
			counter++;
		}
	}
	cout << "counter: " << counter << endl;	// DEBUG
	cout << "positions size: " << positions.size() << endl;
}
