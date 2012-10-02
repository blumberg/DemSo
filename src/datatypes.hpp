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

#ifndef DATATYPES_H
#define DATATYPES_H
#include <vector>
#include <math.h>
#include "../includes/cutil_math.h"

using std::vector;

class DEMParameters {
	public:
		float timeStep;
		vector<int> followedParticles;
		uint imageDIMy;
};

class DEMEnvironment {
	public:
		float2 dimension;
		float2 gravity;
		float boundaryNormalStiffness;
		float boundaryShearStiffness;
		float boundaryDamping;
		float frictionCoefficient;
};

class DEMParticleType {
	public:
		std::string id;
		std::string name;
		float mass;
		float radius;
		float normalStiffness;
		float shearStiffness;
		float normalDamping;
		float3 color;
		DEMParticleType () {};
};

class DEMProperties {
	public:
		std::vector<DEMParticleType> particleTypes;
		void addParticleType (DEMParticleType);
		int particleTypeIndexById (std::string);
};

// Classe dos valores de inicialização das partículas
class DEMParticles {
	public:
		// Para blocos de partículas
		float2 start;	// Coordenadas do canto inferior-esquerdo
		float2 end;		// Coordenadas do canto superior-direito
		float2 num;		// Número de partículas desejado em x, y

		// Para partículas avulsas
		vector<float2>	pos;	// Posições x, y
		vector<float2>	vel;	// Velocidades x, y
		vector<float>	theta;	// Posições angulares
		vector<float>	omega;	// Velocidades angulares
		vector<int>		type;	// Se para uma partícula i, type[i] == 1, então
								// o tipo desta partícula será o que estiver na posição
								// 1 no vetor particleTypes[]
		//void addParticles (DEMParticles);
		//void generateBlock (int, float3, float3, float3, DEMProperties *);
};

class DEMSimulation {
	public:
		DEMParameters 	parameters;
		DEMEnvironment	environment;
		DEMProperties	properties;
		DEMParticles	particles;
		void loadFromFile (const char *filename);
		void printConfiguration ();
};

#endif /* DATATYPES_H */
