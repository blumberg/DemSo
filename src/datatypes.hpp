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
#include "cutil_math.h"
#include "main.cuh"

using std::vector;

class DEMParameters {
	public:
		float timeStep;
		float simDuration;
		vector<int> followedParticles;
		uint imageDIMy;
};

class DEMEnvironment {
	public:
		float2 dimension;
		float2 gravity;
		float boundaryNormalStiffness;
		float boundaryShearStiffness;
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
		float boundaryDamping;
		float frictionCoefficient;
		float attractCoefficient;
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
		// Para retângulos de partículas
		vector<float2>	start;	// Coordenadas do canto inferior-esquerdo
		vector<float2>	end;	// Coordenadas do canto superior-direito
		vector<uint2>	num;	// Número de partículas desejado em x, y
		vector< vector<int> > types;	// Tipos de partículas à sortear dentro do retângulo

		// Para triângulos de partículas
		vector<float2>	t_pos;	// Posição do meio da aresta inferior
		vector<float>	width;	// Largura da base do triângulo
		vector<uint>	t_num;  // Número de partículas na base do triângulo
		vector< vector<int> > t_types;	// Tipos de partículas à sortear dentro do triângulo

		// Para partículas avulsas
		vector<float2>	pos;	// Posições x, y
		vector<float2>	vel;	// Velocidades x, y
		vector<float>	theta;	// Posições angulares
		vector<float>	omega;	// Velocidades angulares
		vector<int>		type;	// Tipos das partículas
		vector<int>		followedParticles; // Partículas a serem seguidas

#if USE_BIG_PARTICLE
		int controlledType;	// Tipo da particula controlada
#endif
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
