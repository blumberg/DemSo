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

/* vim: set filetype=cpp */

#ifndef MAIN_CUH
#define MAIN_CUH

#include <stdio.h>
#include <vector>

#define USE_TEX 0
#define USE_BIG_PARTICLE 1

#define FPS 31.0f
#define MAX_PARTICLES_TYPES 10


// Estrutura de propriedades das partículas, caso houver mais do que um
// tipo de partícula, essa estrutura será criada com um tamanho maior
// Essa estrutura será carregana na memória de constantes da GPU
struct ParticleProperties {

	float 	radius;
	float 	mass;
	float	inertia;
	float 	normalStiffness;
	float	shearStiffness;
	float 	normalDamping;
	float	boundaryDamping;
	float	frictionCoefficient;
	float 	colorR;
	float 	colorG;
	float 	colorB;

};

// Estrutura com os valores armazenados de cada partícula. Todas as
// variáveis dessa estrutura serão alocadas na GPU.
struct ParticlesValues
{
	float *pos1, *vel1, *acc;
	float *pos2, *vel2;

	float *theta1, *omega1, *alpha;
	float *theta2, *omega2;

	uint *ID1, *type1, *loc1;
	uint *ID2, *type2, *loc2;

	uint *cellStart, *cellEnd;
	uint *gridParticleIndex, *gridParticleHash;
	uint *fixParticleIndex, *ctrlParticleIndex;
	
	float2 controlPos;
	uint controlType;

};

// Estrutura com as propriesdades do sistema. Seu tamanho será fixo e esta
// estrutura inteira será passada para a memória de constantes da GPU
struct SystemProperties {

    uint numParticles;

    float2 cubeDimension;
    uint2 gridSize;
	uint numCells;
    
    float timeStep;
    
    float2 gravity;

	float boundaryNormalStiffness;
	float boundaryShearStiffness;
};

struct RenderParameters {

	int imageDIMx;
    int imageDIMy;
	int dimx;
	int dimy;
	float pRadius;

};

struct TimeControl {

	clock_t start, totalStart;
	int tempo;
	int IPS;
};

// Estrutura principal da simulação. Ela contem as 3 outras subestruturas.
struct DataBlock {

	ParticlesValues partValues;
	SystemProperties sisProps;
	RenderParameters renderPar;
	TimeControl timeCtrl;
	FILE * outputFile;
	std::vector<int> followedParticles;
};

#endif /* MAIN_CUH */
