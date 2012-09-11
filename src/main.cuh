/*
 *   DemSo - 2D Discrete Element Method for Soil application
 *   Copyright (C) 2012  UNICAMP FEM/DMC
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
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

// Estrutura de propriedades das partículas, caso houver mais do que um
// tipo de partícula, essa estrutura será criada com um tamanho maior
// Essa estrutura será carregana na memória de constantes da GPU
struct ParticleProperties {

	float radius;
	float mass;
	float collideStiffness;
	float collideDamping;
	float boundaryDamping;

};

// Estrutura com os valores armazenados de cada partícula. Todas as
// variáveis dessa estrutura serão alocadas na GPU.
struct ParticlesValues {

	float *pos1, *vel1, *acc;
	float *pos2, *vel2;
	uint *cellStart, *cellEnd;
	uint *gridParticleIndex, *gridParticleHash;

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

	ParticleProperties partProps;
	ParticlesValues partValues;
	SystemProperties sisProps;
	RenderParameters renderPar;
	TimeControl timeCtrl;

};
