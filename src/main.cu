// Includes que eu acho que preciso
#include <math.h>
#include <stdio.h>
#include <iostream>

// CUDA includes
#include <cuda.h>
#include <curand.h>
#include "../includes/cuda_by_example.h"
#include "../includes/gpu_anim.h"                   // este inclui o cuda.h
#include "../includes/cutil_math.h"

// THRUST includes
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include "thrust/device_ptr.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"

// Dependece files
#include "particles_kernel.cu"
#include "functions.cu"

#define DIM 800
#define PARTICLES 90000
#define BOX_SIZE 10.0f
#define TIME_STEP 1.0e-3
#define GRAVITY 9.81f
#define BOUNDARYDAMPING -0.5f
#define X_PARTICLES 300
#define Y_PARTICLES 300
#define ALWAYS 1

#define log2( x ) log(x)/log(2)

void PrepareSim( SistemProperties *params, ParticlesValues *particle, ParticleProperties *partProps ){

	params->numParticles = PARTICLES;

	params->cubeDimension.x = params->cubeDimension.y = BOX_SIZE;
	
	params->timeStep = TIME_STEP;
	
	params->gravity = make_float2(0,-GRAVITY);
	
	params->imageDIMx = DIM;
	params->imageDIMy = DIM;
		
	partProps->radius = 20e-3f;
	partProps->mass = 1e-2;
	partProps->collideStiffness = 1e3;
	partProps->collideDamping = 0.1f;
	partProps->boundaryDamping = BOUNDARYDAMPING;
	
	params->dimx = ceil(params->imageDIMx/params->cubeDimension.x*partProps->radius)*2;
	if (params->dimx < 2) params->dimx = 2;
	params->dimy = ceil(params->imageDIMy/params->cubeDimension.y*partProps->radius)*2;
	if (params->dimy < 2) params->dimy = 2;
	params->pRadius = params->imageDIMy/params->cubeDimension.y*partProps->radius;

	// Bloco inicial de esferas
	float corner1[2] = {0.1, 0.1};
	float corner2[2] = {9.9, 9.9};
	float sideLenght[2];

	// Grid dimension
	uint grid = params->cubeDimension.x / (4.0f * partProps[0].radius);
	uint temp = log2(grid);
	uint gridUpdate = pow(2,temp);
	float cellSize = params->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;
	params->gridSize.x = pow(2,temp);
	
	grid = params->cubeDimension.y / (4 * partProps[0].radius);
	temp = log2(grid);
	gridUpdate = pow(2,temp);
	cellSize = params->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;	
	params->gridSize.y = pow(2,temp);

	params->numCells = params->gridSize.x * params->gridSize.y;

	// Posicionando as primeiras particulas
	sideLenght[0] = corner2[0] - corner1[0];
	sideLenght[1] = corner2[1] - corner1[1];
	
	uint side[2] = {X_PARTICLES, Y_PARTICLES};
	
	// alocando vetores na placa de video
	float *d_corner1, *d_sideLenght;
	uint *d_side;

	cudaMalloc((void**)&d_corner1, sizeof(float)*2);
	cudaMalloc((void**)&d_sideLenght, sizeof(float)*2);
	cudaMalloc((void**)&d_side, sizeof(uint)*2);
	cudaMalloc((void**)&particle->pos1, sizeof(float2) * params->numParticles);
	cudaMalloc((void**)&particle->pos2, sizeof(float2) * params->numParticles);
	cudaMalloc((void**)&particle->vel1, sizeof(float2) * params->numParticles);
	cudaMalloc((void**)&particle->vel2, sizeof(float2) * params->numParticles);
	cudaMalloc((void**)&particle->acc, sizeof(float2) * params->numParticles);
	cudaMalloc((void**)&particle->cellStart, sizeof(uint) * params->numCells);
	cudaMalloc((void**)&particle->cellEnd, sizeof(uint) * params->numCells);
	cudaMalloc((void**)&particle->particleIndex, sizeof(uint) * params->numParticles);
	cudaMalloc((void**)&particle->particleHash, sizeof(uint) * params->numParticles);
	cudaMemcpy(d_corner1, corner1, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sideLenght, sideLenght, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_side, side, sizeof(uint)*2, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(simPropD, params, sizeof(SistemProperties));
	cudaMemcpyToSymbol(partPropD, partProps , sizeof(ParticleProperties) * 1);

	initializeParticlePosition(particle->pos1,
							   particle->vel1,
							   particle->acc,
							   d_corner1,
							   d_sideLenght,
							   d_side,
							   side);

    cudaFree( d_corner1 );
    cudaFree( d_sideLenght );
    cudaFree( d_side );


	// Screen output	
	printf("Number of Spheres = %d\n",params->numParticles);
	printf("grid %d x %d\n",params->gridSize.x,params->gridSize.y);
}

void SimLooping( uchar4 *pixels, DataBlock *simBlock, int ticks ) {

    SistemProperties *sisProps = &simBlock->sisProps;
//    ParticleProperties *partProps = &simBlock->partProps;
    ParticlesValues *partValues = &simBlock->partValues;

	float2 *oldPos, *oldVel, *sortPos, *sortVel;

	if ((ticks % 2))
	{	
		oldPos = partValues->pos1;
		sortPos = partValues->pos2;
		oldVel = partValues->vel1;
		sortVel = partValues->vel2;
	} else {
		oldPos = partValues->pos2;
		sortPos = partValues->pos1;
		oldVel = partValues->vel2;
		sortVel = partValues->vel1;
	}

	// Define a celula de cada particula
	calcHash(oldPos,
			 partValues->particleIndex,
			 partValues->particleHash,
			 sisProps->numParticles);

	// Ordena o grid pela posicao das particulas
	sortParticles(partValues->particleHash,
				  partValues->particleIndex,
				  sisProps->numParticles);

	// Encontra as particulas de inicializacao e de finalizacao
	reorderDataAndFindCellStart(partValues->cellStart,
								partValues->cellEnd,
								sortPos,
								sortVel,
								partValues->particleHash,
								partValues->particleIndex,
								oldPos,
								oldVel,
								sisProps->numParticles,
								sisProps->numCells);

	// Detecta a colizao das particulas
	collide(sortPos,
			sortVel,
			partValues->acc,
			partValues->particleIndex,
			partValues->cellStart,
			partValues->cellEnd,
			sisProps->numParticles,
			sisProps->numCells);

	// Integracao no tempo (atualizacao das posicoes e velocidades)
	integrateSystem(sortPos,
		 	  		sortVel,
		 	  		partValues->acc,
		 	  		sisProps->numParticles);

#if !ALWAYS
	if (ticks % 10 == 1){
#endif
		// Saida grarica quando necessario
		plotParticles(pixels,
					  sortPos,
					  sisProps->numParticles,
					  sisProps->imageDIMx,
					  sisProps->imageDIMy);

	//printf("Fim %d\n\n",ticks);
#if !ALWAYS
	}
#endif
}

void FinalizingSim( DataBlock *simBlock) {

    // Limpe aqui o que tiver que ser limpo
    
    SistemProperties *d1 = &simBlock->sisProps;
    ParticleProperties *d3 = &simBlock->partProps;
    ParticlesValues *d2 = &simBlock->partValues;
    
    cudaFree( d2->pos1 );
    cudaFree( d2->pos2 );
    cudaFree( d2->vel1 );
    cudaFree( d2->vel2 );
    cudaFree( d2->acc );
    cudaFree( d2->cellStart );
    cudaFree( d2->cellEnd );
    cudaFree( d2->particleIndex );
    cudaFree( d2->particleHash );
    cudaFree( d1 );
    cudaFree( d2 );
    cudaFree( d3 );

}


int main() {

    DataBlock simBlock;
    
    SistemProperties *sisProps = &simBlock.sisProps;
    ParticleProperties *partProps = &simBlock.partProps;
    ParticlesValues *partValues = &simBlock.partValues;
    
    GPUAnimBitmap bitmap(DIM, DIM, &simBlock );

	// Utilizar ARGC e ARGV para pegar propriedades na linha de comando
	// ler esses comandos de um arquivo TXT externo
	// Criar uma rotina para fazer este tipo de leitura
	
	// Prepara a simulacao, define as condicoes iniciais do problema
	PrepareSim(sisProps, partValues, partProps);

    bitmap.anim_and_exit(
        (void (*)(uchar4*,void*,int))SimLooping, (void (*)(void*))FinalizingSim );

}
