// Input and Output include
#include <iostream>								  // entra e saída de dados

// CUDA includes
#include <curand_kernel.h>  		  // bib. randomica para kernel em CUDA
#include "../includes/gpu_anim.h" 	  // bib. de vizualização em tempo real
#include "../includes/cutil_math.h"       // funções matemáticas de vetores

// THRUST includes
#include "thrust/device_ptr.h"  		   // thrust para utilizar ponteiro
#include "thrust/sort.h" 					   // thrust para ordenar vetor

// Dependece files
#include "particles_kernel.cu" 				   // funções executadas na GPU
#include "functions.cu" 	 // arquivo de funções de preparação para a GPU

#define DIM 800
#define PARTICLES 90000
#define BOX_SIZE 10.0f
#define TIME_STEP 1.0e-3
#define GRAVITY 9.81f
#define BOUNDARYDAMPING -0.5f
#define X_PARTICLES 300
#define Y_PARTICLES 300
#define ALWAYS 0
#define NUM_FRAMES_PER_PLOT 2

#define log2( x ) log(x)/log(2)

void PrepareSim( SistemProperties *sisProps, ParticlesValues *partValues, ParticleProperties *partProps ){

	sisProps->numParticles = PARTICLES;

	sisProps->cubeDimension.x = sisProps->cubeDimension.y = BOX_SIZE;
	
	sisProps->timeStep = TIME_STEP;
	
	sisProps->gravity = make_float2(0,-GRAVITY);
	
	sisProps->imageDIMx = DIM;
	sisProps->imageDIMy = DIM;
		
	partProps->radius = 15e-3f;
	partProps->mass = 1e-2;
	partProps->collideStiffness = 1e3;
	partProps->collideDamping = 0.2f;
	partProps->boundaryDamping = BOUNDARYDAMPING;
	
	sisProps->dimx = ceil(sisProps->imageDIMx/sisProps->cubeDimension.x*partProps->radius)*2;
	if (sisProps->dimx < 2) sisProps->dimx = 2;
	sisProps->dimy = ceil(sisProps->imageDIMy/sisProps->cubeDimension.y*partProps->radius)*2;
	if (sisProps->dimy < 2) sisProps->dimy = 2;
	sisProps->pRadius = sisProps->imageDIMy/sisProps->cubeDimension.y*partProps->radius;

	// Bloco inicial de esferas
	float corner1[2] = {0.1, 0.1};
	float corner2[2] = {9.9, 9.9};
	float sideLenght[2];

	// Grid dimension
	uint grid = sisProps->cubeDimension.x / (4.0f * partProps[0].radius);
	uint temp = log2(grid);
	uint gridUpdate = pow(2,temp);
	float cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;
	sisProps->gridSize.x = pow(2,temp);
	
	grid = sisProps->cubeDimension.y / (4 * partProps[0].radius);
	temp = log2(grid);
	gridUpdate = pow(2,temp);
	cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;	
	sisProps->gridSize.y = pow(2,temp);

	sisProps->numCells = sisProps->gridSize.x * sisProps->gridSize.y;

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
	cudaMalloc((void**)&partValues->pos1, sizeof(float2) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->pos2, sizeof(float2) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->vel1, sizeof(float2) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->vel2, sizeof(float2) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->acc, sizeof(float2) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->cellStart, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->cellEnd, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->gridParticleIndex, sizeof(uint) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->gridParticleHash, sizeof(uint) * sisProps->numParticles);
	cudaMemcpy(d_corner1, corner1, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sideLenght, sideLenght, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_side, side, sizeof(uint)*2, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(sisPropD, sisProps, sizeof(SistemProperties));
	cudaMemcpyToSymbol(partPropD, partProps , sizeof(ParticleProperties) * 1);

	initializeParticlePosition(partValues->pos1,
							   partValues->vel1,
							   partValues->acc,
							   d_corner1,
							   d_sideLenght,
							   d_side,
							   side);

    cudaFree( d_corner1 );
    cudaFree( d_sideLenght );
    cudaFree( d_side );


	// Screen output	
	printf("Number of Particles = %d\n",sisProps->numParticles);
	printf("grid %d x %d\n",sisProps->gridSize.x,sisProps->gridSize.y);
}

void SimLooping( uchar4 *image, DataBlock *simBlock, int ticks ) {

    SistemProperties *sisProps = &simBlock->sisProps;
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
			 partValues->gridParticleIndex,
			 partValues->gridParticleHash,
			 sisProps->numParticles);

	// Ordena o grid pela posicao das particulas
	sortParticles(partValues->gridParticleHash,
				  partValues->gridParticleIndex,
				  sisProps->numParticles);

	// Encontra as particulas de inicializacao e de finalizacao
	reorderDataAndFindCellStart(partValues->cellStart,
								partValues->cellEnd,
								sortPos,
								sortVel,
								partValues->gridParticleHash,
								partValues->gridParticleIndex,
								oldPos,
								oldVel,
								sisProps->numParticles,
								sisProps->numCells);

	// Detecta a colizao das particulas
	collide(sortPos,
			sortVel,
			partValues->acc,
			partValues->gridParticleIndex,
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
	if (ticks % NUM_FRAMES_PER_PLOT == 1){
#endif
		// Saida grarica quando necessario
		plotParticles(image,
					  sortPos,
					  sisProps->numParticles,
					  sisProps->imageDIMx,
					  sisProps->imageDIMy);

#if !ALWAYS
	}
#endif
}

void FinalizingSim( DataBlock *simBlock) {

    // Limpe aqui o que tiver que ser limpo
    
    SistemProperties *sisProps = &simBlock->sisProps;
    ParticleProperties *partProps = &simBlock->partProps;
    ParticlesValues *partValues = &simBlock->partValues;
    
    cudaFree( partValues->pos1 );
    cudaFree( partValues->pos2 );
    cudaFree( partValues->vel1 );
    cudaFree( partValues->vel2 );
    cudaFree( partValues->acc );
    cudaFree( partValues->cellStart );
    cudaFree( partValues->cellEnd );
    cudaFree( partValues->gridParticleIndex );
    cudaFree( partValues->gridParticleHash );
    cudaFree( sisProps );
    cudaFree( partValues );
    cudaFree( partProps );

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
