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

#include <iostream>
#include <vector>
#include <thrust/device_ptr.h>  		   // thrust para utilizar ponteiro
#include <thrust/sort.h> 					   // thrust para ordenar vetor
#include <thrust/extrema.h>
#include "main.cuh"
#include "functions.cuh"
#include "particles_kernel.cuh"

using std::cout;
using std::endl;
using std::vector;

// Esse arquivo prepara as funções que serão executadas na GPU. Ele define
// o tamanho do Grid e o número de Threads.


// Função para retornar o maior inteiro da divisão a/b
inline uint iDivUp(uint a, uint b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

// << Operator overloadings for cout
std::ostream& operator<<(std::ostream& output, const uint2& p)
{
	output << "(" << p.x << ", " << p.y << ")";
	return output;
}
std::ostream& operator<<(std::ostream& output, const float2& p)
{
	output << "(" << p.x << ", " << p.y << ")";
	return output;
}
std::ostream& operator<<(std::ostream& output, const float3& p)
{
	output << "(" << p.x << ", " << p.y << ", " << p.z << ")";
	return output;
}

// compute grid and thread block size for a given number of elements
void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
{
    numThreads = min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}

// Aloca espaço na memória da GPU e copia as propriedades para a memória de
// constantes.
void allocateVectors(ParticleProperties* partProps,
					 ParticlesValues* partValues,
					 SystemProperties* sisProps,
					 RenderParameters* renderPar)
{
	uint numParticles = sisProps->numParticles;
	// alocando vetores na placa de video
	// cudaMalloc --> aloca espaço na placa de vídeo
	// cudaMemcpy --> transfere dados entre a CPU (Host) e GPU (Device)
	// cudaMemcpyToSymbol --> copia variável para a memória de constante
	cudaMalloc((void**)&partValues->type1, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->type2, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->ID1, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->ID2, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->loc1, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->loc2, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->pos1, sizeof(float) * numParticles * 2);
	cudaMalloc((void**)&partValues->pos2, sizeof(float) * numParticles * 2);
	cudaMalloc((void**)&partValues->vel1, sizeof(float) * numParticles * 2);
	cudaMalloc((void**)&partValues->vel2, sizeof(float) * numParticles * 2);
	cudaMalloc((void**)&partValues->acc, sizeof(float) * numParticles * 2);
	cudaMalloc((void**)&partValues->theta1, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->theta2, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->omega1, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->omega2, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->alpha, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->cellStart, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->cellEnd, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->gridParticleIndex, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->gridParticleHash, sizeof(uint) * numParticles);
	cudaMalloc((void**)&partValues->pressure, sizeof(float) * numParticles);

#if USE_BIG_PARTICLE
#if USE_ATOMIC
	cudaMalloc((void**)&partValues->controlForce, sizeof(float2));
	cudaMalloc((void**)&partValues->controlMoment, sizeof(float));
#else
	uint numBlocks, numThreads;
	computeGridSize(numParticles, 512, numBlocks, numThreads);
	cudaMalloc((void**)&partValues->controlForceVecX, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->controlForceVecY, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->controlMomentVec, sizeof(float) * numParticles);
	cudaMalloc((void**)&partValues->hCFVx, sizeof(float) * numBlocks);
	cudaMalloc((void**)&partValues->hCFVy, sizeof(float) * numBlocks);
	cudaMalloc((void**)&partValues->hCMV, sizeof(float) * numBlocks);
	cudaMemset(partValues->hCFVx, 0, sizeof(float) * numBlocks);
	cudaMemset(partValues->hCFVy, 0, sizeof(float) * numBlocks);
	cudaMemset(partValues->hCMV, 0, sizeof(float) * numBlocks);
#endif // USE_ATOMIC
	partValues->ctrlF = (float2*)malloc(sizeof(float2));
	partValues->ctrlM = (float*)malloc(sizeof(float));
#endif // USE_BIG_PARTICLE

	
	// Definindo 0 como valor inicial de todos os vetores alocados acima
	cudaMemset(partValues->type1, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->type2, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->ID1, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->ID2, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->loc1, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->loc2, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->pos1, 0, sizeof(float) * numParticles * 2);
	cudaMemset(partValues->pos2, 0, sizeof(float) * numParticles * 2);
	cudaMemset(partValues->vel1, 0, sizeof(float) * numParticles * 2);
	cudaMemset(partValues->vel2, 0, sizeof(float) * numParticles * 2);
	cudaMemset(partValues->acc, 0, sizeof(float) * numParticles * 2);
	cudaMemset(partValues->theta1, 0, sizeof(float) * numParticles);
	cudaMemset(partValues->theta2, 0, sizeof(float) * numParticles);
	cudaMemset(partValues->omega1, 0, sizeof(float) * numParticles);
	cudaMemset(partValues->omega2, 0, sizeof(float) * numParticles);
	cudaMemset(partValues->alpha, 0, sizeof(float) * numParticles);
	cudaMemset(partValues->cellStart, 0, sizeof(uint) * sisProps->numCells);
	cudaMemset(partValues->cellEnd, 0, sizeof(uint) * sisProps->numCells);
	cudaMemset(partValues->gridParticleIndex, 0, sizeof(uint) * numParticles);
	cudaMemset(partValues->gridParticleHash, 0, sizeof(uint) * numParticles);


	cudaMemcpyToSymbol(sisPropD, sisProps, sizeof(SystemProperties));
	cudaMemcpyToSymbol(renderParD, renderPar, sizeof(RenderParameters));
	cudaMemcpyToSymbol(partPropD, partProps , sizeof(ParticleProperties) * MAX_PARTICLES_TYPES);

}



// Desaloca o espaço reservado na GPU
void desAllocateVectors(ParticlesValues* partValues)
{
	cudaFree( partValues->type1 );
	cudaFree( partValues->type2 );
	cudaFree( partValues->ID1 );
	cudaFree( partValues->ID2 );
	cudaFree( partValues->loc1 );
	cudaFree( partValues->loc2 );
    cudaFree( partValues->pos1 );
    cudaFree( partValues->pos2 );
    cudaFree( partValues->vel1 );
    cudaFree( partValues->vel2 );
    cudaFree( partValues->acc );
    cudaFree( partValues->theta1 );
    cudaFree( partValues->theta2 );
    cudaFree( partValues->omega1 );
    cudaFree( partValues->omega2 );
    cudaFree( partValues->alpha );
    cudaFree( partValues->cellStart );
    cudaFree( partValues->cellEnd );
    cudaFree( partValues->gridParticleIndex );
    cudaFree( partValues->gridParticleHash );
    cudaFree( partValues->pressure );

#if USE_BIG_PARTICLE
#if USE_ATOMIC
	cudaFree( partValues->controlForce );
	cudaFree( partValues->controlMoment );
#else
	cudaFree( partValues->controlForceVecX );
	cudaFree( partValues->controlForceVecY );
	cudaFree( partValues->controlMomentVec );
	cudaFree( partValues->hCFVx );
	cudaFree( partValues->hCFVy );
	cudaFree( partValues->hCMV );
#endif // USE_ATOMIC
	free( partValues->ctrlF );
	free( partValues->ctrlM );
#endif // USE_BIG_PARTICLE

}

// cria a posição inicial das partículas. Esse kernel é executado em um
// grid 2D com um número máximo de 16 threads por bloco
void createRectangles (float* 		pos,
					  uint*			ID,
					  uint*			loc,
					  uint*			type,
					  float2		start,
					  float2		sideLenght,
					  uint2			side,
					  uint 			startID,
					  uint 			numParticleTypes,
					  uint*			particleTypeVec,
					  unsigned long	seed)
{
	// alocando vetores na placa de video
	// cudaMalloc --> aloca espaço na placa de vídeo
	// cudaMemcpy --> transfere dados entre a CPU (Host) e GPU (Device)
	// cudaMemcpyToSymbol --> copia variável para a memória de constante
	uint *d_particleTypeVec;

	cudaMalloc((void**)&d_particleTypeVec, sizeof(uint)*numParticleTypes);

	cudaMemcpy(d_particleTypeVec, particleTypeVec,
			   sizeof(uint)*numParticleTypes, cudaMemcpyHostToDevice);

	uint numBlocksx, numBlocksy, numThreadsx, numThreadsy;
	computeGridSize(side.x, 16, numBlocksx, numThreadsx);
	computeGridSize(side.y, 16, numBlocksy, numThreadsy);
	
	dim3 numBlocks(numBlocksx,numBlocksy);
	dim3 numThreads(numThreadsx,numThreadsy);

	createRectanglesD<<<numBlocks,numThreads>>>((float2*)pos,
											   ID,
											   loc,
											   type,
											   start,
											   sideLenght,
											   side,
											   startID,
											   numParticleTypes,
											   d_particleTypeVec,
											   seed);

	// Desalocando espaço na placa de vídeo (Não mais necessário)
    cudaFree( d_particleTypeVec );													  
}

void createTriangles (float*	pos,
					  uint*		ID,
					  uint*		loc,
					  uint*		type,
					  float2	start,
					  uint		N,
					  uint		numParticleTypes,
					  uint*		particleTypeVec,
					  float 	space,
					  float		height,
					  uint 		startID,
					  uint		numParticles)
{
	uint *d_particleTypeVec;

	cudaMalloc((void**)&d_particleTypeVec, sizeof(uint)*numParticleTypes);

	cudaMemcpy(d_particleTypeVec, particleTypeVec,
			   sizeof(uint)*numParticleTypes, cudaMemcpyHostToDevice);

	uint numBlocks, numThreads;
	computeGridSize(numParticles,256,numBlocks,numThreads);
	
	createTrianglesD<<<numBlocks,numThreads>>>((float2*)pos,
											   ID,
											   loc,
											   type,
											   start,
											   N,
											   numParticleTypes,
											   d_particleTypeVec,
											   space,
											   height,
											   startID,
											   numParticles);
	
    cudaFree( d_particleTypeVec );													  
}

void createSingleParticles (float*	pos,
							float*	vel,
							float*	theta,
							float*	omega,
							uint*	ID,
							uint*	loc,
							uint*	type,
							float2*	usrPos,
							float2* usrVel,
							float*	usrTheta,
							float*	usrOmega,
							uint*	usrType,
							uint	numParticles,
							uint	startID){

	uint IDvec[numParticles];
	for (int i = 0; i < numParticles ; i++){
		IDvec[i] = i + startID;
	}
	
	cudaMemcpy(pos,usrPos,sizeof(float)*numParticles*2,cudaMemcpyHostToDevice);
	cudaMemcpy(vel,usrVel,sizeof(float)*numParticles*2,cudaMemcpyHostToDevice);
	cudaMemcpy(theta,usrTheta,sizeof(float)*numParticles,cudaMemcpyHostToDevice);
	cudaMemcpy(omega,usrOmega,sizeof(float)*numParticles,cudaMemcpyHostToDevice);
	cudaMemcpy(type,usrType,sizeof(uint)*numParticles,cudaMemcpyHostToDevice);
	cudaMemcpy(ID,IDvec,sizeof(uint)*numParticles,cudaMemcpyHostToDevice);
	cudaMemcpy(loc,IDvec,sizeof(uint)*numParticles,cudaMemcpyHostToDevice);
							
}

// Calcula o numero da celula de cada particula. Esse kernel é executado
// em um grid 1D com um número máximo de 256 threads por bloco
void calcHash(float* 	pos,
			  uint* 	gridParticleIndex,
			  uint* 	gridParticleHash,
			  uint 		numParticles)
{
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 256, numBlocks, numThreads);

    // execute the kernel
    calcHashD<<< numBlocks, numThreads >>>(gridParticleHash,
                                           gridParticleIndex,
                                           (float2*)pos);
}

// Ordena as partículas com base no número do Hash. Essa rotina é executada
// pela biblioteca THRUST.
// A função device_ptr permite passar o ponteiro de uma variável alocada na
// GPU para o thrust.
// Em seguida a função sort_by_key organiza o vetor dGridParticleHash em
// ordem crescente e arruma o vetor dGridParticleIndex com base na
// ordenação
void sortParticles(uint* dGridParticleHash, uint* dGridParticleIndex, uint numParticles)
{
    thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash),
                        thrust::device_ptr<uint>(dGridParticleHash + numParticles),
                        thrust::device_ptr<uint>(dGridParticleIndex));
}

// Reordena os vetores de posição e velocidade com base na nova ordem das
// partículas. Em seguida o vetor de inicio e fim de cada célula é criado.
// Esse kernel é executado em um grid 1D com um número máximo de 256
// threads por bloco
void reorderDataAndFindCellStart(uint*  cellStart,
							     uint*  cellEnd,
							     float* sortedPos,
							     float* sortedVel,
								 float* sortedTheta,
								 float* sortedOmega,
							     uint* 	sortedID,
							     uint* 	sortedLoc,
							     uint* 	sortedType,
                                 uint*  gridParticleHash,
                                 uint*  gridParticleIndex,
							     float* oldPos,
							     float* oldVel,
								 float* oldTheta,
								 float* oldOmega,
							     uint*	oldID,
							     uint* 	oldType,
							     uint   numParticles,
							     uint   numCells)
{
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 256, numBlocks, numThreads);

    // set all cells to empty
	cudaMemset(cellStart, 0xffffffff, numCells*sizeof(uint));

    uint smemSize = sizeof(uint)*(numThreads+1);
    reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(
        cellStart,
        cellEnd,
        (float2*)sortedPos,
        (float2*)sortedVel,
		sortedTheta,
		sortedOmega,
        sortedID,
        sortedLoc,
        sortedType,
		gridParticleHash,
		gridParticleIndex,
        (float2*)oldPos,
        (float2*)oldVel,
		oldTheta,
		oldOmega,
        oldID,
        oldType);

}

// Rotina que verifica a colisão entre as partículas e transforma a força
// de colisão em aceleração. Esse kernel é executado em um grid 1D com um
// número máximo de 64 threads por bloco
void collide(float* 	oldPos,
             float* 	oldVel,
             float* 	newAcc,
			 float*		oldOmega,
			 float*		newAlpha,
             uint*		oldType,
             uint*  	cellStart,
             uint*  	cellEnd,
             uint   	numParticles,
             uint 		numCells,
#if USE_BIG_PARTICLE
			 float2		controlPos,
			 float2		controlVel,
			 float 		controlTheta,
			 float		controlOmega,
			 uint		controlType,
#if USE_ATOMIC
			 float2*	controlForce,
			 float*		controlMoment,
#else
			 float*		controlForceVecX,
			 float*		hCFVx,
			 float*		controlForceVecY,
			 float*		hCFVy,
			 float*		controlMomentVec,
			 float*		hCMV,
#endif // USE_ATOMIC
			 float2*	ctrlF,
			 float*		ctrlM,
#endif // USE_BIG_PARTICLE
			 float*		pressure)
{

    // thread per particle
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 64, numBlocks, numThreads);


#if (USE_ATOMIC*USE_BIG_PARTICLE)
	cudaMemset(controlForce, 0, sizeof(float2));
	cudaMemset(controlMoment, 0, sizeof(float));
#endif

	cudaMemset(pressure,0,sizeof(float)*numParticles);

    // execute the kernel
    collideD<<< numBlocks, numThreads >>>((float2*)oldPos,
                                          (float2*)oldVel,
                                          (float2*)newAcc,
										  oldOmega,
										  newAlpha,
                                          oldType,
                                          cellStart,
                                          cellEnd,
#if USE_BIG_PARTICLE
			 							  controlPos,
			 							  controlVel,
			 							  controlTheta,
			 							  controlOmega,
										  controlType,
#if USE_ATOMIC
										  controlForce,
										  controlMoment,
#else
										  controlForceVecX,
										  controlForceVecY,
										  controlMomentVec,
#endif // USE_ATOMIC
#endif // USE_BIG_PARTICLE
										  pressure);


#if USE_BIG_PARTICLE
#if USE_ATOMIC
	cudaMemcpy(ctrlF,controlForce,sizeof(float2),cudaMemcpyDeviceToHost);
	cudaMemcpy(ctrlM,controlMoment,sizeof(float),cudaMemcpyDeviceToHost);	
#else
	computeGridSize(numParticles, 512, numBlocks, numThreads);
	uint smemSize = sizeof(float)*numThreads;
	
	cudaStream_t stream0, stream1, stream2;
	float cpu_hCFVx[numBlocks], cpu_hCFVy[numBlocks], cpu_hCMV[numBlocks];
	
	cudaStreamCreate( &stream0 );
	cudaStreamCreate( &stream1 );
	cudaStreamCreate( &stream2 );

	reduceVecD<<<numBlocks, numThreads, smemSize, stream0>>>(controlForceVecX,
															 hCFVx);

	cudaMemcpyAsync(&cpu_hCFVx,hCFVx,sizeof(float)*numBlocks,cudaMemcpyDeviceToHost,stream0);

	reduceVecD<<<numBlocks, numThreads, smemSize, stream0>>>(controlForceVecY,
															 hCFVy);

	cudaMemcpyAsync(&cpu_hCFVy,hCFVy,sizeof(float)*numBlocks,cudaMemcpyDeviceToHost,stream1);

	reduceVecD<<<numBlocks, numThreads, smemSize, stream0>>>(controlMomentVec,
															 hCMV);																										

	cudaMemcpyAsync(&cpu_hCMV,hCMV,sizeof(float)*numBlocks,cudaMemcpyDeviceToHost,stream2);

	cudaStreamSynchronize( stream0 );
	cudaStreamSynchronize( stream1 );
	cudaStreamSynchronize( stream2 );

	cudaStreamDestroy( stream0 );
	cudaStreamDestroy( stream1 );
	cudaStreamDestroy( stream2 );

	ctrlF->x = 0;
	ctrlF->y = 0;
	ctrlM[0] = 0;

	for (int i = 0; i < numBlocks; i++){
		ctrlF->x += cpu_hCFVx[i];
		ctrlF->y += cpu_hCFVy[i];
		ctrlM[0] += cpu_hCMV[i];
	}

#endif // USE_ATOMIC
#endif // USE_BIG_PARTICLE

}

// Realiza a integração numérica do sistema. Essa é uma integração linear,
// onde:
// Velocidade = Velocidade + Aceleração * DeltaTempo
// Posicão    =  Posição   + Velocidade * DeltaTempo
// Esse kernel é executado em um grid 1D com um número máximo de 256
// threads por bloco.
void integrateSystem(float*	pos,
					 float*	vel,
					 float*	acc,
					 float* theta,
					 float* omega,
					 float* alpha,
					 uint*	type,
					 uint	numParticles)
{
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	integrateSystemD<<<numBlocks,numThreads>>>((float2*)pos,
				 							   (float2*)vel,
				 							   (float2*)acc,
											   theta,
											   omega,
											   alpha,
				 							   type);
}

// Desenha as partículas em uma imagem de DIMx x DIMy pixels e mostra na
// tela. O fundo da imagem é definido como preto e as partículas são
// brancas. Esse kernel é executado em um grid 1D com um número máximo de
// 256 threads por bloco.
void plotParticles(uchar4*	ptr,
				   float* 	pos,
				   float*	theta,
				   uint*	type,
				   uint 	numParticles,
				   int 		DIMx,
				   int		DIMy,
				   int		bgColor,
#if USE_BIG_PARTICLE
				   float2 	controlPos,
				   uint		controlType,
				   int		dimx,
				   int		dimy,
#endif
				   float*	pressure){

	// pinta o fundo de preto
	cudaMemset(ptr, bgColor, DIMx*DIMy*sizeof(uchar4));
	
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	plotSpheresD<<<numBlocks,numThreads>>>(ptr,
									 	   (float2*)pos,
										   theta,
									 	   type,
										   pressure);

#if USE_BIG_PARTICLE
	uint numBlocksx, numBlocksy, numThreadsx, numThreadsy;
	computeGridSize(dimx, 16, numBlocksx, numThreadsx);
	computeGridSize(dimy, 16, numBlocksy, numThreadsy);
	
	dim3 numBlocks2(numBlocksx,numBlocksy);
	dim3 numThreads2(numThreadsx,numThreadsy);

	// execute the kernel
	plotControlParticleD<<<numBlocks2,numThreads2>>>(ptr,
													 controlPos,
													 controlType);
#endif
}

// Escreve no arquivo de saída os dados desejados.
// O arquivo de saída é do tipo texto. Na primeira linha encontra-se
// o valor do timeStep. Em seguida, cada linha apresenta, separados
// por vírgulas, o número da iteracão, e cada um dos dados de saída
// desejados.
// TODO:
//  - Checar: Aparecimento de varios NaN quando shearStiffness = 1000
void writeOutputFile (FILE*			outputFile,
					  vector<int>	chosenOnes,
					  float			h_elapsedTime,
					  float2* 		pos,
					  float2*		vel,
					  float2*		acc,
					  float*		theta,
					  float*		omega,
					  float*		alpha,
					  uint*			ID,
					  uint*			type,
					  uint*			loc,
					  float*		pressure)
{
	// Print current elapsed time
	fprintf (outputFile, "%f,", h_elapsedTime); // Don't print newline

	// Loop over followed particles
	for (register int i = 0; i < chosenOnes.size(); i++)
	{
		// Chosen particle's location in arrays
		uint h_loc;
		cudaMemcpy (&h_loc, &loc[chosenOnes[i]], sizeof(uint), cudaMemcpyDeviceToHost);

		// Particle Data
		float2 h_pos, h_vel, h_acc;
		float h_theta, h_omega, h_alpha, h_pressure;
		uint h_id, h_type;

		// Geting particle data from GPU
		cudaMemcpy (&h_pos.x, &pos[h_loc].x, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_pos.y, &pos[h_loc].y, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_vel.x, &vel[h_loc].x, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_vel.y, &vel[h_loc].y, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_acc.x, &acc[h_loc].x, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_acc.y, &acc[h_loc].y, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_theta, &theta[h_loc], sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_omega, &omega[h_loc], sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_alpha, &alpha[h_loc], sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_id,	  &ID[h_loc],    sizeof(uint),  cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_type,  &type[h_loc],  sizeof(uint),  cudaMemcpyDeviceToHost);
		cudaMemcpy (&h_pressure, &pressure[h_loc], sizeof(float), cudaMemcpyDeviceToHost);

		// Print particle data
		fprintf (outputFile, "%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", h_id, h_type,
				 h_pos.x, h_pos.y, h_vel.x, h_vel.y, h_acc.x, h_acc.y,
				 h_theta, h_omega, h_alpha, h_pressure);

		// If we're not at the last particle, print a comma
		if (i != chosenOnes.size()-1) fprintf (outputFile, ",");
	}
	// Go to next line for next time step
	fprintf (outputFile, "\n");
}

void set_gravity (SystemProperties *sisProps, float2 gravity)
{
	sisProps->gravity = gravity;
	cout << "Gravity updated. New value = " << gravity << endl;
	cudaMemcpyToSymbol (sisPropD, sisProps, sizeof(SystemProperties));
}

void set_viewRotations (RenderParameters *renderPar, bool viewRotations)
{
	renderPar->viewRotations = viewRotations;
	cout << "Display of rotations turned " << ((viewRotations) ? "ON" : "OFF") << endl;
	cudaMemcpyToSymbol (renderParD, renderPar, sizeof(RenderParameters));
}

void set_colorByPressure (RenderParameters *renderPar, bool colorByPressure)
{
	renderPar->colorByPressure = colorByPressure;
	cout << "Coloration of particles by pressure value turned " << ((colorByPressure) ? "ON" : "OFF") << endl;
	cudaMemcpyToSymbol (renderParD, renderPar, sizeof(RenderParameters));
}

void updatePressureScale (RenderParameters *renderPar, float *pressure, int numParts)
{
	thrust::device_ptr<float> result = thrust::max_element(thrust::device_ptr<float>(pressure),
										thrust::device_ptr<float>(pressure + numParts));
	renderPar->maxPressure = (*result > 0.0) ? *result : 1.0;

	cout << "Pressure scale updated. New maximum value = " << renderPar->maxPressure << endl;

	cudaMemcpyToSymbol (renderParD, renderPar, sizeof(RenderParameters));
}
