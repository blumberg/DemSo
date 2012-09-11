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

// Esse arquivo prepara as funções que serão executadas na GPU. Ele define
// o tamanho do Grid e o número de Threads.

// Função para retornar o maior inteiro da divisão a/b
uint iDivUp(uint a, uint b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

// compute grid and thread block size for a given number of elements
void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
{
    numThreads = min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}

// cria a posição inicial das partículas. Esse kernel é executado em um
// grid 2D com um número máximo de 16 threads por bloco
void initializeParticlePosition (float* 		pos,
								 float* 		vel,
								 float* 		acc,
								 float*			corner1,
								 float*			comp,
								 uint*			side,
								 uint*			hostSide,
								 unsigned long	seed){

	uint numBlocksx, numBlocksy, numThreadsx, numThreadsy;
	computeGridSize(hostSide[0], 16, numBlocksx, numThreadsx);
	computeGridSize(hostSide[1], 16, numBlocksy, numThreadsy);
	
	dim3 numBlocks(numBlocksx,numBlocksy);
	dim3 numThreads(numThreadsx,numThreadsy);

initializeParticlePositionD<<<numBlocks,numThreads>>>((float2*)pos,
													  (float2*)vel,
													  (float2*)acc,
													  corner1,
													  comp,
													  side,
													  seed);
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
void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
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
                                 uint*  gridParticleHash,
                                 uint*  gridParticleIndex,
							     float* oldPos,
							     float* oldVel,
							     uint   numParticles,
							     uint   numCells)
{
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 256, numBlocks, numThreads);

    // set all cells to empty
	cudaMemset(cellStart, 0xffffffff, numCells*sizeof(uint));

	// Declarando como memória de textura
	#if USE_TEX
		cudaBindTexture(0, oldPosTex, oldPos, numParticles*sizeof(float2));
		cudaBindTexture(0, oldVelTex, oldVel, numParticles*sizeof(float2));
	#endif

    uint smemSize = sizeof(uint)*(numThreads+1);
    reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(
        cellStart,
        cellEnd,
        (float2*)sortedPos,
        (float2*)sortedVel,
		gridParticleHash,
		gridParticleIndex,
        (float2*)oldPos,
        (float2*)oldVel);
    
    // Retirando da memória de textura 
	#if USE_TEX
		cudaUnbindTexture(oldPosTex);
		cudaUnbindTexture(oldVelTex);
	#endif

}

// Rotina que verifica a colisão entre as partículas e transforma a força
// de colisão em aceleração. Esse kernel é executado em um grid 1D com um
// número máximo de 64 threads por bloco
void collide(float* 	oldPos,
             float* 	oldVel,
             float* 	newAcc,
             uint*  	cellStart,
             uint*  	cellEnd,
             uint   	numParticles,
             uint 		numCells)
{
	// Declarando como memória de textura
	#if USE_TEX
		cudaBindTexture(0, oldPosTex, sortedPos, numParticles*sizeof(float2));
		cudaBindTexture(0, oldVelTex, sortedVel, numParticles*sizeof(float2));
		cudaBindTexture(0, cellStartTex, cellStart, numCells*sizeof(uint));
		cudaBindTexture(0, cellEndTex, cellEnd, numCells*sizeof(uint));    
	#endif

    // thread per particle
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 64, numBlocks, numThreads);

    // execute the kernel
    collideD<<< numBlocks, numThreads >>>((float2*)oldPos,
                                          (float2*)oldVel,
                                          (float2*)newAcc,
                                          cellStart,
                                          cellEnd);
    // Retirando da memória de textura 
	#if USE_TEX
		cudaUnbindTexture(oldPosTex);
		cudaUnbindTexture(oldVelTex);
		cudaUnbindTexture(cellStartTex);
		cudaUnbindTexture(cellEndTex);
	#endif
}

// Realiza a integração numérica do sistema. Essa é uma integração linear,
// onde:
// Velocidade = Velocidade + Aceleração * DeltaTempo
// Posicão    =  Posição   + Velocidade * DeltaTempo
// Esse kernel é executado em um grid 1D com um número máximo de 256
// threads por bloco.
void integrateSystem(float *pos,
					 float *vel,
					 float *acc,
					 uint numParticles)
{
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	integrateSystemD<<<numBlocks,numThreads>>>((float2*)pos,
				 							   (float2*)vel,
				 							   (float2*)acc);
}

// Desenha as partículas em uma imagem de DIMx x DIMy pixels e mostra na
// tela. O fundo da imagem é definido como preto e as partículas são
// brancas. Esse kernel é executado em um grid 1D com um número máximo de
// 256 threads por bloco.
void plotParticles(uchar4*	ptr,
				   float* 	pos,
				   uint 	numParticles,
				   int 		DIMx,
				   int		DIMy){

	// pinta o fundo de preto
	cudaMemset(ptr, 0, DIMx*DIMy*sizeof(uchar4));
	
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	plotSpheresD<<<numBlocks,numThreads>>>(ptr,
									 	   (float2*)pos);

}
