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
void initializeParticlePosition (float2* 		pos,
								 float2* 		vel,
								 float2* 		acc,
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

initializeParticlePositionD<<<numBlocks,numThreads>>>(pos,
													  vel,
													  acc,
													  corner1,
													  comp,
													  side,
													  seed);
}

// Calcula o numero da celula de cada particula. Esse kernel é executado
// em um grid 1D com um número máximo de 256 threads por bloco
void calcHash(float2* 	pos,
			  uint* 	gridParticleIndex,
			  uint* 	gridParticleHash,
			  uint 		numParticles)
{
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 256, numBlocks, numThreads);

    // execute the kernel
    calcHashD<<< numBlocks, numThreads >>>(gridParticleHash,
                                           gridParticleIndex,
                                           pos);
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
							     float2* sortedPos,
							     float2* sortedVel,
                                 uint*  gridParticleHash,
                                 uint*  gridParticleIndex,
							     float2* oldPos,
							     float2* oldVel,
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
        sortedPos,
        sortedVel,
		gridParticleHash,
		gridParticleIndex,
        oldPos,
        oldVel);

}

// Rotina que verifica a colisão entre as partículas e transforma a força
// de colisão em aceleração. Esse kernel é executado em um grid 1D com um
// número máximo de 64 threads por bloco
void collide(float2* 	sortPos,
             float2* 	sortVel,
             float2* 	newAcc,
             uint*  	gridParticleIndex,
             uint*  	cellStart,
             uint*  	cellEnd,
             uint   	numParticles,
             uint   	numCells)
{
    // thread per particle
    uint numThreads, numBlocks;
    computeGridSize(numParticles, 64, numBlocks, numThreads);

    // execute the kernel
    collideD<<< numBlocks, numThreads >>>(sortPos,
                                          sortVel,
                                          newAcc,
                                          gridParticleIndex,
                                          cellStart,
                                          cellEnd);

}

// Realiza a integração numérica do sistema. Essa é uma integração linear,
// onde:
// Velocidade = Velocidade + Aceleração * DeltaTempo
// Posicão    =  Posição   + Velocidade * DeltaTempo
// Esse kernel é executado em um grid 1D com um número máximo de 256
// threads por bloco.
void integrateSystem(float2 *pos,
					 float2 *vel,
					 float2 *acc,
					 uint numParticles)
{
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	integrateSystemD<<<numBlocks,numThreads>>>(pos, vel, acc);

}

// Desenha as partículas em uma imagem de DIMx x DIMy pixels e mostra na
// tela. O fundo da imagem é definido como preto e as partículas são
// brancas. Esse kernel é executado em um grid 1D com um número máximo de
// 256 threads por bloco.
void plotParticles(uchar4*	ptr,
				   float2* 	pos,
				   uint 	numParticles,
				   int 		DIMx,
				   int		DIMy){

	// pinta o fundo de preto
	cudaMemset(ptr, 0, DIMx*DIMy*sizeof(uchar4));
	
	uint numThreads, numBlocks;
	computeGridSize(numParticles, 256, numBlocks, numThreads);
	
	// execute the kernel
	plotSpheresD<<<numBlocks,numThreads>>>(ptr,
									 	   pos);

}
