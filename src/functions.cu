

uint iDivUp(uint a, uint b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

// compute grid and thread block size for a given number of elements
void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
{
    numThreads = min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}

void initializeParticlePosition (float2* 	pos,
								 float2* 	vel,
								 float2* 	acc,
								 float*		corner1,
								 float*		comp,
								 uint*		side,
								 uint*		hostSide){

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
													  side);
}

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

void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
{
    thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash),
                        thrust::device_ptr<uint>(dGridParticleHash + numParticles),
                        thrust::device_ptr<uint>(dGridParticleIndex));
}

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

void plotParticles(uchar4*	ptr,
				   float2* 	pos,
				   uint 	numParticles,
				   int 		DIMx,
				   int		DIMy){

	uint numThreadsx, numBlocksx, numThreadsy, numBlocksy;
	computeGridSize(DIMx, 16, numBlocksx, numThreadsx);
	computeGridSize(DIMy, 16, numBlocksy, numThreadsy);

	dim3 numBlocks(numBlocksx, numBlocksy);
	dim3 numThreads(numThreadsx, numThreadsy);
	
	plotBackgroundD<<<numBlocks,numThreads>>>(ptr);
	
	uint numThreads2, numBlocks2;
	computeGridSize(numParticles, 256, numBlocks2, numThreads2);
	
	// execute the kernel
	plotSpheresD<<<numBlocks2,numThreads2>>>(ptr,
											 pos);

}
