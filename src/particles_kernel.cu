// Esse arquivo contem todas as funções executadas na placa de vídeo

// Define se pega a variável da memória global ou da memória de textura
#if USE_TEX
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif

#if USE_TEX
// textures for particle position and velocity
texture<float2, 1, cudaReadModeElementType> oldPosTex;
texture<float2, 1, cudaReadModeElementType> oldVelTex;

//texture<uint, 1, cudaReadModeElementType> gridParticleHashTex;
texture<uint, 1, cudaReadModeElementType> cellStartTex;
texture<uint, 1, cudaReadModeElementType> cellEndTex;
#endif

// Declarando as variáveis da memória de constante
__constant__ SistemProperties sisPropD;
__constant__ ParticleProperties partPropD;

__global__ void initializeParticlePositionD(float2*			pos,
											float2*			vel,
											float2*			acc,
											float*			corner1,
											float*			comp,
											uint*			side,
											unsigned long 	seed) {
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x >= side[0]) return;
	if (y >= side[1]) return;
	
	uint particle = x + y*side[0];

    if (particle >= sisPropD.numParticles) return;

	curandState state;
	curand_init( seed, particle, 0, &state );
	
//	float comp[2];
//	
//	comp[0] = corner2[0] - corner1[0];
//	comp[1] = corner2[1] - corner1[1];

	pos[particle].x = corner1[0] + comp[0]/(side[0]-1) * (x + (curand_normal(&state)-0.5f)/100);
	pos[particle].y = corner1[1] + comp[1]/(side[1]-1) * (y + (curand_normal(&state)-0.5f)/100);
	vel[particle] = make_float2( 0 );
	acc[particle] = make_float2( 0 );

}

// calculate position in uniform grid
__device__ int2 calcGridPos(float2 p)
{
    int2 gridPos;
    gridPos.x = floor(p.x * sisPropD.gridSize.x / sisPropD.cubeDimension.x);
    gridPos.y = floor(p.y * sisPropD.gridSize.y / sisPropD.cubeDimension.y);
    return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHash(int2 gridPos)
{
    gridPos.x = gridPos.x & (sisPropD.gridSize.x-1);  // wrap grid, assumes size is power of 2
    gridPos.y = gridPos.y & (sisPropD.gridSize.y-1);
    return gridPos.y * sisPropD.gridSize.x + gridPos.x;
}

// calculate grid hash value for each particle
__global__
void calcHashD(uint*   gridParticleHash,  // output
               uint*   gridParticleIndex,
               float2* pos)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (index >= sisPropD.numParticles) return;
    
    float2 p = pos[index];

//    // get address in grid
    int2 gridPos = calcGridPos(p);
    uint hash = calcGridHash(gridPos);

//    // store grid hash and particle index
    gridParticleHash[index] = hash;
    gridParticleIndex[index] = index;
}


// rearrange particle data into sorted order, and find the start of each cell
// in the sorted hash array
__global__
void reorderDataAndFindCellStartD(uint*   cellStart,        // output: cell start index
							      uint*   cellEnd,          // output: cell end index
							      float2* sortedPos,        // output: sorted positions
  							      float2* sortedVel,        // output: sorted velocities
                                  uint *  gridParticleHash, // input: sorted grid hashes
                                  uint *  gridParticleIndex,// input: sorted particle indices
				                  float2* oldPos,           // input: sorted position array
							      float2* oldVel)           // input: sorted velocity array
{
	extern __shared__ uint sharedHash[];    // blockSize + 1 elements
    uint index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;
	
    uint hash;
    // handle case when no. of particles not multiple of block size
    if (index < sisPropD.numParticles) {
        hash = gridParticleHash[index];

        // Load hash data into shared memory so that we can look 
        // at neighboring particle's hash value without loading
        // two hash values per thread
	    sharedHash[threadIdx.x+1] = hash;

	    if (index > 0 && threadIdx.x == 0)
	    {
		    // first thread in block must load neighbor particle hash
		    sharedHash[0] = gridParticleHash[index-1];
	    }
	}

	__syncthreads();
	
	if (index < sisPropD.numParticles) {
		// If this particle has a different cell index to the previous
		// particle then it must be the first particle in the cell,
		// so store the index of this particle in the cell.
		// As it isn't the first particle, it must also be the cell end of
		// the previous particle's cell

	    if (index == 0 || hash != sharedHash[threadIdx.x])
	    {
		    cellStart[hash] = index;
            if (index > 0)
                cellEnd[sharedHash[threadIdx.x]] = index;
	    }

        if (index == sisPropD.numParticles - 1)
        {
            cellEnd[hash] = index + 1;
        }

	    // Now use the sorted index to reorder the pos and vel data
	    uint sortedIndex = gridParticleIndex[index];
	    float2 pos = FETCH(oldPos, sortedIndex);       // macro does either global read or texture fetch
        float2 vel = FETCH(oldVel, sortedIndex);       // see particles_kernel.cuh

        sortedPos[index] = pos;
        sortedVel[index] = vel;
	}
}

// collide two spheres using DEM method
__device__
float2 collideSpheres(float2 posA, float2 posB,
                      float2 velA, float2 velB,
                      float radiusA, float radiusB)
{
	// calculate relative position
    float2 relPos = posB - posA;

    float dist = length(relPos);
    float collideDist = radiusA + radiusB;

    float2 force = make_float2(0.0f);
    if (dist < collideDist) {
        float2 norm = relPos / dist;

		// relative velocity
        float2 relVel = velB - velA;

        // relative tangential velocity
//        float2 tanVel = relVel - (dot(relVel, norm) * norm);

        // spring force
        force = -partPropD.collideStiffness*(collideDist - dist) * norm;
        // dashpot (damping) force
        force += partPropD.collideDamping*relVel;
        // tangential shear force
//        force += params.shear*tanVel;
    }

    return force;
}



// collide a particle against all other particles in a given cell
__device__
float2 collideCell(int2    gridPos,
                   uint    index,
                   float2  pos,
                   float2  vel,
                   float2* oldPos, 
                   float2* oldVel,
                   uint*   cellStart,
                   uint*   cellEnd)
{
    uint gridHash = calcGridHash(gridPos);

    // get start of bucket for this cell
    uint startIndex = FETCH(cellStart, gridHash);

    float2 force = make_float2(0.0f);
    if (startIndex != 0xffffffff) {        // cell is not empty
        // iterate over particles in this cell
        uint endIndex = FETCH(cellEnd, gridHash);
        for(uint j=startIndex; j<endIndex; j++) {
            if (j != index) {              // check not colliding with self
	            float2 pos2 = FETCH(oldPos, j);
                float2 vel2 = FETCH(oldVel, j);

                // collide two spheres
                force += collideSpheres(pos, pos2, vel, vel2, partPropD.radius, partPropD.radius);
            }
        }
    }
    return force;
}


__global__
void collideD(float2* oldPos,               // input: sorted positions
              float2* oldVel,               // input: sorted velocities
              float2* newAcc,                // output: new acceleration
              uint*   cellStart,
              uint*   cellEnd)
{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;    
    
    // read particle data from sorted arrays
	float2 pos = FETCH(oldPos, index);
    float2 vel = FETCH(oldVel, index);

    // get address in grid
    int2 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float2 force = make_float2(0.0f);
    for(int y=-1; y<=1; y++) {
        for(int x=-1; x<=1; x++) {
            int2 neighbourPos = gridPos + make_int2(x, y);
            force += collideCell(neighbourPos, index, pos, vel, oldPos, oldVel, cellStart, cellEnd);
        }
    }

	newAcc[index] = force / partPropD.mass;
}

__global__
void integrateSystemD(float2* pos, float2* vel, float2* acc)
{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;  
    
    vel[index] += (sisPropD.gravity + acc[index]) * sisPropD.timeStep;
    pos[index] += vel[index] * sisPropD.timeStep;
    
        // set this to zero to disable collisions with cube sides
#if 1
        if (pos[index].x > sisPropD.cubeDimension.x - partPropD.radius) {
        	pos[index].x = sisPropD.cubeDimension.x - partPropD.radius;
        	vel[index].x *= partPropD.boundaryDamping; }
        if (pos[index].x < partPropD.radius){
        	pos[index].x = partPropD.radius;
        	vel[index].x *= partPropD.boundaryDamping;}
        if (pos[index].y > sisPropD.cubeDimension.x - partPropD.radius) { 
        	pos[index].y = sisPropD.cubeDimension.x - partPropD.radius;
        	vel[index].y *= partPropD.boundaryDamping; }
#endif  
        if (pos[index].y < partPropD.radius) {
        	pos[index].y = partPropD.radius;
        	vel[index].y *= partPropD.boundaryDamping;}
}


// Pinta a esfera.
__global__
void plotSpheresD(uchar4*	ptr,
				  float2* 	sortPos)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (index >= sisPropD.numParticles) return;
	
	float2 pos = sortPos[index];
	
	// calcula a posição do centro da partícula em pixel
	int cPixelx = sisPropD.imageDIMx/sisPropD.cubeDimension.x*pos.x;
	int cPixely = sisPropD.imageDIMy/sisPropD.cubeDimension.y*pos.y;
	
	// percorre o quadrado ocupado pela partícula (em pixel)
	for (int x = -sisPropD.dimx/2; x < sisPropD.dimx/2; x++ ) {
		for (int y = -sisPropD.dimy/2; y < sisPropD.dimy/2; y++) {
			if (x*x + y*y < sisPropD.pRadius*sisPropD.pRadius) {
				
				// posição do ponto atual (em pixel)
				uint gPixelx = cPixelx + x;
				uint gPixely = cPixely + y;
				
				// Cria o efeito 3D da partícula (escurece as bordas)
				float fscale = sqrtf((sisPropD.pRadius*sisPropD.pRadius - x*x - y*y)/(sisPropD.pRadius*sisPropD.pRadius));
				
				// posição do pixel no vetor da imagem
				uint pixel = gPixelx + gPixely*sisPropD.imageDIMx;
				if (pixel >= sisPropD.imageDIMx*sisPropD.imageDIMy) pixel = sisPropD.imageDIMx*sisPropD.imageDIMy-1;
				
				// define a cor do pixel
				ptr[pixel].x = 255.0f * fscale;
				ptr[pixel].y = 255.0f * fscale;
				ptr[pixel].z = 255.0f * fscale;
				ptr[pixel].w = 255.0f * fscale;
				
			}
		}
	}
}
