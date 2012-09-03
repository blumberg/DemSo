#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include "../includes/cutil_math.h"


__constant__ SistemProperties simPropD;
__constant__ ParticleProperties partPropD;

__global__ void initializeParticlePositionD(float2 *pos,
											float2 *vel,
											float2 *acc,
											float *corner1,
											float *corner2,
											uint *side) {
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
    uint particle = x + y * blockDim.x * gridDim.x;
	
	uint numParticles = simPropD.numParticles;
	
    if (particle < numParticles){
		
		curandState state;
		curand_init( 1234, particle, 0, &state );
		
		float comp[2];
		
		comp[0] = corner2[0] - corner1[0];
		comp[1] = corner2[1] - corner1[1];
    
    	pos[particle].x = corner1[0] + comp[0]/(side[0]-1) * (x + (curand_normal(&state)-0.5f)/50);
    	pos[particle].y = corner1[1] + comp[1]/(side[1]-1) * (y + (curand_normal(&state)-0.5f)/50);
    	vel[particle] = make_float2( 0 );
    	acc[particle] = make_float2( 0 );
    }

}

// calculate position in uniform grid
__device__ int2 calcGridPos(float2 p)
{
    int2 gridPos;
    gridPos.x = floor(p.x * simPropD.gridSize.x / simPropD.cubeDimension.x);
    gridPos.y = floor(p.y * simPropD.gridSize.y / simPropD.cubeDimension.y);
    return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHash(int2 gridPos)
{
    gridPos.x = gridPos.x & (simPropD.gridSize.x-1);  // wrap grid, assumes size is power of 2
    gridPos.y = gridPos.y & (simPropD.gridSize.y-1);
    return gridPos.y * simPropD.gridSize.x + gridPos.x;
}

// calculate grid hash value for each particle
__global__
void calcHashD(uint*   gridParticleHash,  // output
               uint*   gridParticleIndex,
               float2* pos)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    uint numParticles = simPropD.numParticles;
    
    if (index >= numParticles) return;
    
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
	uint numParticles = simPropD.numParticles;
	
    uint hash;
    // handle case when no. of particles not multiple of block size
    if (index < numParticles) {
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
	
	if (index < numParticles) {
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

        if (index == numParticles - 1)
        {
            cellEnd[hash] = index + 1;
        }

	    // Now use the sorted index to reorder the pos and vel data
	    uint sortedIndex = gridParticleIndex[index];
	    float2 pos = oldPos[sortedIndex];
        float2 vel = oldVel[sortedIndex];

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
    uint startIndex = cellStart[gridHash];

    float2 force = make_float2(0.0f);
    if (startIndex != 0xffffffff) {        // cell is not empty
        // iterate over particles in this cell
        uint endIndex = cellEnd[gridHash];
        for(uint j=startIndex; j<endIndex; j++) {
            if (j != index) {              // check not colliding with self
	            float2 pos2 = oldPos[j];
                float2 vel2 = oldVel[j];

                // collide two spheres
                force += collideSpheres(pos, pos2, vel, vel2, partPropD.radius, partPropD.radius);
            }
        }
    }
    return force;
}


__global__
void collideD(float2* sortPos,               // input: sorted positions
              float2* sortVel,               // input: sorted velocities
              float2* newAcc,                // output: new acceleration
              uint*   gridParticleIndex,     // input: sorted particle indices
              uint*   cellStart,
              uint*   cellEnd)
{
	uint numParticles = simPropD.numParticles;
	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= numParticles) return;    
    
    // read particle data from sorted arrays
	float2 pos = sortPos[index];
    float2 vel = sortVel[index];

    // get address in grid
    int2 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float2 force = make_float2(0.0f);
    for(int y=-1; y<=1; y++) {
        for(int x=-1; x<=1; x++) {
            int2 neighbourPos = gridPos + make_int2(x, y);
            force += collideCell(neighbourPos, index, pos, vel, sortPos, sortVel, cellStart, cellEnd);
        }
    }

	newAcc[index] = force / partPropD.mass;
}

__global__
void integrateSystemD(float2* pos, float2* vel, float2* acc)
{
	uint numParticles = simPropD.numParticles;
	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= numParticles) return;  
    
    vel[index] += (simPropD.gravity + acc[index]) * simPropD.timeStep;
    pos[index] += vel[index] * simPropD.timeStep;
    
        // set this to zero to disable collisions with cube sides
#if 1
        if (pos[index].x > simPropD.cubeDimension.x - partPropD.radius) {
        	pos[index].x = simPropD.cubeDimension.x - partPropD.radius;
        	vel[index].x *= partPropD.boundaryDamping; }
        if (pos[index].x < partPropD.radius){
        	pos[index].x = partPropD.radius;
        	vel[index].x *= partPropD.boundaryDamping;}
        if (pos[index].y > simPropD.cubeDimension.x - partPropD.radius) { 
        	pos[index].y = simPropD.cubeDimension.x - partPropD.radius;
        	vel[index].y *= partPropD.boundaryDamping; }
#endif  
        if (pos[index].y < partPropD.radius) {
        	pos[index].y = partPropD.radius;
        	vel[index].y *= partPropD.boundaryDamping;}
}

__global__
void plotBackgroundD(uchar4 *ptr){
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
    uint index = x + y * blockDim.x * gridDim.x;
	
	ptr[index] = make_uchar4(0,0,0,0);
}

__global__
void plotSpheresD(uchar4*	ptr,
				  float2* 	sortPos)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	
	uint numParticles = simPropD.numParticles;
	
	if (index < numParticles) {
	
		float2 pos = sortPos[index];
	
		int cPixelx = simPropD.imageDIMx/simPropD.cubeDimension.x*pos.x;
		int cPixely = simPropD.imageDIMy/simPropD.cubeDimension.y*pos.y;
	
		for (int x = -simPropD.dimx/2; x < simPropD.dimx/2; x++ ) {
			for (int y = -simPropD.dimy/2; y < simPropD.dimy/2; y++) {
				if (x*x + y*y < simPropD.pRadius*simPropD.pRadius) {
			
					uint gPixelx = cPixelx + x;
					uint gPixely = cPixely + y;
				
					float fscale = sqrtf((simPropD.pRadius*simPropD.pRadius - x*x - y*y)/(simPropD.pRadius*simPropD.pRadius));
				
					uint pixel = gPixelx + gPixely*simPropD.imageDIMx;
				
					if (pixel >= simPropD.imageDIMx*simPropD.imageDIMy) pixel = simPropD.imageDIMx*simPropD.imageDIMy-1;
					
					ptr[pixel].x = 255.0f * fscale;
					ptr[pixel].y = 255.0f * fscale;
					ptr[pixel].z = 255.0f * fscale;
					ptr[pixel].w = 255.0f * fscale;
					
				}
			}
		}
	}
}
