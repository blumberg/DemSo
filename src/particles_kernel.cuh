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

#ifndef PARTICLES_KERNEL_CUH
#define PARTICLES_KERNEL_CUH

#include <curand_kernel.h>  		  // bib. randomica para kernel em CUDA
#include "cutil_math.h"       // funções matemáticas de vetores
#include "main.cuh"

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
texture<uint, 1, cudaReadModeElementType> oldIDTex;
texture<uint, 1, cudaReadModeElementType> oldTypeTex;

texture<uint, 1, cudaReadModeElementType> cellStartTex;
texture<uint, 1, cudaReadModeElementType> cellEndTex;
#endif

// Declarando as variáveis da memória de constante
__constant__ SystemProperties sisPropD;
__constant__ ParticleProperties partPropD[MAX_PARTICLES_TYPES];
__constant__ RenderParameters renderParD;

__global__
void initializeParticlePositionD(float2*			pos,
								 float2*			vel,
								 float2*			acc,
								 uint*				ID,
								 uint*				type,
								 float*				corner1,
								 float*				comp,
								 uint*				side,
								 unsigned long 		seed,
								 int 				numParticleTypes) {
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x >= side[0]) return;
	if (y >= side[1]) return;
	
	uint particle = x + y*side[0];

    if (particle >= sisPropD.numParticles) return;

	curandState state;
	curand_init( seed, particle, 0, &state );

	pos[particle].x = corner1[0] + comp[0]/(side[0]-1) * (x + (curand_normal(&state)-0.5f)/100);
	pos[particle].y = corner1[1] + comp[1]/(side[1]-1) * (y + (curand_normal(&state)-0.5f)/100);
	vel[particle] = make_float2( 0 );
	acc[particle] = make_float2( 0 );
	ID[particle] = particle;
	type[particle] = (particle+y) % numParticleTypes;
//	type[particle] = 0;
}

// calculate position in uniform grid
__device__
int2 calcGridPos(float2 p)
{
    int2 gridPos;
    gridPos.x = floor(p.x * sisPropD.gridSize.x / sisPropD.cubeDimension.x);
    gridPos.y = floor(p.y * sisPropD.gridSize.y / sisPropD.cubeDimension.y);
    return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__
uint calcGridHash(int2 gridPos)
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
void reorderDataAndFindCellStartD(uint*   	cellStart,        // output: cell start index
							      uint*   	cellEnd,          // output: cell end index
							      float2*	sortedPos,        // output: sorted positions
  							      float2* 	sortedVel,        // output: sorted velocities
  							      uint*		sortedID,		  // output: sorted Identifications
  							      uint*		sortedType,		  // output: sorted Type
                                  uint*  	gridParticleHash, // input: sorted grid hashes
                                  uint*  	gridParticleIndex,// input: sorted particle indices
				                  float2* 	oldPos,           // input: sorted position array
							      float2* 	oldVel,           // input: sorted velocity array
							      uint*		oldID,			  // input: sorted Identifications array
							      uint*		oldType)		  // input: sorted Type array
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
        uint ID = FETCH(oldID, sortedIndex);
        uint type = FETCH(oldType, sortedIndex);

        sortedPos[index] = pos;
        sortedVel[index] = vel;
        sortedID[index] = ID;
        sortedType[index] = type;
	}
}

// collide two spheres using DEM method
__device__
float2 collideSpheres(float2 posA, float2 posB,
                      float2 velA, float2 velB,
                      uint typeA, uint typeB)
{
	// Getting radius
	float radiusA = partPropD[typeA].radius;
	float radiusB = partPropD[typeB].radius;
	
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

		float collideStiffness = (partPropD[typeA].collideStiffness*partPropD[typeB].collideStiffness)/(partPropD[typeA].collideStiffness+partPropD[typeB].collideStiffness);
		float collideDamping = (partPropD[typeA].collideDamping*partPropD[typeB].collideDamping)/(partPropD[typeA].collideDamping+partPropD[typeB].collideDamping);

        // spring force
        force = -collideStiffness*(collideDist - dist) * norm;
        // dashpot (damping) force
        force += collideDamping*relVel;
        // tangential shear force
//        force += params.shear*tanVel;
    }

    return force;
}



// collide a particle against all other particles in a given cell
__device__
float2 collideCell(int2		gridPos,
                   uint		index,
                   float2	pos,
                   float2	vel,
                   uint		type,
                   float2*	oldPos,
                   float2*	oldVel,
                   uint*	oldType,
                   uint*	cellStart,
                   uint*	cellEnd)
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
                uint type2 = FETCH(oldType, j);

                // collide two spheres
                force += collideSpheres(pos, pos2, vel, vel2, type, type2);
            }
        }
    }
    return force;
}


__global__
void collideD(float2* oldPos,               // input: sorted positions
              float2* oldVel,               // input: sorted velocities
              float2* newAcc,                // output: new acceleration
              uint*	  oldType,
              uint*   cellStart,
              uint*   cellEnd)
{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;    
    
    // read particle data from sorted arrays
	float2 pos = FETCH(oldPos, index);
    float2 vel = FETCH(oldVel, index);
    uint type = FETCH(oldType, index);

    // get address in grid
    int2 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float2 force = make_float2(0.0f);
    for(int y=-1; y<=1; y++) {
        for(int x=-1; x<=1; x++) {
            int2 neighbourPos = gridPos + make_int2(x, y);
            force += collideCell(neighbourPos, index, pos, vel, type, oldPos, oldVel, oldType, cellStart, cellEnd);
        }
    }

	newAcc[index] = force / partPropD[type].mass;
}

__global__
void integrateSystemD(float2* pos, float2* vel, float2* acc, uint* type)
{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;  
    
    vel[index] += (sisPropD.gravity + acc[index]) * sisPropD.timeStep;
    pos[index] += vel[index] * sisPropD.timeStep;
    
    float radius = partPropD[type[index]].radius;
    float boundaryDamping = partPropD[type[index]].boundaryDamping;
    
        // set this to zero to disable collisions with cube sides
#if 1
        if (pos[index].x > sisPropD.cubeDimension.x - radius) {
        	pos[index].x = sisPropD.cubeDimension.x - radius;
        	vel[index].x *= boundaryDamping; }
        if (pos[index].x < radius){
        	pos[index].x = radius;
        	vel[index].x *= boundaryDamping;}
        if (pos[index].y > sisPropD.cubeDimension.x - radius) { 
        	pos[index].y = sisPropD.cubeDimension.x - radius;
        	vel[index].y *= boundaryDamping; }
#endif  
        if (pos[index].y < radius) {
        	pos[index].y = radius;
        	vel[index].y *= boundaryDamping;}
}


// Pinta a esfera.
__global__
void plotSpheresD(uchar4*	ptr,
				  float2* 	sortPos,
				  uint*		type)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (index >= sisPropD.numParticles) return;
	
	float2 pos = sortPos[index];
	
	uint currentType = type[index];
	float pRadius = renderParD.imageDIMy/sisPropD.cubeDimension.y*partPropD[currentType].radius;
	
	// calcula a posição do centro da partícula em pixel
	int cPixelx = renderParD.imageDIMx/sisPropD.cubeDimension.x*pos.x;
	int cPixely = renderParD.imageDIMy/sisPropD.cubeDimension.y*pos.y;
	
	// percorre o quadrado ocupado pela partícula (em pixel)
	for (int x = -renderParD.dimx/2; x < renderParD.dimx/2; x++ ) {
		for (int y = -renderParD.dimy/2; y < renderParD.dimy/2; y++) {
			if (x*x + y*y < pRadius*pRadius) {
				
				// posição do ponto atual (em pixel)
				uint gPixelx = cPixelx + x;
				uint gPixely = cPixely + y;
				
				// Cria o efeito 3D da partícula (escurece as bordas)
				float fscale = sqrtf((pRadius*pRadius - x*x - y*y)/(pRadius*pRadius));
				
				// posição do pixel no vetor da imagem
				uint pixel = gPixelx + gPixely*renderParD.imageDIMx;
				if (pixel >= renderParD.imageDIMx*renderParD.imageDIMy) pixel = renderParD.imageDIMx*renderParD.imageDIMy-1;
				
				// define a cor do pixel
				ptr[pixel].x = partPropD[currentType].colorR * fscale;
				ptr[pixel].y = partPropD[currentType].colorG * fscale;
				ptr[pixel].z = partPropD[currentType].colorB * fscale;
				ptr[pixel].w = 255.0f * fscale;
				
			}
		}
	}
}
#endif
