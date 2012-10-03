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
void createRetangleBlockD(float2*			pos,
						  uint*				ID,
						  uint*				loc,
						  uint*				type,
						  float2			start,
						  float2			sideLenght,
						  uint2				side,
						  uint				startID,
						  uint 				numParticleTypes,
						  uint*				particleTypeVec,
						  unsigned long 	seed) {
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x >= side.x) return;
	if (y >= side.y) return;
	
	uint particle = x + y*side.x;

    if (particle >= side.y*side.x) return;

	curandState state;
	curand_init( seed, particle, 0, &state );

	pos[particle].x = start.x + sideLenght.x/(side.x-1) * (x + (curand_normal(&state)-0.5f)/100);
	pos[particle].y = start.y + sideLenght.y/(side.y-1) * (y + (curand_normal(&state)-0.5f)/100);
	ID[particle] = startID + particle;
	loc[particle] = startID + particle;
	int vecPos = (x+y) % numParticleTypes;
	type[particle] = particleTypeVec[vecPos];
}

__global__
void createTriangleBlockD (float2*	pos,
						   uint*	ID,
						   uint*	loc,
						   uint*	type,
						   float2	start,
						   uint		N,
						   uint 	numParticleTypes,
						   uint*	particleTypeVec,
						   float	space,
						   float 	height,
						   uint 	startID,
						   uint 	numParticles){
	
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (index >= numParticles) return;
	
	int L1 = (1 + 2*sqrtf(N*N + N -2*(index+1) + 1)) / 2;
	int L = N - L1;
	
	int s = ( L1 + N + N*N - L1*L1 ) / 2;
	
	int C = index - s + L1;
	
	int A = ( N - 1 );
	
	pos[index].x = start.x + (L + 2*C - A)*space/2;
	pos[index].y = start.y + L*height;
	
	ID[index] = startID + index;
	loc[index] = startID + index;
	int vecPos = (C + L) % numParticleTypes;
	type[index] = particleTypeVec[vecPos];

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

    // get address in grid
    int2 gridPos = calcGridPos(p);
    uint hash = calcGridHash(gridPos);

    // store grid hash and particle index
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
								  float*	sortedTheta,	  // output: sorted angular positions
								  float*	sortedOmega,	  // output: sorted angular velocities
  							      uint*		sortedID,		  // output: sorted Identifications
  							      uint*		sortedLoc,		  // output: sorted Localization
  							      uint*		sortedType,		  // output: sorted Type
                                  uint*  	gridParticleHash, // input: sorted grid hashes
                                  uint*  	gridParticleIndex,// input: sorted particle indices
				                  float2* 	oldPos,           // input: sorted position array
							      float2* 	oldVel,           // input: sorted velocity array
								  float*	oldTheta,		  // input: sorted angular position array
								  float*	oldOmega,		  // input: sorted angular velocity array
							      uint*		oldID,			  // input: sorted Identifications array
							      uint*		oldType)		  // input: sorted Type array
{
	extern __shared__ uint sharedHash[];    // blockSize + 1 elements
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	
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
		float theta = FETCH(oldTheta, sortedIndex);
		float omega = FETCH(oldOmega, sortedIndex);
        uint ID = FETCH(oldID, sortedIndex);
        uint type = FETCH(oldType, sortedIndex);

        sortedPos[index] = pos;
        sortedVel[index] = vel;
		sortedTheta[index] = theta;
		sortedOmega[index] = omega;
        sortedID[index] = ID;
        sortedLoc[ID] = index;
        sortedType[index] = type;
	}
}

// collide two spheres using DEM method
__device__
float2 collideSpheres(float2 posA, float2 posB,
                      float2 velA, float2 velB,
					  float omegaA, float omegaB,
                      uint typeA, uint typeB, float &moment)
{
	// Getting radius
	float radiusA = partPropD[typeA].radius;
	float radiusB = partPropD[typeB].radius;
	
	// calculate relative position
    float2 relPos = posB - posA;

    float dist = length(relPos);
    float collideDist = radiusA + radiusB;

    float2 force = make_float2(0.0f);
    if (dist < collideDist)
	{
		// Normal contact vector
        float2 norm = relPos / dist;

		// Tangential contact vector
		float2 tang;
		tang.x = -norm.y;
		tang.y = norm.x;

		// relative velocity
        float2 relVel = velB - velA;
		float  relVel_n = dot(relVel, norm);
        float2 relVel_t = relVel - relVel_n*norm;

		// Series association of normal damping and stiffness
		float normalStiffness = (partPropD[typeA].normalStiffness*partPropD[typeB].normalStiffness)
							   /(partPropD[typeA].normalStiffness+partPropD[typeB].normalStiffness);
		float shearStiffness = (partPropD[typeA].shearStiffness*partPropD[typeB].shearStiffness)
							  /(partPropD[typeA].shearStiffness+partPropD[typeB].shearStiffness);
		float normalDamping = (partPropD[typeA].normalDamping*partPropD[typeB].normalDamping)
							 /(partPropD[typeA].normalDamping+partPropD[typeB].normalDamping);
		float frictionCoefficient = (partPropD[typeA].boundaryDamping+partPropD[typeB].boundaryDamping)/2;

        // spring force
        force = -normalStiffness*(collideDist - dist) * norm;

        // dashpot (damping) force (not present when particles are moving away from each-other)
        if (relVel_n < 0.0f) force += normalDamping * relVel_n*norm;

        // tangential shear force
		float2 contactVel_t = relVel_t - (radiusA*omegaA + radiusB*omegaB) * tang;
		float2 Ftrial = shearStiffness * contactVel_t * sisPropD.timeStep;
	
		// Max tangential friction force
		float Ftmax = frictionCoefficient*length(force);

		float2 Ft = make_float2(0.0f);
		if (length(Ftrial) < Ftmax) Ft = Ftrial;
		else Ft = Ftmax * contactVel_t / length(contactVel_t);

		// Shear force
		force += Ft;
		// Moment
		moment += radiusA * dot(Ft, tang); 
    }

    return force;
}



// collide a particle against all other particles in a given cell
__device__
float2 collideCell(int2		gridPos,
                   uint		index,
                   float2	pos,
                   float2	vel,
				   float	omega,
                   uint		type,
                   float2*	oldPos,
                   float2*	oldVel,
				   float*	oldOmega,
                   uint*	oldType,
                   uint*	cellStart,
                   uint*	cellEnd,
				   float	&moment)
{
    uint gridHash = calcGridHash(gridPos);

    // get start of bucket for this cell
    uint startIndex = FETCH(cellStart, gridHash);

    float2 force = make_float2(0.0f);
    if (startIndex != (uint)-1) {        // cell is not empty
        // iterate over particles in this cell
        uint endIndex = FETCH(cellEnd, gridHash);
        for(uint j=startIndex; j<endIndex; j++) {
            if (j != index) {              // check not colliding with self
	            float2 pos2 = FETCH(oldPos, j);
                float2 vel2 = FETCH(oldVel, j);
				float  omega2 = FETCH(oldOmega, j);
                uint type2 = FETCH(oldType, j);

                // collide two spheres
                force += collideSpheres(pos, pos2, vel, vel2, omega, omega2, type, type2, moment);
            }
        }
    }
    return force;
}


__device__
float2 collideBoundary(float2 &pos, float2 &vel, float omega,
                       uint type, float2 boundPos, float &moment)
{
	// Getting radius
	float radius = partPropD[type].radius;
	
	// calculate relative position
    float2 relPos = boundPos - pos;

    float dist = length(relPos);

    float2 force = make_float2(0.0f);
    if (dist < radius)
	{
		// Normal contact vector
        float2 norm = relPos / dist;

		// Tangential contact vector
		float2 tang;
		tang.x = -norm.y;
		tang.y = norm.x;

		// relative velocity
        float2 relVel = -vel;
		float  relVel_n = dot(relVel, norm);
        float2 relVel_t = relVel - relVel_n*norm;

		// Series association of normal damping and stiffness
		float normalStiffness = (partPropD[type].normalStiffness*sisPropD.boundaryNormalStiffness)
							   /(partPropD[type].normalStiffness+sisPropD.boundaryNormalStiffness);
		float shearStiffness = (partPropD[type].shearStiffness*sisPropD.boundaryShearStiffness)
							  /(partPropD[type].shearStiffness+sisPropD.boundaryShearStiffness);
		float normalDamping = (partPropD[type].normalDamping*partPropD[type].boundaryDamping)
							 /(partPropD[type].normalDamping+partPropD[type].boundaryDamping);

        // spring force
        force = -normalStiffness*(radius - dist) * norm;

        // dashpot (damping) force (not present when particles are moving away from each-other)
        if (relVel_n < 0.0f) force += normalDamping * relVel_n*norm;

        // tangential shear force
		float2 contactVel_t = relVel_t - radius*omega*tang;
		float2 Ftrial = shearStiffness * contactVel_t * sisPropD.timeStep;
	
		// Max tangential friction force
		float Ftmax = partPropD[type].frictionCoefficient*length(force);

		float2 Ft = make_float2(0.0f);
		if (length(Ftrial) < Ftmax) Ft = Ftrial;
		else Ft = Ftmax * contactVel_t / length(contactVel_t);

		// Shear force
		force += Ft;
		// Moment
		moment += radius * dot(Ft, tang); 

		// Fixing position and velocity
		if (pos.x >= sisPropD.cubeDimension.x || pos.x <= 0.0f) vel.y *= -1;
		if (pos.y >= sisPropD.cubeDimension.y || pos.y <= 0.0f) vel.x *= -1;
		pos = boundPos;
    }

    return force;
}


__global__
void collideD(float2*	oldPos,               // input: sorted positions
              float2*	oldVel,               // input: sorted velocities
              float2*	newAcc,                // output: new acceleration
			  float*  oldOmega,				// input: sorted angular velocities
			  float*  newAlpha,				// output: new angular acceleration
              uint*		oldType,
              uint*		cellStart,
              uint*		cellEnd
#if USE_BIG_PARTICLE
			  , float2	controlPos,
			  uint		controlType
#endif
			  )

{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;    
    
    // read particle data from sorted arrays
	float2 pos = FETCH(oldPos, index);
    float2 vel = FETCH(oldVel, index);
	float omega = FETCH(oldOmega, index);
    uint type = FETCH(oldType, index);

    // get address in grid
    int2 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float2 force = make_float2(0.0f);
	float moment = 0.0f;
    for(int y=-1; y<=1; y++) {
        for(int x=-1; x<=1; x++) {
            int2 neighbourPos = gridPos + make_int2(x, y);
            force += collideCell(neighbourPos, index, pos, vel, omega,
								 type, oldPos, oldVel, oldOmega,
								 oldType, cellStart, cellEnd, moment);
        }
    }

	// Check if cell is next to boundary
//	if (gridPos.x <= 0)
//		force += collideBoundary(pos, vel, omega, type,
//								 make_float2(0.0f, pos.y), moment);
//	else if(gridPos.x >= sisPropD.gridSize.x-1)
//		force += collideBoundary(pos, vel, omega, type,
//								 make_float2(sisPropD.cubeDimension.x, pos.y), moment);
//
//	if (gridPos.y <= 0)
//		force += collideBoundary(pos, vel, omega, type,
//								 make_float2(pos.x, 0.0f), moment);
//	else if(gridPos.y >= sisPropD.gridSize.y-1)
//		force += collideBoundary(pos, vel, omega, type,
//								 make_float2(pos.x, sisPropD.cubeDimension.y), moment);

    
#if USE_BIG_PARTICLE
    force += collideSpheres(pos, controlPos, vel, make_float2(0), omega, 0, type, controlType, moment);
#endif

	newAcc[index] = force / partPropD[type].mass;
	// moment / momentOfInertia
	newAlpha[index] = moment / partPropD[type].inertia;
}

__global__
void integrateSystemD(float2* pos, float2* vel, float2* acc,
					  float* theta, float* omega, float* alpha, uint* type)
{	
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= sisPropD.numParticles) return;  
    
    vel[index] += (sisPropD.gravity + acc[index]) * sisPropD.timeStep;
    pos[index] += vel[index] * sisPropD.timeStep;

	omega[index] += alpha[index] * sisPropD.timeStep;
	theta[index] += omega[index] * sisPropD.timeStep;

	// Correct the value of angular position
	if (theta[index] > 2*M_PI) theta[index] -= 2*M_PI*floor(theta[index]/(2*M_PI));
	if (theta[index] < 0) theta[index] -= 2*M_PI*floor(theta[index]/(2*M_PI));

	//TODO: implementar o contato com a borda aqui. É o único jeito de
	//garantir que as partículas não vão sair do universo...
	//
	//Seria como um fator de "correcão". Se a partícula sair, faz ela
	//colidir com a borda, calcula as forcas/momentos, atualiza as
	//aceleracões mas reverte as velocidades e trava a posicão
	float radius = partPropD[type[index]].radius;
#if 1
	if (pos[index].x > sisPropD.cubeDimension.x - radius) {
		pos[index].x = sisPropD.cubeDimension.x - radius;
		vel[index].x *= -partPropD[type[index]].boundaryDamping; }
	if (pos[index].x < radius) {
		pos[index].x = radius;
		vel[index].x *= -partPropD[type[index]].boundaryDamping; }
	if (pos[index].y > sisPropD.cubeDimension.y - radius) {
		pos[index].y = sisPropD.cubeDimension.y - radius;
		vel[index].y *= -partPropD[type[index]].boundaryDamping; }
	if (pos[index].y < radius) {
		pos[index].y = radius;
		vel[index].y *= -partPropD[type[index]].boundaryDamping; }
#endif
}


// Pinta a esfera.
__global__
void plotSpheresD(uchar4*	ptr,
				  float2* 	sortPos,
				  float*	sortTheta,
				  uint*		type)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (index >= sisPropD.numParticles) return;
	
	float2 pos = sortPos[index];
	float theta = sortTheta[index];
	
	uint currentType = type[index];
	float pRadius = renderParD.imageDIMy/sisPropD.cubeDimension.y*partPropD[currentType].radius;
	
	// calcula a posição do centro da partícula em pixel
	int cPixelx = renderParD.imageDIMx/sisPropD.cubeDimension.x*pos.x;
	int cPixely = renderParD.imageDIMy/sisPropD.cubeDimension.y*pos.y;
	
	// percorre o quadrado ocupado pela partícula (em pixel)
//	for (int x = -renderParD.dimx/2; x < renderParD.dimx/2; x++ ) {
	for (int x = -ceil(pRadius); x < ceil(pRadius); x++ ) {
//		for (int y = -renderParD.dimy/2; y < renderParD.dimy/2; y++) {
		for (int y = -ceil(pRadius); y < ceil(pRadius); y++) {
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
	// Drawing line to view rotations
	for(register int r = 0; r < pRadius; r++)
	{
		// posição do ponto atual (em pixel)
		uint gPixelx = cPixelx + r * cosf(theta);
		uint gPixely = cPixely + r * sinf(theta);
		
		// posição do pixel no vetor da imagem
		uint pixel = gPixelx + gPixely*renderParD.imageDIMx;
		if (pixel >= renderParD.imageDIMx*renderParD.imageDIMy) pixel = renderParD.imageDIMx*renderParD.imageDIMy-1;

		// define a cor do pixel
		ptr[pixel].x = 0.0f;
		ptr[pixel].y = 0.0f;
		ptr[pixel].z = 0.0f;
		ptr[pixel].w = 0.0f;
	}
	
}

__global__
void plotControlParticleD(uchar4*	ptr,
				  		  float2 	pos,
				  		  uint		type)
{
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
	
	int halfSidex = blockDim.x*gridDim.x/2.0f;
	int halfSidey = blockDim.y*gridDim.y/2.0f;
	
	float rx = x - halfSidex;
	float ry = y - halfSidey;
	
	float pRadius = renderParD.imageDIMy/sisPropD.cubeDimension.y*partPropD[type].radius;
	
	if ((rx*rx + ry*ry) < (halfSidex*halfSidey)){
		
		int cPixelx = renderParD.imageDIMx/sisPropD.cubeDimension.x*pos.x;
		int cPixely = renderParD.imageDIMy/sisPropD.cubeDimension.y*pos.y;
	
		uint gPixelx = cPixelx + x - blockDim.x*gridDim.x/2;
		uint gPixely = cPixely + y - blockDim.y*gridDim.y/2;

		if (gPixelx >= renderParD.imageDIMx) return;
		if (gPixely >= renderParD.imageDIMy) return;

		float fscale = sqrtf((pRadius*pRadius - rx*rx - ry*ry)/(pRadius*pRadius));
		
		uint pixel = gPixelx + gPixely*renderParD.imageDIMx;
		if (pixel >= renderParD.imageDIMx*renderParD.imageDIMy) pixel = renderParD.imageDIMx*renderParD.imageDIMy-1;
		
		// define a cor do pixel
		ptr[pixel].x = partPropD[type].colorR * fscale;
		ptr[pixel].y = partPropD[type].colorG * fscale;
		ptr[pixel].z = partPropD[type].colorB * fscale;
		ptr[pixel].w = 255.0f * fscale;
	
	}
}
#endif
