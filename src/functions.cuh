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

#ifndef FUNCTIONS_CUH
#define FUNCTIONS_CUH
#include <vector>

void allocateVectors(ParticleProperties*,
					 ParticlesValues*,
					 SystemProperties*,
					 RenderParameters*);

void desAllocateVectors(ParticlesValues*);

uint iDivUp(uint, uint);

void computeGridSize(uint, uint, uint &, uint &);

void createRetangleBlock (float*, uint*, uint*, uint*, float2, float2, 
						  uint2, uint, uint, uint*, unsigned long);

void createTriangleBlock (float*, uint*, uint*, uint*, float2, uint, 
						  uint, uint*, float, float, uint, uint);

void createUserDefineBlock (float*, float*, float*, float*, uint*, 
							uint*, uint*, float2*, float2*, float*, 
							float*, uint*, uint, uint);

void calcHash(float*, uint*, uint*, uint);

void sortParticles(uint*, uint*, uint);

void reorderDataAndFindCellStart(uint*, uint*, float*, float*, float*,
                                 float*, uint*, uint*, uint*, uint*, uint*,
                                 float*, float*, float*, float*, uint*,
                                 uint*, uint, uint);

void collide(float*, float*, float*, float*, float*, uint*, uint*, uint*, uint, uint
#if USE_BIG_PARTICLE
	, float2, float2, float, float, uint, float2*, float*
#endif
	);

void integrateSystem(float*, float*, float*, float*, float*, float*, uint*, uint);

void plotParticles(uchar4*, float*, float*, uint*, uint, int, int
#if USE_BIG_PARTICLE
				  , float2, uint, int, int
#endif
				  );

void writeOutputFile (FILE*, std::vector<int>, float, float2*, float2*, float2*, float*, float*,
					  float*, uint*, uint*, uint*);

#endif /* FUNCTIONS_CUH */
