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
#include <iostream>

std::ostream& operator<<(std::ostream&, const uint2&);
std::ostream& operator<<(std::ostream&, const float2&);
std::ostream& operator<<(std::ostream&, const float3&);

void allocateVectors(ParticleProperties*,
					 ParticlesValues*,
					 SystemProperties*,
					 RenderParameters*);

void desAllocateVectors(ParticlesValues*);

uint iDivUp(uint, uint);

void computeGridSize(uint, uint, uint &, uint &);

void createRectangles (float*, uint*, uint*, uint*, float2, float2, 
					   uint2, uint, uint, uint*, unsigned long);

void createTriangles (float*, uint*, uint*, uint*, float2, uint, 
					  uint, uint*, float, float, uint, uint);

void createSingleParticles (float*, float*, float*, float*, uint*, 
							uint*, uint*, float2*, float2*, float*, 
							float*, uint*, uint, uint);

void calcHash(float*, uint*, uint*, uint);

void sortParticles(uint*, uint*, uint);

void reorderDataAndFindCellStart(uint*, uint*, float*, float*, float*,
                                 float*, uint*, uint*, uint*, uint*, uint*,
                                 float*, float*, float*, float*, uint*,
                                 uint*, uint, uint);

void collide(float*, float*, float*, float*, float*, uint*, uint*, uint*, uint, uint,
#if USE_BIG_PARTICLE
	float2, float2, float, float, uint,
#if USE_ATOMIC
	float2*, float*,
#else
	float*, float*, float*, float*, float*, float*,
#endif // USE_ATOMIC
	float2*, float*,
#endif // USE_BIG_PARTICLE
	float*);


void integrateSystem(float*, float*, float*, float*, float*, float*, uint*, float*, uint);

void plotParticles(uchar4*, float*, float*, uint*, uint, int, int, int,
#if USE_BIG_PARTICLE
				  float2, uint, int, int,
#endif
				  float*);

void writeOutputFile (FILE*, std::vector<int>, float, float2*, float2*, float2*, float*, float*,
					  float*, uint*, uint*, uint*, float*);

void set_gravity (SystemProperties *, float2);

void set_viewRotations (RenderParameters *, bool);

void set_colorByPressure (RenderParameters *, bool);

void updatePressureScale (RenderParameters *, float *, int);

#endif /* FUNCTIONS_CUH */
