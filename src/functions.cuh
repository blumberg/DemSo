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

void allocateVectors(ParticleProperties*,
					 ParticlesValues*,
					 SystemProperties*,
					 RenderParameters*);

void desAllocateVectors(ParticlesValues*);

uint iDivUp(uint, uint);

void computeGridSize(uint, uint, uint &, uint &);

void initializeParticlePosition (float*, float*, float*, float*, float*, float*,
								 uint*, uint*, float*, float*, uint*, unsigned long, int);

void calcHash(float*, uint*, uint*, uint);

void sortParticles(uint*, uint*, uint);

void reorderDataAndFindCellStart(uint*, uint*, float*, float*, float*, float*,
								 uint*, uint*, uint*, uint*, float*, float*,
                                 float*, float*, uint*, uint*, uint, uint);

void collide(float*, float*, float*, float*, float*, uint*, uint*, uint*, uint, uint);

void integrateSystem(float*, float*, float*, float*, float*, float*, uint*, uint);

void plotParticles(uchar4*, float*, float*, uint*, uint, int, int);

void writeOutputFile (DataBlock *, int);

#endif /* FUNCTIONS_CUH */
