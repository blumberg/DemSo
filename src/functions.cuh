#ifndef FUNCTIONS_CUH
#define FUNCTIONS_CUH

void allocateVectors(ParticleProperties*,
					 ParticlesValues*,
					 SystemProperties*,
					 RenderParameters*);

void desAllocateVectors(ParticlesValues*);

uint iDivUp(uint, uint);

void computeGridSize(uint, uint, uint &, uint &);

void initializeParticlePosition (float*, float*, float*, uint*, uint*,
								 float*, float*, uint*, unsigned long);

void calcHash(float*, uint*, uint*, uint);

void sortParticles(uint*, uint*, uint);

void reorderDataAndFindCellStart(uint*, uint*, float*, float*, uint*,
                                 uint*, uint*, uint*, float*, float*,
                                 uint*, uint*, uint, uint);

void collide(float*, float*, float*, uint*, uint*, uint*, uint, uint);

void integrateSystem(float*, float*, float*, uint*, uint);

void plotParticles(uchar4*, float*, uint*, uint, int, int);

#endif /* FUNCTIONS_CUH */
