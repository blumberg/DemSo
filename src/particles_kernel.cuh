#ifndef PARTICLES_KERNEL_CUH
#define PARTICLES_KERNEL_CUH

// Define se pega a variável da memória global ou da memória de textura
#if USE_TEX
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif

#if USE_TEX
// textures for particle position and velocity
extern texture<float2, 1, cudaReadModeElementType> oldPosTex;
extern texture<float2, 1, cudaReadModeElementType> oldVelTex;

//texture<uint, 1, cudaReadModeElementType> gridParticleHashTex;
extern texture<uint, 1, cudaReadModeElementType> cellStartTex;
extern texture<uint, 1, cudaReadModeElementType> cellEndTex;
#endif

// Declarando as variáveis da memória de constante
__constant__ SystemProperties sisPropD;
__constant__ ParticleProperties partPropD;

__global__ void initializeParticlePositionD(float2*, float2*, float2*,
											float*, float*, uint*,
											unsigned long);

__device__ int2 calcGridPos(float2);

__device__ uint calcGridHash(int2);

__global__ void calcHashD(uint*, uint*, float2*);

__global__ void reorderDataAndFindCellStartD(uint*, uint*, float2*,
											float2*, uint*, uint*, float2*, float2*);

__device__ float2 collideSpheres(float2, float2, float2, float2, float, float);

__device__ float2 collideCell(int2, uint, float2, float2, float2*, float2*, uint*, uint*);

__global__ void collideD(float2*, float2*, float2*, uint*, uint*);

__global__ void integrateSystemD(float2*, float2*, float2*);

__global__ void plotSpheresD(uchar4*, float2*);

#endif /* PARTICLES_KERNEL_CUH */
