
struct ParticleProperties {

	float radius;
	float mass;
	float collideStiffness;
	float collideDamping;
	float boundaryDamping;

};

struct ParticlesValues {

	float2 *pos1, *vel1, *acc;
	float2 *pos2, *vel2;
	uint *cellStart, *cellEnd;
	uint *gridParticleIndex, *gridParticleHash;

};

struct SistemProperties {

    uint numParticles;

    float2 cubeDimension;
    uint2 gridSize;
	uint numCells;
    
    float timeStep;
    
    float2 gravity;
    
    int imageDIMx;
    int imageDIMy;
	int dimx;
	int dimy;
	float pRadius;

};

struct DataBlock {

	ParticleProperties partProps;
	ParticlesValues partValues;
	SistemProperties sisProps;

};
