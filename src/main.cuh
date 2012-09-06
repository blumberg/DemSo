
// Estrutura de propriedades das partículas, caso houver mais do que um
// tipo de partícula, essa estrutura será criada com um tamanho maior
// Essa estrutura será carregana na memória de constantes da GPU
struct ParticleProperties {

	float radius;
	float mass;
	float collideStiffness;
	float collideDamping;
	float boundaryDamping;

};

// Estrutura com os valores armazenados de cada partícula. Todas as
// variáveis dessa estrutura serão alocadas na GPU.
struct ParticlesValues {

	float *pos1, *vel1, *acc;
	float *pos2, *vel2;
	uint *cellStart, *cellEnd;
	uint *gridParticleIndex, *gridParticleHash;

};

// Estrutura com as propriesdades do sistema. Seu tamanho será fixo e esta
// estrutura inteira será passada para a memória de constantes da GPU
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

// Estrutura principal da simulação. Ela contem as 3 outras subestruturas.
struct DataBlock {

	ParticleProperties partProps;
	ParticlesValues partValues;
	SistemProperties sisProps;
	
	clock_t start, totalStart;
	int tempo;
	int IPS;

};
