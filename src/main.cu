// Input and Output include
#include <iostream>								  // entra e saída de dados
#include <time.h> 		   // biblioteca de tempo para criar o seed do rand

// CUDA includes
#include "gpu_anim.h" 				  // bib. de vizualização em tempo real
#include "cutil_math.h" 		      // funções matemáticas de vetores

// Dependece files
#include "main.cuh"
#include "functions.cuh" 	 // arquivo de funções de preparação para a GPU


#define log2( x ) log(x)/log(2)

void PrepareSim( SystemProperties *sisProps,
				 ParticlesValues *partValues,
				 ParticleProperties *partProps ) {

	sisProps->numParticles = PARTICLES;

	sisProps->cubeDimension.x = sisProps->cubeDimension.y = BOX_SIZE;
	
	sisProps->timeStep = TIME_STEP;
	
	sisProps->gravity = make_float2(0,-GRAVITY);
	
	sisProps->imageDIMx = DIM;
	sisProps->imageDIMy = DIM;
		
	partProps->radius = 16e-3f;
	partProps->mass = 1e-2;
	partProps->collideStiffness = 1e3;
	partProps->collideDamping = 0.2f;
	partProps->boundaryDamping = BOUNDARYDAMPING;
	
	// tamanho do quadrado que contem a esfera em PIXEL (para a saida grafica)
	sisProps->dimx = ceil(sisProps->imageDIMx/sisProps->cubeDimension.x*partProps->radius)*2;
	if (sisProps->dimx < 2) sisProps->dimx = 2;
	sisProps->dimy = ceil(sisProps->imageDIMy/sisProps->cubeDimension.y*partProps->radius)*2;
	if (sisProps->dimy < 2) sisProps->dimy = 2;
	
	// raio da esfera em PIXEL (para a saída grafica)
	sisProps->pRadius = sisProps->imageDIMy/sisProps->cubeDimension.y*partProps->radius;

	// Bloco inicial de esferas
	float corner1[2] = {0.1, 0.1}; 				 // canto inferior esquerdo
	float corner2[2] = {9.9, 9.9}; 				  // canto superior direito
	float sideLenght[2];
	sideLenght[0] = corner2[0] - corner1[0]; 			   // dimensao em X
	sideLenght[1] = corner2[1] - corner1[1]; 			   // dimensao em Y
	
	uint side[2] = {X_PARTICLES, Y_PARTICLES}; // numero de partículas em X
											  // e Y (deve ser maior que 2)

	// Calcula o tamanho do grid arredondando para um valor que seja
	// potencia de 2. O grid deve ser de 1.2 a 3 vezes o diametro da esfera
	uint grid = sisProps->cubeDimension.x / (4.0f * partProps[0].radius);
	uint temp = log2(grid);
	uint gridUpdate = 1 << temp;
	float cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;
	sisProps->gridSize.x = 1 << temp;
	
	grid = sisProps->cubeDimension.y / (4 * partProps[0].radius);
	temp = log2(grid);
	gridUpdate = 1 << temp;
	cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * partProps[0].radius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * partProps[0].radius ) temp += 1;	
	sisProps->gridSize.y = 1 << temp;

	sisProps->numCells = sisProps->gridSize.x * sisProps->gridSize.y;
	
	// alocando vetores na placa de video
	// cudaMalloc --> aloca espaço na placa de vídeo
	// cudaMemcpy --> transfere dados entre a CPU (Host) e GPU (Device)
	// cudaMemcpyToSymbol --> copia variável para a memória de constante
	float *d_corner1, *d_sideLenght;
	uint *d_side;

	cudaMalloc((void**)&d_corner1, sizeof(float)*2);
	cudaMalloc((void**)&d_sideLenght, sizeof(float)*2);
	cudaMalloc((void**)&d_side, sizeof(uint)*2);
	cudaMalloc((void**)&partValues->pos1, sizeof(float) * sisProps->numParticles * 2);
	cudaMalloc((void**)&partValues->pos2, sizeof(float) * sisProps->numParticles * 2);
	cudaMalloc((void**)&partValues->vel1, sizeof(float) * sisProps->numParticles * 2);
	cudaMalloc((void**)&partValues->vel2, sizeof(float) * sisProps->numParticles * 2);
	cudaMalloc((void**)&partValues->acc, sizeof(float) * sisProps->numParticles * 2);
	cudaMalloc((void**)&partValues->cellStart, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->cellEnd, sizeof(uint) * sisProps->numCells);
	cudaMalloc((void**)&partValues->gridParticleIndex, sizeof(uint) * sisProps->numParticles);
	cudaMalloc((void**)&partValues->gridParticleHash, sizeof(uint) * sisProps->numParticles);
	cudaMemcpy(d_corner1, corner1, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sideLenght, sideLenght, sizeof(float)*2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_side, side, sizeof(uint)*2, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(sisPropD, sisProps, sizeof(SystemProperties));
	cudaMemcpyToSymbol(partPropD, partProps , sizeof(ParticleProperties) * 1);

	// Função para definir a posição inicial das esferas
	initializeParticlePosition(partValues->pos1,
							   partValues->vel1,
							   partValues->acc,
							   d_corner1,
							   d_sideLenght,
							   d_side,
							   side,
							   time(NULL));
	
	// Desalocando espaço na placa de vídeo (Não mais necessário)
    cudaFree( d_corner1 );
    cudaFree( d_sideLenght );
    cudaFree( d_side );

	// Screen output	
	printf("\nNumero de Particulas = %d\n",sisProps->numParticles);
	printf("grid %d x %d\n\n",sisProps->gridSize.x,sisProps->gridSize.y);
#ifdef USE_TEX
	printf("Memoria de textura: UTILIZADA\n\n");
#else
	printf("Memoria de textura: NAO\n\n");
#endif 
}

void SimLooping( uchar4 *image, DataBlock *simBlock, int ticks ) {

	// inicia o cronometro
	simBlock->start = clock();

	// Estruturas auxiliares
    SystemProperties *sisProps = &simBlock->sisProps;
    ParticlesValues *partValues = &simBlock->partValues;
	
	// para ordenarmos os vetores de posicao e velocidade sem necessidade
	// de retornarmos a variável para o vetor original, um switch entre os
	// dois vetores de posição alocados na GPU é criado. A cada iteração o
	// vetor de início e o vetor reorganizado são invertidos, reduzindo uma
	// operação de cópia
	float *oldPos, *oldVel, *sortPos, *sortVel;
	
	// Integrando o programa IPS vezes antes de exibir a imagem
	for (int i = 0 ; i < simBlock->IPS ; i++) {

		if ((ticks + i) & 1) // quando par (FALSE) quando impar (TRUE)
		{	
			oldPos = partValues->pos1;
			sortPos = partValues->pos2;
			oldVel = partValues->vel1;
			sortVel = partValues->vel2;
		} else {
			oldPos = partValues->pos2;
			sortPos = partValues->pos1;
			oldVel = partValues->vel2;
			sortVel = partValues->vel1;
		}
		
		// Integracao no tempo (atualizacao das posicoes e velocidades)
		integrateSystem(oldPos,
			 	  		oldVel,
			 	  		partValues->acc,
			 	  		sisProps->numParticles);

		// Define a celula de cada particula, criando os vetores
		// gridParticleIndex e gridParticleHash ordenados pelo Index
		calcHash(oldPos,
				 partValues->gridParticleIndex,
				 partValues->gridParticleHash,
				 sisProps->numParticles);

		// Reordena os vetores baseado no Hash
		sortParticles(partValues->gridParticleHash,
					  partValues->gridParticleIndex,
					  sisProps->numParticles);

		// Reorganiza as variaveis de Pos e Vel para a nova ordem de particulas
		// e cria os vetores indicando a partícula de início e fim de cada
		// celula
		reorderDataAndFindCellStart(partValues->cellStart,
									partValues->cellEnd,
									sortPos,
									sortVel,
									partValues->gridParticleHash,
									partValues->gridParticleIndex,
									oldPos,
									oldVel,
									sisProps->numParticles,
									sisProps->numCells);

		// Detecta a colisao das particulas e transforma a força de colisão em
		// aceleração
		collide(sortPos,
				sortVel,
				partValues->acc,
				partValues->cellStart,
				partValues->cellEnd,
				sisProps->numParticles,
				sisProps->numCells);

//		// Integracao no tempo (atualizacao das posicoes e velocidades)
//		integrateSystem(sortPos,
//			 	  		sortVel,
//			 	  		partValues->acc,
//			 	  		sisProps->numParticles);

		simBlock->tempo++;
	}

	// Saida grarica quando necessario
	plotParticles(image,
				  sortPos,
				  sisProps->numParticles,
				  sisProps->imageDIMx,
				  sisProps->imageDIMy);

	
	// calcula o tempo de exibição do frame
	double time = ((double)clock() - simBlock->start)/CLOCKS_PER_SEC;
	if (time < 0.003f) time = 0.03f;
	
	// Define o número de Interações por segundo para exibir a imagem em 
	// FPS (definida no cabeçalho) frames por segundo.
	// Após a conta, transforma o número em impar para não calcular duas
	// duas vezes a mesma iteração (por causa do switch)
	simBlock->IPS = floor(1.0f/time/FPS*simBlock->IPS);
	simBlock->IPS = simBlock->IPS | 0x0001;

}

void FinalizingSim( DataBlock *simBlock) {

    // Limpe aqui o que tiver que ser limpo
    cudaFree( simBlock->partValues.pos1 );
    cudaFree( simBlock->partValues.pos2 );
    cudaFree( simBlock->partValues.vel1 );
    cudaFree( simBlock->partValues.vel2 );
    cudaFree( simBlock->partValues.acc );
    cudaFree( simBlock->partValues.cellStart );
    cudaFree( simBlock->partValues.cellEnd );
    cudaFree( simBlock->partValues.gridParticleIndex );
    cudaFree( simBlock->partValues.gridParticleHash );
    
   	printf("Integracoes por plot = %d\n\n",simBlock->IPS);
	double elapsedTime = ((double)clock() - simBlock->totalStart)/CLOCKS_PER_SEC;
	double simulationTime = simBlock->tempo * simBlock->sisProps.timeStep;
	
	printf("Duracao da simulacao = %4.2f s\n",elapsedTime);
	printf("Tempo de simulacao = %4.2f s\n\n",simulationTime);
	printf("Razao de tempo (Real/Simulado) = %3.3f\n\n",elapsedTime/simulationTime);
	
}


int main() {
	
	// declarando estrutura de dados principal
    DataBlock simBlock;
    
    // declarando as subestruturas (apenas por facilidade)
    SystemProperties *sisProps = &simBlock.sisProps;
    ParticleProperties *partProps = &simBlock.partProps;
    ParticlesValues *partValues = &simBlock.partValues;
    
    // Definindo que a primeira iteração será exibida
    simBlock.IPS = 1;
	simBlock.totalStart = clock();
	simBlock.tempo = 0;
    
    // função que define o tamanho da imagem e a estrutura que será
    // repassada para dentro do looping
    GPUAnimBitmap bitmap(DIM, DIM, &simBlock );

	// Utilizar ARGC e ARGV para pegar propriedades na linha de comando
	// ler esses comandos de um arquivo TXT externo
	// Criar uma rotina para fazer este tipo de leitura
	
	// Prepara a simulacao, define as condicoes iniciais do problema
	PrepareSim(sisProps, partValues, partProps);
	
	// Executa o looping até que a tecla ESC seja pressionada
    bitmap.anim_and_exit(
        (void (*)(uchar4*,void*,int))SimLooping, (void (*)(void*))FinalizingSim );

}
