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

// Input and Output include
#include <iostream>								  // entra e saída de dados
#include <ctime> 		   // biblioteca de tempo para criar o seed do rand

// CUDA includes
#include "gpu_anim.h" 				  // bib. de vizualização em tempo real
#include "cutil_math.h" 		      // funções matemáticas de vetores

// Dependece files
#include "main.cuh"
#include "functions.cuh" 	 // arquivo de funções de preparação para a GPU
#include "datatypes.hpp"
#include "parser.hpp"


#define log2( x ) log(x)/log(2)

using std::cout;
using std::endl;

void PrepareSim( const char *filename,
				 SystemProperties *sisProps,
				 ParticlesValues *partValues,
				 ParticleProperties *partProps,
				 RenderParameters *renderPar ) {

	/* Usamos a estrutura de dados C++ e carregamos o arquivo de estado */
	DEMSimulation sim;
	sim.loadFromFile(filename);
	sim.printConfiguration();
	/* Agora vamos copiar para a estrutura C */

	sisProps->numParticles = sim.particles.num.x * sim.particles.num.y;

	sisProps->cubeDimension.x = sim.environment.dimension.x;
	sisProps->cubeDimension.y = sim.environment.dimension.y;
	
	sisProps->timeStep = sim.parameters.timeStep;
	
	sisProps->gravity = make_float2(sim.environment.gravity); // Transformando a gravidade de float3 para float2
	
	renderPar->imageDIMx = DIM; //TODO: Fazer uma funcão q pega o ratio do environment e aplica nos imageDIM
	renderPar->imageDIMy = DIM;

	// Zerando o raio das partículas
	for (int i = 0; i < MAX_PARTICLES_TYPES; i++){
		partProps[i].radius = 0;
	}

	//PARSER: copiando as propriedades de partículas
	for (register int i = 0; i < sim.properties.particleTypes.size(); i++)
	{
		partProps[i].mass = sim.properties.particleTypes[i].mass;
		partProps[i].radius = sim.properties.particleTypes[i].radius;
		partProps[i].collideStiffness = sim.properties.particleTypes[i].normalStiffness;
		partProps[i].collideDamping = sim.properties.particleTypes[i].normalDamping;
		partProps[i].boundaryDamping = sim.properties.particleTypes[i].boundaryDamping;
		partProps[i].colorR = sim.properties.particleTypes[i].color.x;
		partProps[i].colorG = sim.properties.particleTypes[i].color.y;
		partProps[i].colorB = sim.properties.particleTypes[i].color.z;
	}
/*
	// Propriedades da partícula 0
	partProps[0].radius = 20e-3f;
	partProps[0].mass = 1e-2;
	partProps[0].collideStiffness = 1e3;
	partProps[0].collideDamping = 0.2f;
	partProps[0].boundaryDamping = BOUNDARYDAMPING;
	partProps[0].colorR = 255;
	partProps[0].colorG = 0;
	partProps[0].colorB = 0;
	
	// Propriedades da partícula 1
	partProps[1].radius = 35e-3f;
	partProps[1].mass = 1e-2;
	partProps[1].collideStiffness = 1e3;
	partProps[1].collideDamping = 0.2f;
	partProps[1].boundaryDamping = BOUNDARYDAMPING;
	partProps[1].colorR = 0;
	partProps[1].colorG = 255;
	partProps[1].colorB = 0;
	
	// Propriedades da partícula 2
	partProps[2].radius = 30e-3f;
	partProps[2].mass = 1e-2;
	partProps[2].collideStiffness = 1e3;
	partProps[2].collideDamping = 0.2f;
	partProps[2].boundaryDamping = BOUNDARYDAMPING;
	partProps[2].colorR = 0;
	partProps[2].colorG = 0;
	partProps[2].colorB = 255;
	
	// Propriedades da partícula 3
	partProps[3].radius = 25e-3f;
	partProps[3].mass = 1e-2;
	partProps[3].collideStiffness = 1e3;
	partProps[3].collideDamping = 0.2f;
	partProps[3].boundaryDamping = BOUNDARYDAMPING;
	partProps[3].colorR = 255;
	partProps[3].colorG = 255;
	partProps[3].colorB = 255;
*/	
	float maxRadius = 0;
	for (int i = 0; i < MAX_PARTICLES_TYPES; i++){
		if (maxRadius < partProps[i].radius) maxRadius = partProps[i].radius;
	}
	
	// tamanho do quadrado que contem a esfera em PIXEL (para a saida grafica)
	renderPar->dimx = ceil(renderPar->imageDIMx/sisProps->cubeDimension.x*maxRadius)*2;
	if (renderPar->dimx < 2) renderPar->dimx = 2;
	renderPar->dimy = ceil(renderPar->imageDIMy/sisProps->cubeDimension.y*maxRadius)*2;
	if (renderPar->dimy < 2) renderPar->dimy = 2;
	
	// raio da esfera em PIXEL (para a saída grafica)
	renderPar->pRadius = renderPar->imageDIMy/sisProps->cubeDimension.y*maxRadius;

	// Bloco inicial de esferas
	float corner1[2] = {1.1, 1.1}; 				 // canto inferior esquerdo
	float corner2[2] = {8.9, 8.9}; 				  // canto superior direito
	float sideLenght[2];
	//sideLenght[0] = corner2[0] - corner1[0]; 			   // dimensao em X
	//sideLenght[1] = corner2[1] - corner1[1]; 			   // dimensao em Y
	sideLenght[0] = sim.particles.end[0] - sim.particles.start[0]; 			   // dimensao em X
	sideLenght[1] = sim.particles.end[1] - sim.particles.start[1]; 			   // dimensao em Y
	
	//uint side[2] = {X_PARTICLES, Y_PARTICLES}; // numero de partículas em X
											  // e Y (deve ser maior que 2)
	uint side[2];
	side[0] = sim.particles.num.x;
	side[1] = sim.particles.num.y;

	// Calcula o tamanho do grid arredondando para um valor que seja
	// potencia de 2. O grid deve ser de 1.2 a 3 vezes o diametro da esfera
	uint grid = sisProps->cubeDimension.x / (4.0f * maxRadius);
	uint temp = log2(grid);
	uint gridUpdate = 1 << temp;
	float cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * maxRadius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * maxRadius ) temp += 1;
	sisProps->gridSize.x = 1 << temp;
	
	grid = sisProps->cubeDimension.y / (4 * maxRadius);
	temp = log2(grid);
	gridUpdate = 1 << temp;
	cellSize = sisProps->cubeDimension.x / gridUpdate;
	if ( cellSize/2.0f <= 1.2f * maxRadius ) temp -= 1;
	else if (cellSize/2.0f >= 3.0f * maxRadius ) temp += 1;	
	sisProps->gridSize.y = 1 << temp;

	sisProps->numCells = sisProps->gridSize.x * sisProps->gridSize.y;
	
	allocateVectors(partProps, partValues, sisProps, renderPar);

	/*
	//PARSER: Copiando partículas da estrutura C++ para a estrutura C 
	std::vector<uint> merda_id;
	std::vector<float> merda_pos, merda_vel, merda_acc;
	merda_id.reserve(sim.particles.positions.size());
	merda_pos.reserve(2*sim.particles.positions.size());
	merda_vel.reserve(2*sim.particles.velocities.size());
	merda_acc.reserve(2*sim.particles.accelerations.size());

	cudaMemcpy (partValues->type1,
				&sim.particles.typeIndexes[0],
				sizeof(int)*sim.particles.typeIndexes.size(),
				cudaMemcpyHostToDevice);
	for (register int i = 0; i < sim.particles.positions.size(); i++)
	{
		merda_id.push_back(i);
	}
	cudaMemcpy (partValues->ID1, &merda_id[0], sizeof(uint)*merda_id.size(), cudaMemcpyHostToDevice);

	srand(time(NULL));
	for (register int i = 0; i < sim.particles.positions.size(); i++)
	{
		merda_pos.push_back (sim.particles.positions[i].x); //FIXME: considerar o tamanho das particulas
		merda_pos.push_back (sim.particles.positions[i].y);
		merda_vel.push_back (sim.particles.velocities[i].x);
		merda_vel.push_back (sim.particles.velocities[i].y);
		merda_acc.push_back (sim.particles.accelerations[i].x);
		merda_acc.push_back (sim.particles.accelerations[i].y);
	}
	cudaMemcpy (partValues->pos1, &merda_pos[0], sizeof(float)*merda_pos.size(), cudaMemcpyHostToDevice);
	cudaMemcpy (partValues->vel1, &merda_vel[0], sizeof(float)*merda_vel.size(), cudaMemcpyHostToDevice);
	cudaMemcpy (partValues->acc, &merda_acc[0], sizeof(float)*merda_acc.size(), cudaMemcpyHostToDevice);

*/
	// Função para definir a posição inicial das esferas
	initializeParticlePosition(partValues->pos1,
							   partValues->vel1,
							   partValues->acc,
							   partValues->ID1,
							   partValues->type1,
							   sim.particles.start,
							   sideLenght,
							   side,
							   time(NULL));

	// Screen output	
	printf("\nNumero de Particulas = %d\n",sisProps->numParticles);
	printf("grid %d x %d\n\n",sisProps->gridSize.x,sisProps->gridSize.y);
#if USE_TEX
	printf("Memoria de textura: UTILIZADA\n\n");
#else
	printf("Memoria de textura: NAO\n\n");
#endif 
}

void SimLooping( uchar4 *image, DataBlock *simBlock, int ticks ) {

	// Estruturas auxiliares
    SystemProperties *sisProps = &simBlock->sisProps;
    ParticlesValues *partValues = &simBlock->partValues;
	RenderParameters *renderPar = &simBlock->renderPar;
	TimeControl *timeCtrl = &simBlock->timeCtrl;

	// inicia o cronometro
	timeCtrl->start = clock();
	
	// para ordenarmos os vetores de posicao e velocidade sem necessidade
	// de retornarmos a variável para o vetor original, um switch entre os
	// dois vetores de posição alocados na GPU é criado. A cada iteração o
	// vetor de início e o vetor reorganizado são invertidos, reduzindo uma
	// operação de cópia
	float  *oldPos,  *oldVel;
	float *sortPos, *sortVel;
	uint  *oldID,  *oldType;
	uint *sortID, *sortType;
	
	// Integrando o programa IPS vezes antes de exibir a imagem
	for (int i = 0 ; i < timeCtrl->IPS ; i++) {

		if ((ticks + i) & 1) // quando par (FALSE) quando impar (TRUE)
		{	
			oldPos = partValues->pos1;
			oldVel = partValues->vel1;
			oldID = partValues->ID1;
			oldType = partValues->type1;
			sortPos = partValues->pos2;
			sortVel = partValues->vel2;
			sortID = partValues->ID2;
			sortType = partValues->type2;
		} else {
			oldPos = partValues->pos2;
			oldVel = partValues->vel2;
			oldID = partValues->ID2;
			oldType = partValues->type2;
			sortPos = partValues->pos1;
			sortVel = partValues->vel1;
			sortID = partValues->ID1;
			sortType = partValues->type1;
		}
		
		// Integracao no tempo (atualizacao das posicoes e velocidades)
		integrateSystem(oldPos,
			 	  		oldVel,
			 	  		partValues->acc,
			 	  		oldType,
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
									sortID,
									sortType,
									partValues->gridParticleHash,
									partValues->gridParticleIndex,
									oldPos,
									oldVel,
									oldID,
									oldType,
									sisProps->numParticles,
									sisProps->numCells);

		// Detecta a colisao das particulas e transforma a força de colisão em
		// aceleração
		collide(sortPos,
				sortVel,
				partValues->acc,
				sortType,
				partValues->cellStart,
				partValues->cellEnd,
				sisProps->numParticles,
				sisProps->numCells);

//		// Integracao no tempo (atualizacao das posicoes e velocidades)
//		integrateSystem(sortPos,
//			 	  		sortVel,
//			 	  		partValues->acc,
//						sortType,
//			 	  		sisProps->numParticles);

		timeCtrl->tempo++;
	}

	// Saida grarica quando necessario
	plotParticles(image,
				  sortPos,
				  sortType,
				  sisProps->numParticles,
				  renderPar->imageDIMx,
				  renderPar->imageDIMy);

	
	// calcula o tempo de exibição do frame
	double time = ((double)clock() - timeCtrl->start)/CLOCKS_PER_SEC;
	if (time < 0.003f) time = 0.03f;
	
	// Define o número de Interações por segundo para exibir a imagem em 
	// FPS (definida no cabeçalho) frames por segundo.
	// Após a conta, transforma o número em impar para não calcular duas
	// duas vezes a mesma iteração (por causa do switch)
	timeCtrl->IPS = floor(1.0f/time/FPS*timeCtrl->IPS);
	timeCtrl->IPS = timeCtrl->IPS | 0x0001;

}

void FinalizingSim( DataBlock *simBlock) {

	TimeControl *timeCtrl = &simBlock->timeCtrl;

    // Limpe aqui o que tiver que ser limpo
	desAllocateVectors( &simBlock->partValues );
    
   	printf("Integracoes por plot = %d\n\n",timeCtrl->IPS);
	double elapsedTime = ((double)clock() - timeCtrl->totalStart)/CLOCKS_PER_SEC;
	double simulationTime = timeCtrl->tempo * simBlock->sisProps.timeStep;
	
	printf("Duracao da simulacao = %4.2f s\n",elapsedTime);
	printf("Tempo de simulacao = %4.2f s\n\n",simulationTime);
	printf("Razao de tempo (Real/Simulado) = %3.3f\n\n",elapsedTime/simulationTime);
	
}


int main(int argc, char **argv) {
	
	// declarando estrutura de dados principal
    DataBlock simBlock;
    
    // declarando as subestruturas (apenas por facilidade)
    SystemProperties *sisProps = &simBlock.sisProps;
//    ParticleProperties *partProps = &simBlock.partProps;
    ParticlesValues *partValues = &simBlock.partValues;
	RenderParameters *renderPar = &simBlock.renderPar;
	TimeControl *timeCtrl = &simBlock.timeCtrl;
    
    ParticleProperties partProps[MAX_PARTICLES_TYPES];
    
    // Definindo que a primeira iteração será exibida
    timeCtrl->IPS = 1;
	timeCtrl->totalStart = clock();
	timeCtrl->tempo = 0;
    
    // função que define o tamanho da imagem e a estrutura que será
    // repassada para dentro do looping
    GPUAnimBitmap bitmap(DIM, DIM, &simBlock );

	// Utilizar ARGC e ARGV para pegar propriedades na linha de comando
	// ler esses comandos de um arquivo TXT externo
	// Criar uma rotina para fazer este tipo de leitura
	
	// Prepara a simulacao, define as condicoes iniciais do problema
	PrepareSim(argv[1], sisProps, partValues, partProps, renderPar);
	
	// Executa o looping até que a tecla ESC seja pressionada
    bitmap.anim_and_exit(
        (void (*)(uchar4*,void*,int))SimLooping, (void (*)(void*))FinalizingSim );

}
