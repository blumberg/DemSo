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
#include <stdio.h>
#include <ctime> 		   // biblioteca de tempo para criar o seed do rand

// CUDA includes
#include "gpu_anim.h" 				  // bib. de vizualização em tempo real
#include "cutil_math.h" 		      // funções matemáticas de vetores

// Dependece files
#include "initialStruct.cuh"
#include "main.cuh"
#include "functions.cuh" 	 // arquivo de funções de preparação para a GPU
#include "datatypes.hpp"
#include "parser.hpp"


#define log2( x ) log(x)/log(2)

using std::cout;
using std::endl;

bool file_exist (const char *filename)
{
	if (FILE *file = fopen(filename, "r"))
	{
		fclose(file);
		return true;
	} else {
		return false;
	}
}

void PrepareSim (const char *filename,
				 DataBlock *simBlock,
				 ParticleProperties *partProps)
{
	// Estruturas auxiliares
    SystemProperties *sisProps = &simBlock->sisProps;
    ParticlesValues *partValues = &simBlock->partValues;
	RenderParameters *renderPar = &simBlock->renderPar;

	/* Usamos a estrutura de dados C++ e carregamos o arquivo de estado */
	DEMSimulation sim;
	sim.loadFromFile(filename);
	sim.printConfiguration();
	/* Agora vamos copiar para a estrutura C */

/*************************************************************************/
/*************************************************************************/
// Escrevendo parametros que serão carregados pelo XML posteriormente

	sisProps->numParticles = 0;

	Quantity qtd;
	
	// verifica quantos blocos de cada tipo existe
	qtd.rectangle = sim.particles.start.size();
	qtd.triangle = sim.particles.t_pos.size();
	qtd.singleParticles = sim.particles.pos.size();
	qtd.controlParticle = USE_BIG_PARTICLE; // Quantidade máxima igual a 1
	
	Rectangle rectangle[qtd.rectangle];
	Triangle triangle[qtd.triangle];
	SingleParticles singleParts;
	ControlParticle ctrlParticle;

	// Load rectangles' properties
	for (register int i = 0; i < qtd.rectangle; i++)
	{
		rectangle[i].start = sim.particles.start[i];
		rectangle[i].end = sim.particles.end[i];
		rectangle[i].num = sim.particles.num[i];
		rectangle[i].types = sim.particles.types[i].size();
		rectangle[i].typeVec = (uint*)malloc(sizeof(uint) * rectangle[i].types);
		memcpy (rectangle[i].typeVec,
				&sim.particles.types[i][0],
				sizeof(int) * rectangle[i].types);

		sisProps->numParticles += rectangle[i].num.x * rectangle[i].num.y;
	}
	
	// Load triangles' properties
	for (register int i = 0; i < qtd.triangle; i++)
	{
		triangle[i].pos = sim.particles.t_pos[i];
		triangle[i].side = sim.particles.width[i];
		triangle[i].num = sim.particles.t_num[i];
		triangle[i].types = sim.particles.t_types[i].size();
		triangle[i].typeVec = (uint*)malloc(sizeof(uint) * triangle[i].types);
		memcpy (triangle[i].typeVec,
				&sim.particles.t_types[i][0],
				sizeof(int) * triangle[i].types);
		
		sisProps->numParticles += triangle[i].num*(triangle[i].num+1)/2;
	}

	// Load single particles' properties
	if (qtd.singleParticles > 0)
	{	
		singleParts.num = qtd.singleParticles;
		
		singleParts.pos = (float2*)malloc(sizeof(float2) * singleParts.num);
		singleParts.vel = (float2*)malloc(sizeof(float2) * singleParts.num);
		singleParts.theta = (float*)malloc(sizeof(float) * singleParts.num);
		singleParts.omega = (float*)malloc(sizeof(float) * singleParts.num);
		singleParts.type = (uint*)malloc(sizeof(uint) * singleParts.num);

		memcpy(singleParts.pos,	  &sim.particles.pos[0],   sizeof(float2) * singleParts.num);
		memcpy(singleParts.vel,   &sim.particles.vel[0],   sizeof(float2) * singleParts.num);
		memcpy(singleParts.theta, &sim.particles.theta[0], sizeof(float)  * singleParts.num);
		memcpy(singleParts.omega, &sim.particles.omega[0], sizeof(float)  * singleParts.num);
		memcpy(singleParts.type,  &sim.particles.type[0],  sizeof(int)	  * singleParts.num);
		
		sisProps->numParticles += singleParts.num;
	}

#if USE_BIG_PARTICLE
	if (qtd.controlParticle > 0)
	{
		ctrlParticle.pos = make_float2( 5,9.5);
		ctrlParticle.vel = make_float2( 0 , 0);
		ctrlParticle.theta = 0;
		ctrlParticle.omega = 0;
		ctrlParticle.type = sim.particles.controlledType;
	}
#endif

/*************************************************************************/
/*************************************************************************/

	sisProps->cubeDimension = sim.environment.dimension;

	sisProps->timeStep = sim.parameters.timeStep;

	simBlock->followedParticles = sim.particles.followedParticles;

	// Appending particles declared in <parameters> tag to the list of followed particles
	simBlock->followedParticles.insert(simBlock->followedParticles.end(),
									   sim.parameters.followedParticles.begin(),
									   sim.parameters.followedParticles.end());
	
	sisProps->gravity = sim.environment.gravity;

	sisProps->boundaryNormalStiffness = sim.environment.boundaryNormalStiffness;
	sisProps->boundaryShearStiffness = sim.environment.boundaryShearStiffness;
	
	renderPar->imageDIMy = sim.parameters.imageDIMy;
	renderPar->imageDIMx = sisProps->cubeDimension.x/sisProps->cubeDimension.y*renderPar->imageDIMy;

	// PARSER: copiando as propriedades de partículas
	for (register int i = 0; i < sim.properties.particleTypes.size(); i++)
	{
		partProps[i].mass = sim.properties.particleTypes[i].mass;
		partProps[i].radius = sim.properties.particleTypes[i].radius;
		partProps[i].normalStiffness = sim.properties.particleTypes[i].normalStiffness;
		partProps[i].shearStiffness = sim.properties.particleTypes[i].shearStiffness;
		partProps[i].normalDamping = sim.properties.particleTypes[i].normalDamping;
		partProps[i].boundaryDamping = sim.properties.particleTypes[i].boundaryDamping;
		partProps[i].frictionCoefficient = sim.properties.particleTypes[i].frictionCoefficient;
		partProps[i].colorR = sim.properties.particleTypes[i].color.x;
		partProps[i].colorG = sim.properties.particleTypes[i].color.y;
		partProps[i].colorB = sim.properties.particleTypes[i].color.z;

		partProps[i].inertia = partProps[i].mass*partProps[i].radius*partProps[i].radius / 2;
	}

	// Definindo o maior raio da simulação
	float maxRadius = 0;
	for (int i = 0; i < sim.properties.particleTypes.size(); i++)
	{
#if USE_BIG_PARTICLE
		if (i == sim.particles.controlledType) continue; // Skip controlled particle
#endif
		if (maxRadius < partProps[i].radius) maxRadius = partProps[i].radius;
	}

#if USE_BIG_PARTICLE
	float plotRadius = partProps[sim.particles.controlledType].radius;
	// tamanho do quadrado que contem a esfera controlada em PIXEL (para a saida grafica)
	renderPar->dimx = ceil(renderPar->imageDIMx/sisProps->cubeDimension.x*plotRadius)*2;
	if (renderPar->dimx < 2) renderPar->dimx = 2;
	renderPar->dimy = ceil(renderPar->imageDIMy/sisProps->cubeDimension.y*plotRadius)*2;
	if (renderPar->dimy < 2) renderPar->dimy = 2;
#endif

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


	// Alocando os vetores necessarios na CPU e GPU	
	allocateVectors(partProps, partValues, sisProps, renderPar);
	
/*************************************************************************/
/*************************************************************************/
// criando o vetor de POS, VEL, THETA, OMEGA, TYPE, ID e LOC globais
	
	uint startParticle = 0;
	uint blockSize;

	// Criando partículas avulsas
	if (qtd.singleParticles > 0)
	{
		createSingleParticles(partValues->pos1 + startParticle*2,
							  partValues->vel1 + startParticle*2,
							  partValues->theta1 + startParticle,
							  partValues->omega1 + startParticle,
							  partValues->ID1 + startParticle,
							  partValues->loc1 + startParticle,
							  partValues->type1	+ startParticle,
							  singleParts.pos,
							  singleParts.vel,
							  singleParts.theta,
							  singleParts.omega,
							  singleParts.type,
							  singleParts.num,
							  startParticle);
		
		startParticle += singleParts.num;
	}

	// Criando retangulos
	if (qtd.rectangle > 0) {
		float2 sideLenght;
		
		for (int i = 0; i < qtd.rectangle; i++)
		{
			sideLenght = rectangle[i].end - rectangle[i].start;
			
			blockSize = rectangle[i].num.x * rectangle[i].num.y;

			// Função para definir a posição inicial das esferas
			createRectangles(partValues->pos1 + startParticle*2,
							partValues->ID1 + startParticle,
							partValues->loc1 + startParticle,
							partValues->type1 + startParticle,
							rectangle[i].start,
							sideLenght,
							rectangle[i].num,
							startParticle,
							rectangle[i].types,
							rectangle[i].typeVec,
							time(NULL));

			startParticle += blockSize;
		}
	}
	
	// Criando triangulos
	if (qtd.triangle > 0)
	{
		float space;
		float height;
		
		for (int i = 0 ; i < qtd.triangle ; i++ )
		{
			blockSize = triangle[i].num * ( triangle[i].num + 1 ) / 2;
			space = triangle[i].side / ( triangle[i].num - 1 );
			height = space * sqrt(3) / 2;
			
			createTriangles(partValues->pos1 + startParticle*2,
							partValues->ID1 + startParticle,
							partValues->loc1 + startParticle,
							partValues->type1 + startParticle,
							triangle[i].pos,
							triangle[i].num,
							triangle[i].types,
							triangle[i].typeVec,
							space,
							height,
							startParticle,
							blockSize);

			startParticle += blockSize;
		}
	}
	
/*************************************************************************/
/*************************************************************************/	

#if USE_BIG_PARTICLE
	partValues->controlPos = ctrlParticle.pos;
	partValues->controlType = ctrlParticle.type;
#endif

	// Screen output	
	printf("\nNumero de Particulas = %d\n", sisProps->numParticles);
	printf("grid %d x %d\n\n", sisProps->gridSize.x, sisProps->gridSize.y);

/*************************************************************************/
/*************************************************************************/	
// Liberando as variáveis alocadas

if (qtd.rectangle > 0){
	for (int i = 0 ; i < qtd.rectangle ;  i++){
		free( rectangle[i].typeVec );
	}
}

if (qtd.triangle > 0){
	for (int i = 0 ; i < qtd.triangle ;  i++){
		free( triangle[i].typeVec );
	}
}

if (qtd.singleParticles > 0){
	free( singleParts.pos );
	free( singleParts.vel );
	free( singleParts.theta );
	free( singleParts.omega );
	free( singleParts.type );
}
/*************************************************************************/
/*************************************************************************/	

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
	float *oldTheta, *oldOmega;
	float *sortTheta, *sortOmega;
	uint  *oldID,  *oldType;
	uint *sortID, *sortType, *sortLoc;

	
	// Integrando o programa IPS vezes antes de exibir a imagem
	for (int i = 0 ; i < timeCtrl->IPS ; i++) {

		if ((ticks + i) & 1) // quando par (FALSE) quando impar (TRUE)
		{	
			oldPos = partValues->pos1;
			oldVel = partValues->vel1;
			oldTheta = partValues->theta1;
			oldOmega = partValues->omega1;
			oldID = partValues->ID1;
			oldType = partValues->type1;
			sortPos = partValues->pos2;
			sortVel = partValues->vel2;
			sortTheta = partValues->theta2;
			sortOmega = partValues->omega2;
			sortID = partValues->ID2;
			sortLoc = partValues->loc2;
			sortType = partValues->type2;
		} else {
			oldPos = partValues->pos2;
			oldVel = partValues->vel2;
			oldTheta = partValues->theta2;
			oldOmega = partValues->omega2;
			oldID = partValues->ID2;
			oldType = partValues->type2;
			sortPos = partValues->pos1;
			sortVel = partValues->vel1;
			sortTheta = partValues->theta1;
			sortOmega = partValues->omega1;
			sortID = partValues->ID1;
			sortLoc = partValues->loc1;
			sortType = partValues->type1;
		}
		
//		// Integracao no tempo (atualizacao das posicoes e velocidades)
//		integrateSystem(oldPos,
//			 	  		oldVel,
//			 	  		partValues->acc,
//						oldTheta,
//						oldOmega,
//						partValues->alpha,
//			 	  		oldType,
//			 	  		sisProps->numParticles);

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
									sortTheta,
									sortOmega,
									sortID,
									sortLoc,
									sortType,
									partValues->gridParticleHash,
									partValues->gridParticleIndex,
									oldPos,
									oldVel,
									oldTheta,
									oldOmega,
									oldID,
									oldType,
									sisProps->numParticles,
									sisProps->numCells);

		// Detecta a colisao das particulas e transforma a força de colisão em
		// aceleração
		collide(sortPos,
				sortVel,
				partValues->acc,
				sortOmega,
				partValues->alpha,
				sortType,
				partValues->cellStart,
				partValues->cellEnd,
				sisProps->numParticles,
				sisProps->numCells,
#if USE_BIG_PARTICLE
				partValues->controlPos,
				partValues->controlType,
#endif
				partValues->pressure);

		// Integracao no tempo (atualizacao das posicoes e velocidades)
		integrateSystem(sortPos,
			 	  		sortVel,
			 	  		partValues->acc,
						sortTheta,
						sortOmega,
						partValues->alpha,
			 	  		sortType,
			 	  		sisProps->numParticles);


#if USE_BIG_PARTICLE
		partValues->controlPos.y += -.005;
		partValues->controlPos.x += .003;
		if (partValues->controlPos.y < -1) partValues->controlPos.y = 10.5;
		if (partValues->controlPos.x > sisProps->cubeDimension.x + 0.5) partValues->controlPos.x = -0.5;
#endif

		timeCtrl->tempo++;
	}

	// Saida grarica quando necessario
	plotParticles(image,
				  sortPos,
				  sortTheta,
				  sortType,
				  sisProps->numParticles,
				  renderPar->imageDIMx,
				  renderPar->imageDIMy,
#if USE_BIG_PARTICLE
				  partValues->controlPos,
				  partValues->controlType,
				  renderPar->dimx,
				  renderPar->dimy,
#endif
				  partValues->pressure);

	// Escreve no arquivo de output os dados de saída
	if (!simBlock->followedParticles.empty())
		writeOutputFile (simBlock->outputFile,
						 simBlock->followedParticles,
						 sisProps->timeStep * timeCtrl->tempo, // Current elapsed time
						 (float2*)sortPos,
						 (float2*)sortVel,
						 (float2*)partValues->acc,
						 sortTheta,
						 sortOmega,
						 partValues->alpha,
						 sortID,
						 sortType,
						 sortLoc,
						 partValues->pressure);
	
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


int main(int argc, char **argv)
{	
	// Verificando arquivo de entrada (Parametros da simulacao)
	char *filename;
	
	if (argc == 2){
		if (file_exist(argv[1])){
			printf("\nUsing %s parameters file\n\n",argv[1]);
			filename = argv[1];
		}else if (file_exist("exemplos/default.dsml")){
			printf("\nFile %s does not exist, using exemplos/default.dsml file\n\n",argv[1]);
			filename = "exemplos/default.dsml";
		}else{
			printf("\nFile %s and exemlos/default.dsml does not exist.\n\nClosing simulation\n\n",argv[1]);
			return 0;
		}
	}else if (argc == 1){
		if (file_exist("exemplos/default.dsml")){
			printf("\nUsing default parameters file (exemplos/default.dsml)\n\n");
			filename = "exemplos/default.dsml";
		}else{
			printf("\nDefault file exemlos/default.dsml does not exist.\n\nClosing simulation\n\n");
			return 0;
		}
	}else{
		printf("\nToo many arguments.\n\nClosing simulation\n\n");
		return 0;
	}
	
	// declarando estrutura de dados principal
    DataBlock simBlock;
    
    // declarando as subestruturas (apenas por facilidade)
    //SystemProperties *sisProps = &simBlock.sisProps;
    //ParticlesValues *partValues = &simBlock.partValues;
	RenderParameters *renderPar = &simBlock.renderPar;
	TimeControl *timeCtrl = &simBlock.timeCtrl;
    
    ParticleProperties partProps[MAX_PARTICLES_TYPES];
    
    // Definindo que a primeira iteração será exibida
    timeCtrl->IPS = 1;
	timeCtrl->totalStart = clock();
	timeCtrl->tempo = 0;
	
	// Prepara a simulacao, define as condicoes iniciais do problema
	PrepareSim(filename, &simBlock, partProps);
	

    // função que define o tamanho da imagem e a estrutura que será
    // repassada para dentro do looping
    GPUAnimBitmap bitmap(renderPar->imageDIMx, renderPar->imageDIMy, &simBlock );
	

	// Abre arquivo de output
	simBlock.outputFile = fopen ("output.csv", "w");

	// Executa o looping até que a tecla ESC seja pressionada
    bitmap.anim_and_exit(
        (void (*)(uchar4*,void*,int))SimLooping, (void (*)(void*))FinalizingSim );

	// Fecha arquivo de output
	fclose (simBlock.outputFile);
	
	return 0;
}
