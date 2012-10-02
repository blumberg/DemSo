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

#ifndef INITIAL_STRUCT_CUH
#define INITIAL_STRUCT_CUH

// Estrutura para criar um retangulo de partículas
struct retangle {

	uint2 num; 		// Numero de partículas em X e Y
	float2 start;	// Coordenada do canto inferior esquerdo do quadrado
	float2 end;		// Coordenada do canto superior direito do quadrado

	uint *typeVec; 	// Vetor com os tipos de partículas associados

};

// Estrutura para criar um triangulo de partículas
struct triangle {

	uint num;		// Numero de partículas no lado do triangulo
	float2 pos;		// Posicao central inferior do triangulo
	float2 side;	// Comprimento do lado do triangulo

	uint *typeVec; 	// Vetor com os tipos de partículas associados
	
};

// Entrada na mão de partículas
struct userDefine {

	uint num;		// Numero total de particulas adicionadas
	float2 *pos;	// Posicao das partículas
	float2 *vel;	// Velocidade das partículas
	float *theta;	// posição angular das partículas
	float *omega;	// velocidade angular das partículas
	uint *type; 	// vetor de tipo das partículas

};

// Partícula controlada externamente (Não integrada)
struct controlParticle {

	float2 pos;
	float2 vel;
	float theta;
	float omega;
	uint type;

}


#endif
