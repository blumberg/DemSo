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

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"
#include "rapidxml_utils.hpp"
#include "datatypes.hpp"
#include "parser.hpp"

using namespace std;
using namespace rapidxml;

void DEMParser::readFile (const char *filename)
{
	file<> inputFile (filename);			// Abre o arquivo e lê seu conteúdo
	xml_document<> doc;
	doc.parse<0>(inputFile.data());			// Intrepreta o XML (0 = default parse flags)
	
	rootTag = doc.first_node("simulation");	// Inicializamos o rootTag com o nó <simulation>
											// para estarmos prontos para carregar as configs.
}

DEMParser::DEMParser (const char *filename) { readFile(filename); }

DEMParameters DEMParser::loadParameters (void)
{
	DEMParameters params;
	xml_node<> *root = rootTag->first_node("parameters");

	params.followedParticle = -1;	// Valor indicando que não há partícula a ser seguida
	// Percorrento os nós filhos do nó principal (<parameters>)
	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("timestep")) params.timeStep = atof(node->value());
		else if (node->name() == string("imageHeight")) params.imageDIMy = atoi(node->value());
		else if (node->name() == string("follow")) params.followedParticle = atoi(node->value());
	}

	return params;
}

DEMEnvironment DEMParser::loadEnvironment (void)
{
	DEMEnvironment env;
	xml_node<> *root = rootTag->first_node("environment");

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("dimensions")) env.dimension = readVector(node);
		else if (node->name() == string("gravity")) env.gravity = readVector(node);
		else if (node->name() == string("stiffness"))
		{
			xml_attribute<> *attr = node->first_attribute("dir");
			if (attr->value() == string("normal")) env.boundaryNormalStiffness = atof(node->value());
			else if (attr->value() == string("shear")) env.boundaryShearStiffness = atof(node->value());
			else throw string("Unrecognized stiffness direction");
		}
		else if (node->name() == string("damping")) env.boundaryDamping = atof(node->value());
		else if (node->name() == string("friction")) env.frictionCoefficient = atof(node->value());
	}
	return env;
}

DEMProperties DEMParser::loadProperties (void)
{
	DEMProperties prop;
	xml_node<> *root = rootTag->first_node("properties");

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("particletype")) prop.addParticleType (loadParticleType (node));
	}
	return prop;
}

DEMParticleType DEMParser::loadParticleType (xml_node<> *root)
{
	DEMParticleType ptype;
	
	if (root->first_attribute("id")) ptype.id = root->first_attribute("id")->value();
	if (root->first_attribute("name")) ptype.name = root->first_attribute("name")->value();
	if (root->first_attribute("color"))
	{
		string color = root->first_attribute("color")->value();
		char red[4];
		char green[4];
		char blue[4];
		strcpy (red,"0x");
		strcpy (green,"0x");
		strcpy (blue,"0x");
		strcat (red,color.substr(0,2).c_str());
		strcat (green,color.substr(2,2).c_str());
		strcat (blue,color.substr(4,2).c_str());
		
		ptype.color.x = atof(red);
		ptype.color.y = atof(green);
		ptype.color.z = atof(blue);
	}
	else ptype.color = make_float3(255,255,255);

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("mass")) ptype.mass = atof(node->value());
		else if (node->name() == string("radius")) ptype.radius = atof(node->value());
		else if (node->name() == string("stiffness"))
		{
			xml_attribute<> *attr = node->first_attribute("dir");
			if (attr->value() == string("normal")) ptype.normalStiffness = atof(node->value());
			else if (attr->value() == string("shear")) ptype.shearStiffness = atof(node->value());
			else throw string("Unrecognized stiffness direction");
		}
		else if (node->name() == string("damping")) ptype.normalDamping = atof(node->value());
		else throw string("Unrecognized tag inside <particletype>");
	}
	return ptype;
}

DEMParticles DEMParser::loadParticles (DEMProperties *properties)
{
	DEMParticles parts;
	xml_node<> *root = rootTag->first_node("particles");

	// Se não foram definidos típos de partículas, abortar.
	if (properties->particleTypes.empty())
		throw string("No particle types defined, aborting.");

	// Percorre os nós filhos de <particles> de 1 em 1
	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		// Caso for encontrado um bloco de partículas
		if (node->name() == string("block"))
		{
			parts.start = make_float2(0);
			parts.end = make_float2(0);
			parts.num = make_float2(0);
			if (node->first_node("start")) parts.start = readVector(node->first_node("start"));
			if (node->first_node("end")) parts.end = readVector(node->first_node("end"));
			if (node->first_node("num")) parts.num = readVector(node->first_node("num"));
		}

		// Caso for encontrada uma partícula avulsa (tag <particle>)
		if (node->name() == string("particle"))
		{
			// Tipo da partícula
			int type;

			// Tenta recuperar o atributo particletype
			xml_attribute<> *attr = node->first_attribute("particletype");

			// Se ele está presente, seu valor contém o ID do tipo da partícula,
			// então, procuramos qual é o index do típo de partícula a partir do ID
			if (attr) type = properties->particleTypeIndexById(attr->value());

			// Senão utilizamos o primeiro típo de partículas definido
			else {
				type = 0;
				cout << "Particle type not specified, defaulting to first one: "
					 << properties->particleTypes[type].id << endl;
			}
			parts.type.push_back(type);

			// Se a posição inicial (tag <pos>) foi definida, adiciona ao vetor
			// de posições, senão, aborta.
			xml_node<> *node_pos = node->first_node("pos");
			if (node_pos)
				parts.pos.push_back(readVector(node_pos));			
			else
				throw "Missing <pos> tag on particle definition.";

			// Lê a velocidade inicial (tag <vel>)
			xml_node<> *node_vel = node->first_node("vel");
			if (node_vel)
				parts.vel.push_back(readVector(node_vel));
			else
				parts.vel.push_back(make_float2(0.0f));

			// Lê a posição angular inicial (tag <theta>)
			xml_node<> *node_theta = node->first_node("theta");
			if (node_theta)
				parts.theta.push_back(atof(node_theta->value()));
			else
				parts.theta.push_back(0.0f);

			// Lê a velocidade angular inicial (tag <omega>)
			xml_node<> *node_omega = node->first_node("omega");
			if (node_omega)
				parts.omega.push_back(atof(node_omega->value()));
			else
				parts.omega.push_back(0.0f);
		}
	}

	return parts;
}

float2 DEMParser::readVector (xml_node<> *root)
{
	float2 vect = make_float2(0.0f);
	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("x")) vect.x = atof(node->value());
		else if (node->name() == string("y")) vect.y = atof(node->value());
	}
	return vect;
}

/*int main (int argc, char **argv)
{
	DEMSimulation sim;
	sim.loadFromFile (argv[1]);
	sim.printConfiguration();

	return 0;
}*/
