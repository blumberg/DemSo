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
#include <sstream>
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

	// Percorrento os nós filhos do nó principal (<parameters>)
	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == string("timestep")) params.timeStep = atof(node->value());
		else if (node->name() == string("imageHeight")) params.imageDIMy = atoi(node->value());
		else if (node->name() == string("follow")) params.followedParticles.push_back(atoi(node->value()));
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
			if (attr)
			{
				if (attr->value() == string("normal")) ptype.normalStiffness = atof(node->value());
				else if (attr->value() == string("shear")) ptype.shearStiffness = atof(node->value());
				else throw string("Unrecognized stiffness direction");
			} else throw string("Stiffness direction not specified");
		}
		else if (node->name() == string("damping"))
		{
			xml_attribute<> *attr = node->first_attribute("type");
			if (attr)
			{
				if (attr->value() == string("boundary")) ptype.boundaryDamping = atof(node->value());
				else throw string("Unrecognized damping type");
			}
			else ptype.normalDamping = atof(node->value());
		}
		else if (node->name() == string("friction")) ptype.frictionCoefficient = atof(node->value());
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
		// Caso for encontrado um retângulo
		if (node->name() == string("rectangle"))
		{
			// Tenta recuperar o atributo particletype
			xml_attribute<> *attr = node->first_attribute("particletype");

			// Se ele for encontrado, lê os tipos e guarda no vetor
			if (attr)
			{
				// Le o ID dos tipos e guarda num vetor de IDs
				vector<string> typeids = readCSVLine(attr->value());

				// Percorre o vetor de IDs e popula o vetor dos respectivos índices
				vector<int> typeindexes;
				for (register int i = 0; i < typeids.size(); i++)
					typeindexes.push_back(properties->particleTypeIndexById(typeids[i]));

				// Adiciona o vetor de tipos à lista
				parts.types.push_back(typeindexes);
			}
			// Senão usa um tipo só: o primeiro a ser definido
			else {
				vector<int> typeindexes(1,0);
				parts.types.push_back(typeindexes);
				cout << "Particle type not specified for rectangle, defaulting to first one: "
					 << properties->particleTypes[0].id << endl;
			}

			// Leitura das propriedades do retângulo ou término do programa caso falte alguma
			if (node->first_node("start")) parts.start.push_back(readVector(node->first_node("start")));
			else throw string("<start> tag not found inside rectangle");

			if (node->first_node("end")) parts.end.push_back(readVector(node->first_node("end")));
			else throw string("<end> tag not found inside rectangle");

			if (node->first_node("num")) parts.num.push_back(make_uint2(readVector(node->first_node("num"))));
			else throw string("<num> tag not found inside rectangle");
		}

		// Caso for encontrado um triângulo
		if (node->name() == string("triangle"))
		{
			// Tenta recuperar o atributo particletype
			xml_attribute<> *attr = node->first_attribute("particletype");

			// Se ele for encontrado, lê os tipos e guarda no vetor
			if (attr)
			{
				// Le o ID dos tipos e guarda num vetor de IDs
				vector<string> typeids = readCSVLine(attr->value());

				// Percorre o vetor de IDs e popula o vetor dos respectivos índices
				vector<int> typeindexes;
				for (register int i = 0; i < typeids.size(); i++)
					typeindexes.push_back(properties->particleTypeIndexById(typeids[i]));

				// Adiciona o vetor de tipos à lista
				parts.t_types.push_back(typeindexes);
			}
			// Senão usa um tipo só: o primeiro a ser definido
			else {
				vector<int> typeindexes(1,0);
				parts.t_types.push_back(typeindexes);
				cout << "Particle type not specified for triangle, defaulting to first one: "
					 << properties->particleTypes[0].id << endl;
			}

			// Leitura das propriedades do triângulo ou término do programa caso falte alguma
			if (node->first_node("pos")) parts.t_pos.push_back(readVector(node->first_node("pos")));
			else throw string("<pos> tag not found inside triangle");

			if (node->first_node("width")) parts.width.push_back(atof(node->first_node("width")->value()));
			else throw string("<width> tag not found inside triangle");

			if (node->first_node("num")) parts.t_num.push_back(atoi(node->first_node("num")->value()));
			else throw string("<num> tag not found inside triangle");
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

// FIXME: por enquanto espaços em branco não são desconsiderados
vector<string> DEMParser::readCSVLine (string line)
{
	stringstream ss(line);
	vector<string> values;

	while (ss.good())
	{
		string substr;
		getline (ss, substr, ',');
		values.push_back(substr);
	}

	return values;
}

/*int main (int argc, char **argv)
{
	DEMSimulation sim;
	sim.loadFromFile (argv[1]);
	sim.printConfiguration();

	return 0;
}*/
