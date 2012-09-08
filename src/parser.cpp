#include <iostream>
#include <string>
#include "../includes/rapidxml.hpp"
#include "../includes/rapidxml_print.hpp"
#include "../includes/rapidxml_utils.hpp"
#include <cstdlib>
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
		if (node->name() == (string) "timestep") params.timeStep = atof(node->value());
		else if (node->name() == (string) "fps") params.framesPerSecond = atof(node->value());
	}

	return params;
}

DEMEnvironment DEMParser::loadEnvironment (void)
{
	DEMEnvironment env;
	xml_node<> *root = rootTag->first_node("environment");

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "dimensions") env.dimension = readVector(node);
		if (node->name() == (string) "gravity") env.gravity = readVector(node);
	}
	return env;
}

DEMProperties DEMParser::loadProperties (void)
{
	DEMProperties prop;
	xml_node<> *root = rootTag->first_node("properties");

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "particletype") prop.addParticleType (loadParticleType (node));
	}
	return prop;
}

DEMParticleType DEMParser::loadParticleType (xml_node<> *root)
{
	DEMParticleType ptype;
	
	if (root->first_attribute("id")) ptype.id = root->first_attribute("id")->value();
	if (root->first_attribute("name")) ptype.name = root->first_attribute("name")->value();

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "mass") ptype.mass = atof(node->value());
		else if (node->name() == (string) "radius") ptype.radius = atof(node->value());
		else if (node->name() == (string) "stiffness")
		{
			xml_attribute<> *attr = node->first_attribute("dir");
			if (attr->value() == (string) "normal") ptype.normalStiffness = atof(node->value());
			//else if (attr->value() == (string) "tangent") ptype.tangentStiffness = atof(node->value());
			else throw "Unrecognized stiffness direction";
		}
		else if (node->name() == (string) "damping")
		{
			xml_attribute<> *attr = node->first_attribute("type");
			if (attr && attr->value() == (string) "boundary")
				ptype.boundaryDamping = atof(node->value()); //FIXME: test if tangent or normal
			else {
				attr = node->first_attribute("dir");
				if (attr->value() == (string) "normal") ptype.normalDamping = atof(node->value());
				//else if (attr->value() == (string) "tangent") ptype.tangentDamping = atof(node->value());
				else throw "Unrecognized damping direction";
			}
			
		}
		else throw "Unrecognized tag inside <particletype>";
	}
	return ptype;
}

DEMParticles DEMParser::loadParticles (void)
{
	DEMParticles parts;
	xml_node<> *root = rootTag->first_node("particles");

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "block") cout << "Achei um block!\n"; //parts.addParticles(loadBlock(node));
	}

	return parts;
}

float3 DEMParser::readVector (xml_node<> *root)
{
	float3 vect = make_float3(0.0f);

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "x") vect.x = atof(node->value());
		else if (node->name() == (string) "y") vect.y = atof(node->value());
		else if (node->name() == (string) "z") vect.z = atof(node->value());
	}
	return vect;
}

int main (int argc, char **argv)
{
	DEMSimulation sim;
	sim.loadFromFile (argv[1]);
	sim.printConfiguration();

	return 0;
}
