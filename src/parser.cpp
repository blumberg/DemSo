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
		if (node->name() == (string) "dimensions") cout << node->value();//env.dimension = readVector(node);
	}
	return env;
}

float * DEMParser::readVector (rapidxml::xml_node<> *root)
{
	float vect[3] = {0.0f, 0.0f, 0.0f};

	for (xml_node<> *node = root->first_node(); node; node = node->next_sibling())
	{
		if (node->name() == (string) "x") vect[0] = atof(node->value());
		else if (node->name() == (string) "y") vect[1] = atof(node->value());
		else if (node->name() == (string) "z") vect[2] = atof(node->value());
	}
	return vect;
}

int main (int argc, char **argv)
{
	DEMParser parser(argv[1]);

	DEMParameters params = parser.loadParameters();
	DEMEnvironment env = parser.loadEnvironment();

	cout << "-- Parameters" << endl;
	cout << "timeStep: " << params.timeStep << endl;
	cout << "FPS: " << params.framesPerSecond << endl;
	cout << "-- Environment" << endl;
	cout << "dimensions: " << env.dimension << endl;
	return 0;
}
