#include <iostream>
#include "../includes/rapidxml.hpp"
#include "../includes/rapidxml_print.hpp"
#include "../includes/rapidxml_utils.hpp"

using namespace std;
using namespace rapidxml;

void loadUniverseFile (const char *universeFileName)
{
	file<> inputUniverse (universeFileName);
	xml_document<> doc;
	doc.parse<0>(inputUniverse.data()); // 0 = default parse flags
	//cout << doc;
	



	for (xml_node<> *node = doc.first_node(); node; node = node->next_sibling())
	{
		cout << "Node: " << node->name() << ". ";
		cout << "Value: " << node->value() << endl;
	}
}

int main (int argc, char **argv)
{
	loadUniverseFile(argv[1]);
	return 0;
}
