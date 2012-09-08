#ifndef PARSER_H
#define PARSER_H
#include "../includes/rapidxml.hpp"
#include "datatypes.hpp"

using namespace rapidxml;

class DEMParser {
	private:
		xml_node<> *rootTag;
		float *readVector (xml_node<> *);
	public:
		DEMParser () {};
		DEMParser (const char *);

		void readFile (const char *);
		DEMParameters loadParameters (void);
		DEMEnvironment loadEnvironment (void);
		//DEMProperties loadProperties (void);
		//DEMParticles loadParticles (void);
};

#endif /* PARSER_H */
