#ifndef PARSER_H
#define PARSER_H
#include "../includes/rapidxml.hpp"
#include "../includes/cutil_math.h"
#include "datatypes.hpp"

using namespace rapidxml;

class DEMParser {
	private:
		xml_node<> *rootTag;
		float3 readVector (xml_node<> *);
		DEMParticleType loadParticleType (xml_node<> *);
		DEMParticles loadBlock (xml_node<> *, DEMProperties *);
	public:
		DEMParser () {};
		DEMParser (const char *);

		void readFile (const char *);
		DEMParameters loadParameters (void);
		DEMEnvironment loadEnvironment (void);
		DEMProperties loadProperties (void);
		DEMParticles loadParticles (DEMProperties *);
};

#endif /* PARSER_H */
