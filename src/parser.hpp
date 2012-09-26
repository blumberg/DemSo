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

#ifndef PARSER_H
#define PARSER_H
#include "../includes/rapidxml.hpp"
#include "../includes/cutil_math.h"
#include "datatypes.hpp"

using namespace rapidxml;

class DEMParser {
	private:
		xml_node<> *rootTag;
		float2 readVector (xml_node<> *);
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
