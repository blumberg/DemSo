#ifndef DATATYPES_H
#define DATATYPES_H
#include <vector>
#include <math.h>
#include "../includes/cutil_math.h"

class DEMParameters {
	public:
		float timeStep;
		float framesPerSecond;
};

class DEMEnvironment {
	public:
		float3 dimension;
		float3 gravity;
};

class DEMParticleType {
	public:
		std::string id;
		std::string name;
		float mass;
		float radius;
		float normalStiffness;
		float normalDamping;
		float boundaryDamping;
		DEMParticleType () {};
};

class DEMProperties {
	public:
		std::vector<DEMParticleType> particleTypes;
		void addParticleType (DEMParticleType);
		int particleTypeIndexById (std::string);
};

class DEMParticles {
	public:
		std::vector<float3> positions;
		std::vector<float3> velocities;
		std::vector<float3> accelerations;
		std::vector<int> typeIndexes;	// Se para uma partícula i, typeIndexes[i] == 1, então
										// o tipo desta partícula será o que estiver na posição
										// 1 no vetor particleTypes[]
		void addParticles (DEMParticles);
		void generateBlock (int, float3, float3, float3, DEMProperties *);
};

class DEMSimulation {
	private:
		DEMParameters 	parameters;
		DEMEnvironment	environment;
		DEMProperties	properties;
		DEMParticles	particles;
	public:
		void loadFromFile (const char *filename);
		void printConfiguration ();
};

#endif /* DATATYPES_H */
