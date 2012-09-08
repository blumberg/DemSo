#ifndef DATATYPES_H
#define DATATYPES_H

class DEMParameters {
	public:
		float timeStep;
		float framesPerSecond;
};

class DEMEnvironment {
	public:
		float dimension[3];		//FIXME: deveria ser float3
		float gravity[3];		//FIXME: deveria ser float3
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
};

class DEMProperties {
	private:
		DEMParticleType *particleTypes;
	public:
		void addParticleType (DEMParticleType newType) {};
};

class DEMParticles {
	public:
		float *positions;
		float *velocities;
		float *accelerations;
		int   *typeIndexes;		// Se para uma partícula i, typeIndexes[i] == 1, então
								// o tipo desta partícula será o que estiver na posição
								// 1 no vetor particleTypes[]
};

class DEMSimulation {
	private:
		DEMParameters 	parameters;
		DEMEnvironment	environment;
		DEMProperties	properties;
		DEMParticles	particles;
	public:
		void loadFromFile (const char *filename);
};

#endif /* DATATYPES_H */
