#include <curand_kernel.h>
#include <cstdint>

struct hut_exp_control{
	double xlim[2];
	double ylim[2];		//dimensions of the domain
	double sig_acc[2];	//parameters for attraction
	double pnet;		//probability of being blocked by the net
	double phut; //prob of nbot exiting the hut
	double eps;			//distance to the host treated as bite
	double mu;			// natural death rate
	float tmax;			//max time spend in taxis, in hours
};

struct RandomGeneratorState{
	curandState_t state_xinit;
	curandState_t state_yinit;
	curandState_t state_active;
	curandState_t state_dead;
	curandState_t state_d1;
	curandState_t state_d2;
	curandState_t state_th1;
	curandState_t state_th2;
	curandState_t state_attr;
	curandState_t state_acc;
};

struct ComputationContext{

	//Pointers to host memory blocks
	char* hostIn;
	char* hostDead;
	char* hostTrap;
	char* hostFed;
	char* hostUnfDead;

	//Pointers to device memory blocks
	RandomGeneratorState* p_rnd_states;
	char* devIn;
	char* devDead;
	char* devTrap;
	char* devFed;
	char* devUnfDead;
};


void initialize(uint32_t repetitions, uint32_t num_experiments, ComputationContext& cmpt_ctx);
void CallKernel(uint32_t repetitions, uint32_t num_experiments, const hut_exp_control& experimental_params, const ComputationContext& cmpt_ctx);
void shutdown(const ComputationContext& cmpt_ctx);
