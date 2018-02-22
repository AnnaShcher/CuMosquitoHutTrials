#include "CuMosquitoHutTrialsITN.h"
#pragma comment(lib, "cudart")

#include <ctime>


#define PI 3.14159265359
#define e 2.718281828459
#define FINAL_TIME 36000//duration of experiment in seconds
#define NITERATIONS 18000
#define NSTEPS 19000
#define HUT_SIZE 1.5	//meters
#define NET_SIZE 0.75
#define STD_A 80/3.0	//std of attraction kernel
#define EPS 2.2204e-16
#define CONC_NET_SURF pow(e,(-0.5*pow(NET_SIZE,2.0)/pow(STD_A,2.0)))	//concentration of CO2 on net surface
#define CUDA_BLOCK_DIM	16	//dimension of a single rectangular CUDA execution block
#define CUDA_BLOCK_SIZE CUDA_BLOCK_DIM*CUDA_BLOCK_DIM	//number of CUDA threads in a single block

template<typename T>
T min(T a, T b){ return a<b ? a : b; }
template<typename T>
T max(T a, T b){ return a>b ? a : b; }


__global__ void setup_kernel(RandomGeneratorState* p_rnd_states, uint64_t seed_spin_up, uint32_t repetitions, uint32_t experiments)
{
	unsigned int id_x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int id_y = blockIdx.y * blockDim.y + threadIdx.y;

	if (id_x >= experiments || id_y >= repetitions) return;

	unsigned int id = id_y * experiments + id_x;

	/* Each thread gets same seed, a different sequence number,
	no offset */
	curand_init(seed_spin_up, id + 0, 0, &p_rnd_states[id].state_xinit);
	curand_init(seed_spin_up, id + 1, 0, &p_rnd_states[id].state_yinit);
	curand_init(seed_spin_up, id + 2, 0, &p_rnd_states[id].state_active);
	curand_init(seed_spin_up, id + 3, 0, &p_rnd_states[id].state_dead);
	curand_init(seed_spin_up, id + 4, 0, &p_rnd_states[id].state_d1);
	curand_init(seed_spin_up, id + 5, 0, &p_rnd_states[id].state_d2);
	curand_init(seed_spin_up, id + 6, 0, &p_rnd_states[id].state_th1);
	curand_init(seed_spin_up, id + 7, 0, &p_rnd_states[id].state_th2);
	curand_init(seed_spin_up, id + 8, 0, &p_rnd_states[id].state_attr);
	curand_init(seed_spin_up, id + 9, 0, &p_rnd_states[id].state_acc);
}


__global__ void hut_exp_kernel(char* in, char* dead, char* trap, char* fed, char* unf_dead, double* Ctot, double* net_cont,
	const hut_exp_ITN xin,
	RandomGeneratorState* p_rnd_states,
	uint32_t repetitions, uint32_t experiments)
{
	unsigned int id_x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int id_y = blockIdx.y * blockDim.y + threadIdx.y;

	if (id_x >= experiments || id_y >= repetitions) return;

	unsigned int id = id_y * experiments + id_x;
    unsigned int id_CTOT, id_ncont;

	bool entered, out;
	int tmax = (int)(3600 * xin.tmax);
	int tindrs = 0;
    int step = 0;
	bool  attr, acc, move;
	double death, p_attr, p_rep;
	double death_rate = xin.mu / 1800 / 34;
	double d1, th1, d2, th2;
	double x, y, xnew, ynew, dold, dnew, cold, cnew;
	double sig_acc, conc_chem = 0.0, count_ncont = 0.0, death_rate_psn;
    double sig_max = xin.sig_acc[0] + xin.sig_acc[1] * 80.0;
	bool cond[7] = { false };

	x = xin.xlim[0] + (xin.xlim[1] - xin.xlim[0]) * curand_uniform(&p_rnd_states[id].state_xinit);
	y = xin.ylim[0] + (xin.ylim[1] - xin.ylim[0]) * curand_uniform(&p_rnd_states[id].state_yinit);
	dold = sqrt(pow(x, 2.0) + pow(y, 2.0));
	while (dold < NET_SIZE)
	{
		x = xin.xlim[0] + (xin.xlim[1] - xin.xlim[0]) * curand_uniform(&p_rnd_states[id].state_xinit);
		y = xin.ylim[0] + (xin.ylim[1] - xin.ylim[0]) * curand_uniform(&p_rnd_states[id].state_yinit);
		dold = sqrt(pow(x, 2.0) + pow(y, 2.0));
	}
	cold = pow(e, (-0.5 * pow(dold, 2.0) / pow(STD_A, 2.0)));
	entered = dold < HUT_SIZE;
	cond[0] = entered;
	//variables, conditions
	//[0      1       2    3   4                5          
	//[In/Out Trapped Dead Fed taxis/kinesis  inside_net  


	for (int n = 2; n <= FINAL_TIME; n += 2)
	{
		//natural death
		death = curand_uniform_double(&p_rnd_states[id].state_dead);
		//death_rate_psn = (1 - death_rate) * pow(e, -1.0 / xin.alpha_p / conc_chem);
		//cond[2] = death < death_rate + death_rate_psn;//select the 'fortune'
        if (conc_chem == 0)
            cond[2] = death < death_rate;
        else
            //cond[2] = death < death_rate  + 1.0/(1.0 + pow(e,-(conc_chem - xin.alpha_p) / xin.alpha_d));//select the 'fortune'
            cond[2] = death < death_rate + xin.alpha_p*conc_chem;//select the 'fortune'
        //cond[2] = death < death_rate + 1.0 / (1.0 + exp(-(conc_chem - xin.alpha_d) / xin.alpha_p));
		//if (cond[2])
		//	break;

		move = !cond[1] && !cond[2];	//not trapped,& not dead, & not resting

		//candidate step
		d1 = 0.4 + 0.1* curand_normal_double(&p_rnd_states[id].state_d1);
		th1 = 2 * PI * curand_uniform(&p_rnd_states[id].state_th1);
		d2 = 0.4 + 0.1* curand_normal_double(&p_rnd_states[id].state_d1);
		th2 = 2 * PI * curand_uniform(&p_rnd_states[id].state_th1);
		xnew = x + d1 * cos(th1) + d2 * cos(th2);
		ynew = y + d1 * sin(th1) + d2 * sin(th2);

		//measuring concentration for new position
		dnew = sqrt(pow(xnew, 2.0) + pow(ynew, 2.0));
		cnew = pow(e, (-0.5 * pow(dnew, 2.0) / pow(STD_A, 2.0)));
       
        sig_acc = xin.sig_acc[0] + xin.sig_acc[1] * dold;//scaling factor for attraction
        //sig_acc = xin.sig_acc[0] + xin.sig_acc[1] * dold;

		p_attr = min(1.0, pow(e, (cnew - cold) / sig_acc));//prob of acceptance
		if (!cond[5])
			p_rep = xin.r * min(1.0, (1.0 - 1.0 / (1.0 + pow(e, -(dnew - xin.d50) / xin.s))));//concentration of repellent

		else
			p_rep = xin.r * min(1.0, (1.0 / (1.0 + pow(e, -(dnew - xin.d50) / xin.s))));//concentration of repellent
		attr = cond[4] && curand_uniform_double(&p_rnd_states[id].state_attr) < (1 - p_rep) ||
			curand_uniform_double(&p_rnd_states[id].state_attr) < p_attr * (1 - p_rep);//attraction or kinesis

        id_CTOT = repetitions * experiments* step + id_y * experiments + id_x;
        id_ncont = id * NITERATIONS + step;
        //Ctot[id_CTOT] = conc_chem;
        //net_cont[id_CTOT] = count_ncont;
        Ctot[id_CTOT] = xin.d50_NetCont;
        net_cont[id_CTOT] = xin.s_NetCont;
        step = step + 1;
        conc_chem = conc_chem -xin.d50_NetCont*conc_chem;

		if (!cond[5] && dnew <= NET_SIZE)
		{
			acc = curand_uniform_double(&p_rnd_states[id].state_acc) < 1.0 - xin.pnet;
			if (move && attr && !acc)
			{
				x = NET_SIZE + EPS;
				y = 0;//hiting the sufice
				cold = CONC_NET_SURF;
				dold = NET_SIZE + EPS;
				conc_chem = conc_chem + xin.r * min(1.0, (1.0 - 1.0 / (1.0 + pow(e, -(dold - xin.d50) / xin.s))));
                count_ncont = count_ncont + 1.0;
			}
		}
		else if (cond[0] && dnew > HUT_SIZE)
			acc = curand_uniform_double(&p_rnd_states[id].state_active) < xin.phut;
		else
			acc = true;//taking into account net barrier
		if (move && attr && acc)
		{
			x = xnew;
			y = ynew;
			cold = cnew; //and resp. site & conc. values
			dold = dnew;
		}
		out = !cond[0];
		cond[0] = cond[0] || (dold <= HUT_SIZE);//inside
		entered = out && cond[0];
		cond[1] = cond[1] || (cond[0] && (dold >= HUT_SIZE));//mark trapped mosquitoes
		cond[3] = cond[3] || (dold < xin.eps);
		cond[5] = dold < NET_SIZE;
		if (cond[0] && !cond[2] && !cond[1])
			tindrs = tindrs + 2;//increment time spent indoors
		cond[4] = cond[4] || (tindrs >= tmax) || cond[3];
	}
	death = curand_uniform_double(&p_rnd_states[id].state_dead);
	//cond[2] = cond[2] || (death < 24 * xin.mu / 34 + xin.alpha_p * pow(conc,0.7) * 24 * 1800);	//select the 'fortune'
	//death_rate_psn = (1 - death_rate) * pow(e, -1.0 / xin.alpha_p / conc_chem);
    //cond[2] = cond[2] || (death < 24 * xin.mu / 34 + death_rate_psn * 24 * 1800);	//select the 'fortune'
    if(conc_chem == 0)
        cond[2] = cond[2] || (death < 24 * xin.mu / 34);	//select the 'fortune'
    else
	    //cond[2] = cond[2] || (death < 24 * xin.mu / 34 + 1.0 / (1.0 + pow(e,-(conc_chem - xin.alpha_p) / xin.alpha_d)) * 24.0 * 1800.0);	//select the 'fortune
        cond[2] = cond[2] || (death < 24 * xin.mu / 34 + +xin.alpha_p*conc_chem * 24.0 * 1800.0);	//select the 'fortune'
	in[id] = cond[0];
	dead[id] = cond[0] && cond[2];
	trap[id] = cond[1];
	fed[id] = cond[3];
	unf_dead[id] = cond[2] && !cond[3];
	/* Copy state back to global memory */
	/* Store results */
}


inline void setupGrid(uint32_t thread_row_num, uint32_t thread_col_num, dim3& grid_dim, dim3& block_dim)
{
	block_dim = dim3(CUDA_BLOCK_DIM, CUDA_BLOCK_DIM);

	uint32_t grid_dim_y = (thread_row_num % CUDA_BLOCK_DIM == 0) ?
		thread_row_num / CUDA_BLOCK_DIM :
		(thread_row_num & ~(CUDA_BLOCK_DIM - 1)) / CUDA_BLOCK_DIM + 1;

	uint32_t grid_dim_x = (thread_col_num % CUDA_BLOCK_DIM == 0) ?
		thread_col_num / CUDA_BLOCK_DIM :
		(thread_col_num & ~(CUDA_BLOCK_DIM - 1)) / CUDA_BLOCK_DIM + 1;
	grid_dim = dim3(grid_dim_x, grid_dim_y);
}


void initialize(uint32_t repetitions, uint32_t num_experiments, ComputationContext& cmpt_ctx)
{
	dim3 grid_dim, block_dim;
	setupGrid(repetitions, num_experiments, grid_dim, block_dim);

	cmpt_ctx.hostIn = (char*)malloc(repetitions * num_experiments);
	cmpt_ctx.hostDead = (char*)malloc(repetitions * num_experiments);
	cmpt_ctx.hostTrap = (char*)malloc(repetitions * num_experiments);
	cmpt_ctx.hostFed = (char*)malloc(repetitions * num_experiments);
	cmpt_ctx.hostUnfDead = (char*)malloc(repetitions * num_experiments);
	cmpt_ctx.hostCtot = (double*)malloc(repetitions * num_experiments * NITERATIONS*sizeof(double));
    cmpt_ctx.hostNetCont = (double*)malloc(repetitions * num_experiments * NITERATIONS * sizeof(double));
	/* Set results to 0 */
	/* Allocate space for prng states on device */

	cudaSetDevice(0);
	cudaDeviceReset();

	cudaMalloc(&cmpt_ctx.p_rnd_states, repetitions * num_experiments * sizeof(RandomGeneratorState));

	/* Setup prng states */
	setup_kernel << <grid_dim, block_dim >> >(cmpt_ctx.p_rnd_states, time(NULL), repetitions, num_experiments);

	cudaMalloc(&cmpt_ctx.devIn, repetitions * num_experiments);
	cudaMalloc(&cmpt_ctx.devDead, repetitions * num_experiments);
	cudaMalloc(&cmpt_ctx.devTrap, repetitions * num_experiments);
	cudaMalloc(&cmpt_ctx.devFed, repetitions * num_experiments);
	cudaMalloc(&cmpt_ctx.devUnfDead, repetitions * num_experiments);
	cudaMalloc(&cmpt_ctx.devCtot, repetitions * num_experiments* NITERATIONS * sizeof(double));
    cudaMalloc(&cmpt_ctx.devNetCont, repetitions * num_experiments * NITERATIONS * sizeof(double));
}


void CallKernel(uint32_t repetitions, uint32_t num_experiments, const hut_exp_ITN& experimental_params, const ComputationContext& cmpt_ctx)
{
	if (cudaGetLastError() != CUDA_SUCCESS) return;

	dim3 grid_dim, block_dim;
	setupGrid(repetitions, num_experiments, grid_dim, block_dim);

	hut_exp_kernel << <grid_dim, block_dim >> >(cmpt_ctx.devIn, cmpt_ctx.devDead, cmpt_ctx.devTrap, cmpt_ctx.devFed, cmpt_ctx.devUnfDead, cmpt_ctx.devCtot, cmpt_ctx.devNetCont,
		experimental_params, cmpt_ctx.p_rnd_states, repetitions, num_experiments);

	cudaMemcpy(cmpt_ctx.hostIn, cmpt_ctx.devIn, repetitions*num_experiments, cudaMemcpyDeviceToHost);
	cudaMemcpy(cmpt_ctx.hostDead, cmpt_ctx.devDead, repetitions*num_experiments, cudaMemcpyDeviceToHost);
	cudaMemcpy(cmpt_ctx.hostTrap, cmpt_ctx.devTrap, repetitions*num_experiments, cudaMemcpyDeviceToHost);
	cudaMemcpy(cmpt_ctx.hostFed, cmpt_ctx.devFed, repetitions*num_experiments, cudaMemcpyDeviceToHost);
	cudaMemcpy(cmpt_ctx.hostUnfDead, cmpt_ctx.devUnfDead, repetitions*num_experiments, cudaMemcpyDeviceToHost);
	cudaMemcpy(cmpt_ctx.hostCtot, cmpt_ctx.devCtot, repetitions*num_experiments* NITERATIONS * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(cmpt_ctx.hostNetCont, cmpt_ctx.devNetCont, repetitions * num_experiments * NITERATIONS * sizeof(double), cudaMemcpyDeviceToHost);
}


void shutdown(const ComputationContext& cmpt_ctx)
{
	//Free host memory
	free(cmpt_ctx.hostIn);
	free(cmpt_ctx.hostDead);
	free(cmpt_ctx.hostTrap);
	free(cmpt_ctx.hostFed);
	free(cmpt_ctx.hostUnfDead);
	free(cmpt_ctx.hostCtot);
    free(cmpt_ctx.hostNetCont);

	//Free GPU memory
	cudaFree(cmpt_ctx.devIn);
	cudaFree(cmpt_ctx.devDead);
	cudaFree(cmpt_ctx.devTrap);
	cudaFree(cmpt_ctx.devFed);
	cudaFree(cmpt_ctx.devUnfDead);
	cudaFree(cmpt_ctx.devCtot);
    cudaFree(cmpt_ctx.devNetCont);
	cudaFree(cmpt_ctx.p_rnd_states);

}


