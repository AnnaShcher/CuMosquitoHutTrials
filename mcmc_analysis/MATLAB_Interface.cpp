#include "CuMosquitoHutTrialsControl.cuh"

#include <mex.h>


#pragma comment(lib,"CuMosquitoHutTrialsControl")
#pragma comment(lib,"libmex")
#pragma comment(lib,"libmx")

#pragma comment(lib,"cudart")


#define IN_REPETITIONS 0
#define IN_EXPERIMENTS 1
#define IN_EXPPARAMS  2

#define OUT_IN      0
#define OUT_DEAD    1
#define OUT_TRAP    2
#define OUT_FED     3
#define OUT_UNFDEAD 4

static uint32_t repetitions = 0, experiments = 0;
static bool is_initialized = false;
static ComputationContext computation_context;

void AtExit(void);

bool mxIsScalar(const mxArray* p_mxArray)
{
	mwSize dim_num = mxGetNumberOfDimensions(p_mxArray);
	const mwSize* dims = mxGetDimensions(p_mxArray);
	if (dim_num > 2 || dims[0] != 1 || dims[1] != 1)
		return false;
	return true;
}

bool mxIsVector(const mxArray* p_mxArray)
{
	mwSize dim_num = mxGetNumberOfDimensions(p_mxArray);
	const mwSize* dims = mxGetDimensions(p_mxArray);

	if (dim_num > 2 || dims[0] > 1 && dims[1] > 1)
		return false;
	return true;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nrhs != 3)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Invalid number of input parameters: must be 3");

	if (nlhs > 5)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidLhsParam", "Number of output parameters can not exceed 5");

	if (!mxIsScalar(prhs[IN_REPETITIONS]) || !mxIsDouble(prhs[IN_REPETITIONS]))
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Number of repetitions (1st input argument) must be a double scalar");

	if (!mxIsScalar(prhs[IN_EXPERIMENTS]) || !mxIsDouble(prhs[IN_REPETITIONS]))
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Number of experiments (2nd input argument) must be a double scalar");

	bool exp_params_validness_test = true;

	mxArray *xlim = 0, *ylim = 0, *sig_acc = 0, *pnet = 0, *phut = 0, *eps = 0, *mu = 0, *tmax = 0;

	if (!mxIsStruct(prhs[IN_EXPPARAMS])) exp_params_validness_test = false;
	else
	{
		if (!mxIsScalar(prhs[IN_EXPPARAMS])) exp_params_validness_test = false;

		int exp_params_xlim, exp_params_ylim, exp_params_sig_acc,
			exp_params_pnet, exp_params_phut, exp_params_eps, exp_params_mu, exp_params_tmax;

		if ((exp_params_xlim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "xlim")) == -1 ||
			(exp_params_ylim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "ylim")) == -1 ||
			(exp_params_sig_acc = mxGetFieldNumber(prhs[IN_EXPPARAMS], "sig_acc")) == -1 ||
			(exp_params_pnet = mxGetFieldNumber(prhs[IN_EXPPARAMS], "pnet")) == -1 ||
			(exp_params_phut = mxGetFieldNumber(prhs[IN_EXPPARAMS], "phut")) == -1 ||
			(exp_params_eps = mxGetFieldNumber(prhs[IN_EXPPARAMS], "eps")) == -1 ||
			(exp_params_mu = mxGetFieldNumber(prhs[IN_EXPPARAMS], "mu")) == -1 ||
			(exp_params_tmax = mxGetFieldNumber(prhs[IN_EXPPARAMS], "tmax")) == -1)
			exp_params_validness_test = false;

		xlim = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_xlim);
		ylim = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_ylim);
		sig_acc = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_sig_acc);
		pnet = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_pnet);
		phut = mxGetFieldByNumber(prhs[IN_EXPPARAMS],0,exp_params_phut);
		eps = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_eps);
		mu = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_mu);
		tmax = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_tmax);

		if (!exp_params_validness_test || !mxIsDouble(xlim) || !mxIsDouble(ylim) || !mxIsDouble(sig_acc) ||
			!mxIsDouble(pnet)|| !mxIsDouble(phut)  || !mxIsDouble(eps) || !mxIsDouble(mu) || !mxIsDouble(tmax))
			exp_params_validness_test = false;

		if (!exp_params_validness_test || !mxIsVector(xlim) || mxGetNumberOfElements(xlim) != 2 ||
			!mxIsVector(ylim) || mxGetNumberOfElements(ylim) != 2 ||
			!mxIsVector(sig_acc) || mxGetNumberOfElements(sig_acc) != 2 ||
			!mxIsScalar(pnet) || !mxIsScalar(phut) || !mxIsScalar(eps) || !mxIsScalar(mu) || !mxIsScalar(tmax))
			exp_params_validness_test = false;

	}

	if (!exp_params_validness_test)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid");

	
	void* p_raw_data = 0;
	uint32_t new_repetitions, new_experiments;
	hut_exp_control exp_params;

	//Initialize repetitions
	p_raw_data = mxGetData(prhs[IN_REPETITIONS]);
	new_repetitions = (uint32_t)*((double*)p_raw_data);

	//Initialize experiments
	p_raw_data = mxGetData(prhs[IN_EXPERIMENTS]);
	new_experiments = (uint32_t)*((double*)p_raw_data);

	//Initialize xlim
	p_raw_data = mxGetData(xlim);
	exp_params.xlim[0] = *((double*)p_raw_data);
	exp_params.xlim[1] = *((double*)p_raw_data + 1);

	//Initialize ylim
	p_raw_data = mxGetData(ylim);
	exp_params.ylim[0] = *((double*)p_raw_data);
	exp_params.ylim[1] = *((double*)p_raw_data + 1);

	//Initialize sig_acc
	p_raw_data = mxGetData(sig_acc);
	exp_params.sig_acc[0] = *((double*)p_raw_data);
	exp_params.sig_acc[1] = *((double*)p_raw_data + 1);

	//Initialize pnet
	p_raw_data = mxGetData(pnet);
	exp_params.pnet = *((double*)p_raw_data);

	//Initialize phut
	p_raw_data = mxGetData(phut);
	exp_params.phut = *((double*)p_raw_data);

	//Initialize eps
	p_raw_data = mxGetData(eps);
	exp_params.eps = *((double*)p_raw_data);

	//Initialize mu
	p_raw_data = mxGetData(mu);
	exp_params.mu = *((double*)p_raw_data);

	//Initialize tmax
	p_raw_data = mxGetData(tmax);
	exp_params.tmax = *((double*)p_raw_data);

	if (new_repetitions == 0 || new_experiments == 0)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidCUDAGridParam", "Dimensions of the CUDA grid must not be 0!");

	if (repetitions != new_repetitions || experiments != new_experiments)
	{
		if (is_initialized)
			shutdown(computation_context);
		else
			is_initialized = true;

		initialize(new_repetitions, new_experiments, computation_context);
		repetitions = new_repetitions;
		experiments = new_experiments;
	}

	CallKernel(repetitions, experiments, exp_params, computation_context);


	//Create output data
	plhs[OUT_IN] = mxCreateLogicalMatrix(experiments, repetitions);
	plhs[OUT_DEAD] = mxCreateLogicalMatrix(experiments, repetitions);
	plhs[OUT_TRAP] = mxCreateLogicalMatrix(experiments, repetitions);
	plhs[OUT_FED] = mxCreateLogicalMatrix(experiments, repetitions);
	plhs[OUT_UNFDEAD] = mxCreateLogicalMatrix(experiments, repetitions);

	
	memcpy(mxGetData(plhs[OUT_IN]), computation_context.hostIn, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_DEAD]), computation_context.hostDead, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_TRAP]), computation_context.hostTrap, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_FED]), computation_context.hostFed, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_UNFDEAD]), computation_context.hostUnfDead, experiments*repetitions);

	//Register exit callback that cleans the thing up...
	mexAtExit(AtExit);
}


void AtExit()
{
	if (is_initialized)
		shutdown(computation_context);
}

/*int main(int argc, char* argv[])
{
	repetitions = 4;
	experiments = 600;

	hut_exp_control exp_params;
	exp_params.pnet = 0.9995;
	exp_params.mu = 0.08;
	exp_params.tmax = 9.0;
	exp_params.sig_acc[0] = 1e-5;
	exp_params.sig_acc[1] = 1e-2;
	exp_params.xlim[0] = -60;
	exp_params.xlim[1] = 60;
	exp_params.ylim[0] = -60;
	exp_params.ylim[1] = 60;
	exp_params.eps = 0.2;
	repetitions = 4;
	experiments = 600;

	initialize(repetitions, experiments, computation_context);
	CallKernel(repetitions, experiments, exp_params, computation_context);
	shutdown(computation_context);

	return EXIT_SUCCESS;
}*/