#include "CuMosquitoHutTrialsITN.h"

#include <mex.h>


#pragma comment(lib,"CuMosquitoHutTrialsITN")
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
#define OUT_CTOT    5
#define OUT_NET_CONT 6

#define EXP_TIME 18000

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

	if (nlhs > 7)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidLhsParam", "Number of output parameters can not exceed 7");

	if (!mxIsScalar(prhs[IN_REPETITIONS]) || !mxIsDouble(prhs[IN_REPETITIONS]))
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Number of repetitions (1st input argument) must be a double scalar");

	if (!mxIsScalar(prhs[IN_EXPERIMENTS]) || !mxIsDouble(prhs[IN_REPETITIONS]))
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Number of experiments (2nd input argument) must be a double scalar");

	bool exp_params_validness_test = true;

	mxArray *xlim = 0, *ylim = 0, *sig_acc = 0, *pnet = 0, *phut = 0, *eps = 0, *mu = 0, *tmax = 0, 
		    *d50 = 0, *r = 0, *s = 0, *alpha_p = 0, *alpha_d = 0, *d50_NetCont = 0, *s_NetCont;

	if (!mxIsStruct(prhs[IN_EXPPARAMS])) 
    {
        exp_params_validness_test = false;
        mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid 70");
    }
	else
	{
        
		int exp_params_xlim, exp_params_ylim, exp_params_sig_acc,
			exp_params_pnet, exp_params_phut, exp_params_eps, exp_params_mu, exp_params_tmax,
			exp_params_d50, exp_params_r, exp_params_s, exp_params_alpha_p, exp_params_alpha_d, exp_params_d50_NetCont, exp_params_s_NetCont;

// 		if ((exp_params_xlim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "xlim")) == -1 ||
// 			(exp_params_ylim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "ylim")) == -1 ||
// 			(exp_params_sig_acc = mxGetFieldNumber(prhs[IN_EXPPARAMS], "sig_acc")) == -1 ||
// 			(exp_params_pnet = mxGetFieldNumber(prhs[IN_EXPPARAMS], "pnet")) == -1 ||
// 			(exp_params_phut = mxGetFieldNumber(prhs[IN_EXPPARAMS], "phut")) == -1 ||
// 			(exp_params_eps = mxGetFieldNumber(prhs[IN_EXPPARAMS], "eps")) == -1 ||
// 			(exp_params_mu = mxGetFieldNumber(prhs[IN_EXPPARAMS], "mu")) == -1 ||
// 			(exp_params_tmax = mxGetFieldNumber(prhs[IN_EXPPARAMS], "tmax")) == -1 ||
// 			(exp_params_d50 = mxGetFieldNumber(prhs[IN_EXPPARAMS], "d50")) == -1 || 
// 			(exp_params_s = mxGetFieldNumber(prhs[IN_EXPPARAMS], "s")) == -1 ||
// 			(exp_params_alpha_p = mxGetFieldNumber(prhs[IN_EXPPARAMS], "alpha_p")) == -1 ||
// 			(exp_params_Cmax = mxGetFieldNumber(prhs[IN_EXPPARAMS], "Cmax")) == -1)
// 			exp_params_validness_test = false;
//         
        if ((exp_params_xlim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "xlim")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid xlim");
        if ((exp_params_ylim = mxGetFieldNumber(prhs[IN_EXPPARAMS], "ylim")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid ylim");
        if ((exp_params_sig_acc = mxGetFieldNumber(prhs[IN_EXPPARAMS], "sig_acc")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid xlim");
        if ((exp_params_pnet = mxGetFieldNumber(prhs[IN_EXPPARAMS], "pnet")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid pnet");
        if ((exp_params_phut = mxGetFieldNumber(prhs[IN_EXPPARAMS], "phut")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid phut");
        if ((exp_params_eps = mxGetFieldNumber(prhs[IN_EXPPARAMS], "eps")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid eps");
        if ((exp_params_mu = mxGetFieldNumber(prhs[IN_EXPPARAMS], "mu")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid mu");
        if ((exp_params_tmax = mxGetFieldNumber(prhs[IN_EXPPARAMS], "tmax")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid tmax");
        if ((exp_params_d50 = mxGetFieldNumber(prhs[IN_EXPPARAMS], "d50")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid d50");
		if ((exp_params_r = mxGetFieldNumber(prhs[IN_EXPPARAMS], "r")) == -1)
			mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid r");
        if ((exp_params_s = mxGetFieldNumber(prhs[IN_EXPPARAMS], "s")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid s");
        if ((exp_params_alpha_p = mxGetFieldNumber(prhs[IN_EXPPARAMS], "alpha_p")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid alpha_p");
        if ((exp_params_alpha_d = mxGetFieldNumber(prhs[IN_EXPPARAMS], "alpha_d")) == -1)         
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid alpha_d");
        if ((exp_params_d50_NetCont = mxGetFieldNumber(prhs[IN_EXPPARAMS], "d50_NetCont")) == -1)
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid d50_NetCont");
        if ((exp_params_s_NetCont = mxGetFieldNumber(prhs[IN_EXPPARAMS], "s_NetCont")) == -1)
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid s_NetCont");

		xlim = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_xlim);
		ylim = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_ylim);
		sig_acc = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_sig_acc);
		pnet = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_pnet);
		phut = mxGetFieldByNumber(prhs[IN_EXPPARAMS],0,exp_params_phut);
		eps = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_eps);
		mu = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_mu);
		tmax = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_tmax);
		d50 = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_d50);
		r = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_r);
		s = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_s);
		alpha_p = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_alpha_p);
        alpha_d = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_alpha_d);
		d50_NetCont = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_d50_NetCont);
        s_NetCont = mxGetFieldByNumber(prhs[IN_EXPPARAMS], 0, exp_params_s_NetCont);

// 		if (!exp_params_validness_test || !mxIsDouble(xlim) || !mxIsDouble(ylim) || !mxIsDouble(sig_acc) ||
// 			!mxIsDouble(pnet)|| !mxIsDouble(phut)  || !mxIsDouble(eps) || !mxIsDouble(mu) || !mxIsDouble(tmax) ||
// 			!mxIsDouble(d50) || !mxIsDouble(s) || !mxIsDouble(alpha_p) || !mxIsDouble(Cmax) )
// 			exp_params_validness_test = false;
        if (!mxIsDouble(xlim))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "xlim is not of type double ");
        if (!mxIsDouble(ylim))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "xlim is not of type double ");
        if (!mxIsDouble(sig_acc))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "sig_acc is not of type double ");
        if (!mxIsDouble(pnet))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "pnet is not of type double ");
        if (!mxIsDouble(phut))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "phut is not of type double ");
        if (!mxIsDouble(eps))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "eps is not of type double ");
        if (!mxIsDouble(mu))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "mu is not of type double ");
        if (!mxIsDouble(tmax))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "tmax is not of type double ");
        if (!mxIsDouble(d50))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "d50 is not of type double ");
		if (!mxIsDouble(r))
			mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "r is not of type double ");
        if (!mxIsDouble(s))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "s is not of type double "); 
        if (!mxIsDouble(alpha_p))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "alpha_p is not of type double ");
        if (!mxIsDouble(alpha_d))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "alpha_d is not of type double ");
        if (!mxIsDouble(d50_NetCont))
             mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "d50_NetCont is not of type double ");
        if (!mxIsDouble(s_NetCont))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "s_NetCont is not of type double ");
// 		if (!exp_params_validness_test || !mxIsVector(xlim) || mxGetNumberOfElements(xlim) != 2 ||
// 			!mxIsVector(ylim) || mxGetNumberOfElements(ylim) != 2 ||
// 			!mxIsVector(sig_acc) || mxGetNumberOfElements(sig_acc) != 2 ||
// 			!mxIsScalar(pnet) || !mxIsScalar(phut) || !mxIsScalar(eps) || !mxIsScalar(mu) || !mxIsScalar(tmax) ||
// 			!mxIsScalar(d50) || !mxIsScalar(s) || !mxIsScalar(alpha_p) || mxIsScalar(Cmax) )
// 			exp_params_validness_test = false;
        
        if (!mxIsVector(xlim) || mxGetNumberOfElements(xlim) != 2)
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "xlim is not a vector ");
        if (!mxIsVector(ylim) || mxGetNumberOfElements(ylim) != 2)
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "ylim is not a vector ");
        if (!mxIsVector(sig_acc) || mxGetNumberOfElements(sig_acc) != 2)
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "sig_acc is not a vector ");
        if (!mxIsScalar(pnet))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "pnet is not a scalar ");
        if (!mxIsScalar(phut))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "phut is not a scalar ");
        if (!mxIsScalar(eps))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "eps is not a scalar ");
        if (!mxIsScalar(mu))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "mu is not a scalar ");
        if (!mxIsScalar(tmax))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "tmax is not a scalar ");
        if (!mxIsScalar(d50))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "d50 is not a scalar ");
		if (!mxIsScalar(r))
			mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "r is not a scalar ");
        if (!mxIsScalar(s))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "s is not a scalar ");
        if (!mxIsScalar(alpha_p))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "alpha_p is not a scalar ");
        if (!mxIsScalar(alpha_d))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "alpha_d is not a scalar ");
        if (!mxIsScalar(d50_NetCont))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "d50_NetCont is not a scalar ");
        if (!mxIsScalar(s_NetCont))
            mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "s_NetCont is not a scalar ");
//        mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid 125");
	}

	if (!exp_params_validness_test)
		mexErrMsgIdAndTxt("CuMosquitoRepellent:InvalidRhsParam", "Experimental parameters descriptor (3rd argument) is invalid");

	
	void* p_raw_data = 0;
	uint32_t new_repetitions, new_experiments;
	hut_exp_ITN exp_params;

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

	//Initialize d50
	p_raw_data = mxGetData(d50);
	exp_params.d50 = *((double*)p_raw_data);

	//Initialize d50
	p_raw_data = mxGetData(r);
	exp_params.r = *((double*)p_raw_data);

	//Initialize s
	p_raw_data = mxGetData(s);
	exp_params.s = *((double*)p_raw_data);

	//Initialize alpha_p
	p_raw_data = mxGetData(alpha_p);
	exp_params.alpha_p = *((double*)p_raw_data);
    
    //Initialize alpha_d
	p_raw_data = mxGetData(alpha_d);
	exp_params.alpha_d = *((double*)p_raw_data);

	//Initialize d50_NetCont
	p_raw_data = mxGetData(d50_NetCont);
	exp_params.d50_NetCont = *((double*)p_raw_data);

    //Initialize d50_NetCont
    p_raw_data = mxGetData(s_NetCont);
    exp_params.s_NetCont = *((double*)p_raw_data);

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
   //plhs[OUT_CTOT] = mxCreateLogicalMatrix(experiments, repetitions);
    //plhs[OUT_CTOT] = mxCreateDoubleMatrix(experiments, repetitions, mxREAL);
    size_t CTOT_dims[] = {experiments,repetitions,EXP_TIME};
    plhs[OUT_CTOT] = mxCreateUninitNumericArray(3, CTOT_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[OUT_NET_CONT] = mxCreateUninitNumericArray(3, CTOT_dims, mxDOUBLE_CLASS, mxREAL);
	
	memcpy(mxGetData(plhs[OUT_IN]), computation_context.hostIn, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_DEAD]), computation_context.hostDead, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_TRAP]), computation_context.hostTrap, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_FED]), computation_context.hostFed, experiments*repetitions);
	memcpy(mxGetData(plhs[OUT_UNFDEAD]), computation_context.hostUnfDead, experiments*repetitions);
    memcpy(mxGetData(plhs[OUT_CTOT]), computation_context.hostCtot, experiments*repetitions*EXP_TIME*sizeof(double));
    memcpy(mxGetData(plhs[OUT_NET_CONT]), computation_context.hostNetCont, experiments*repetitions*EXP_TIME*sizeof(double));
 
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