/*
 * Massive Brownian Particle
 *
 * $\ddot{x} + \gamma\dot{x} = -V'(x) + a\cos(\omega t) + f + \xi(t) + \eta(t)
 *
 * see J. Spiechowicz, J. Luczka and P. Hanggi, J. Stat. Mech. (2013) P02044
 *
 */

/******
 * !!!!!!!!!!!!!!!!!!!!!!!!!
 * NOT FINISHED!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!
 */

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979f

//model
float h_omega;

//simulation
float h_trans;
int h_spp;
long h_paths, h_periods, h_steps, h_trigger;

//vector
float *h_x, *h_v, *h_w;
size_t size_f;

static struct option options[] = {
    {"amp", required_argument, NULL, 'a'},
    {"omega", required_argument, NULL, 'b'},
    {"force", required_argument, NULL, 'c'},
    {"gam", required_argument, NULL, 'd'},
    {"Dg", required_argument, NULL, 'e'},
    {"Dp", required_argument, NULL, 'f'},
    {"lambda", required_argument, NULL, 'g'},
    {"comp", required_argument, NULL, 'h'},
    {"paths", required_argument, NULL, 'k'},
    {"periods", required_argument, NULL, 'l'},
    {"trans", required_argument, NULL, 'm'},
    {"spp", required_argument, NULL, 'n'},
    {"algorithm", required_argument, NULL, 'o'},
};

void usage(char **argv)
{
    printf("Usage: %s <params> \n\n", argv[0]);
    printf("Model params:\n");
    printf("    -a, --amp=FLOAT         set the AC driving amplitude 'amp' to FLOAT\n");
    printf("    -b, --omega=FLOAT       set the AC driving frequency '\\omega' to FLOAT\n");
    printf("    -c, --force=FLOAT       set the external bias 'force' to FLOAT\n");
    printf("    -d, --gam=FLOAT         set the viscosity '\\gamma' to FLOAT\n");
    printf("    -e, --Dg=FLOAT          set the Gaussian noise intensity 'Dg' to FLOAT\n");
    printf("    -f, --Dp=FLOAT          set the Poissonian noise intensity 'Dp' to FLOAT\n");
    printf("    -g, --lambda=FLOAT      set the Poissonian kicks frequency '\\lambda' to FLOAT\n\n");
    printf("    -h, --comp=INT          choose between biased and unbiased Poissonian noise. INT can be one of:\n");
    printf("                            0: biased; 1: unbiased\n");
    printf("Simulation params:\n");
    printf("    -k, --paths=LONG        set the number of paths to LONG\n");
    printf("    -l, --periods=LONG      set the number of periods to LONG\n");
    printf("    -m, --trans=FLOAT       specify fraction FLOAT of periods which stands for transients\n");
    printf("    -n, --spp=INT           specify how many integration steps should be calculated\n");
    printf("                            for a single period of the driving force\n\n");
    printf("    -o, --algorithm=STRING  sets the algorithm. STRING can be one of:\n");
    printf("                            predcorr: simplified weak order 2.0 adapted predictor-corrector\n");
    printf("                            euler: simplified weak order 1.0 regular euler-maruyama\n");
    printf("\n");
}

void parse_cla(int argc, char **argv)
{
    float ftmp;
    int c, itmp;

    while( (c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:y:z:A", options, NULL)) != EOF) {
        switch (c) {
            case 'a':
                ftmp = atof(optarg);
                break;
            case 'b':
                h_omega = atof(optarg);
                break;
            case 'c':
                ftmp = atof(optarg);
                break;
            case 'd':
                ftmp = atof(optarg);
                break;
            case 'e':
                ftmp = atof(optarg);
                break;
            case 'f':
                ftmp = atof(optarg);
                break;
            case 'g':
                ftmp = atof(optarg);
                break;
            case 'h':
                itmp = atoi(optarg);
                break;
            case 'k':
                h_paths = atol(optarg);
                break;
            case 'l':
                h_periods = atol(optarg);
                break;
            case 'm':
                h_trans = atof(optarg);
                break;
            case 'n':
                h_spp = atoi(optarg);
                break;
            case 'o':
                if ( !strcmp(optarg, "predcorr") )
                    itmp = 1;
                else if ( !strcmp(optarg, "euler") )
                    itmp = 0;
                break;
        }
    }
}


float drift(float l_x, float l_v, float l_w, float l_gam, float l_amp, float l_force)
{
    return -l_gam*l_v - 2.0f*PI*cosf(2.0f*PI*l_x) + l_amp*cosf(l_w) + l_force;
}

float diffusion(float l_gam, float l_Dg, float l_dt, int l_2ndorder, curandState *l_state)
{
    if (l_Dg != 0.0f) {
        float r = curand_uniform(l_state);
        if (l_2ndorder) {
            if ( r <= 1.0f/6 ) {
                return -sqrtf(6.0f*l_gam*l_Dg*l_dt);
            } else if ( r > 1.0f/6 && r <= 2.0f/6 ) {
                return sqrtf(6.0f*l_gam*l_Dg*l_dt);
            } else {
                return 0.0f;
            }
        } else {
            if ( r <= 0.5f ) {
                return -sqrtf(2.0f*l_gam*l_Dg*l_dt);
            } else {
                return sqrtf(2.0f*l_gam*l_Dg*l_dt);
            }
        }
    } else {
        return 0.0f;
    }
}

__device__ float adapted_jump(int &npcd, int pcd, float l_lambda, float l_Dp, int l_comp, float l_dt, curandState *l_state)
{
    if (l_Dp != 0.0f) {
        float comp = sqrtf(l_Dp*l_lambda)*l_dt;
        if (pcd <= 0) {
            float ampmean = sqrtf(l_lambda/l_Dp);
           
            npcd = (int) floor( -logf( curand_uniform(l_state) )/l_lambda/l_dt + 0.5f );

            if (l_comp) {
                return -logf( curand_uniform(l_state) )/ampmean - comp;
            } else {
                return -logf( curand_uniform(l_state) )/ampmean;
            }
        } else {
            npcd = pcd - 1;
            if (l_comp) {
                return -comp;
            } else {
                return 0.0f;
            }
        }
    } else {
        return 0.0f;
    }
}

__device__ float regular_jump(float l_lambda, float l_Dp, int l_comp, float l_dt, curandState *l_state)
{
    if (l_Dp != 0.0f) {
        float mu, ampmean, comp, s;
        int i;
        unsigned int n;

        mu = l_lambda*l_dt;
        ampmean = sqrtf(l_lambda/l_Dp);
        comp = sqrtf(l_Dp*l_lambda)*l_dt;
        n = curand_poisson(l_state, mu);
        s = 0.0f;
            for (i = 0; i < n; i++) {
                s += -logf( curand_uniform(l_state) )/ampmean;
            }
        if (l_comp) s -= comp;
        return s;
    } else {
        return 0.0f;
    }
}

__device__ void predcorr(float &corrl_x, float l_x, float &corrl_v, float l_v, float &corrl_w, float l_w, int &npcd, int pcd, curandState *l_state, \
                         float l_amp, float l_omega, float l_force, float l_gam, float l_Dg, int l_2ndorder, float l_Dp, float l_lambda, int l_comp, float l_dt)
/* simplified weak order 2.0 adapted predictor-corrector scheme
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 503, p. 532 )
*/
{
    float l_xt, l_xtt, l_vt, l_vtt, l_wt, l_wtt, predl_x, predl_v, predl_w;

    l_xt = l_v;
    l_vt = drift(l_x, l_v, l_w, l_gam, l_amp, l_force);
    l_wt = l_omega;

    predl_x = l_x + l_xt*l_dt;
    predl_v = l_v + l_vt*l_dt + diffusion(l_gam, l_Dg, l_dt, l_2ndorder, l_state);
    predl_w = l_w + l_wt*l_dt;

    l_xtt = predl_v;
    l_vtt = drift(predl_x, predl_v, predl_w, l_gam, l_amp, l_force);
    l_wtt = l_omega;

    predl_x = l_x + 0.5f*(l_xt + l_xtt)*l_dt;
    predl_v = l_v + 0.5f*(l_vt + l_vtt)*l_dt + diffusion(l_gam, l_Dg, l_dt, l_2ndorder, l_state);
    predl_w = l_w + 0.5f*(l_wt + l_wtt)*l_dt;

    l_xtt = predl_v;
    l_vtt = drift(predl_x, predl_v, predl_w, l_gam, l_amp, l_force);
    l_wtt = l_omega;

    corrl_x = l_x + 0.5f*(l_xt + l_xtt)*l_dt;
    corrl_v = l_v + 0.5f*(l_vt + l_vtt)*l_dt + diffusion(l_gam, l_Dg, l_dt, l_2ndorder, l_state) + adapted_jump(npcd, pcd, l_lambda, l_Dp, l_comp, l_dt, l_state);
    corrl_w = l_w + 0.5f*(l_wt + l_wtt)*l_dt;
}

__device__ void eulermaruyama(float &nl_x, float l_x, float &nl_v, float l_v, float &nl_w, float l_w, curandState *l_state, \
                         float l_amp, float l_omega, float l_force, float l_gam, float l_Dg, int l_2ndorder, float l_Dp, float l_lambda, int l_comp, float l_dt)
/* simplified weak order 1.0 regular euler-maruyama scheme 
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 508, 
  C. Kim, E. Lee, P. Talkner, and P.Hanggi; Phys. Rev. E 76; 011109; 2007 ) 
*/ 
{
    float l_xt, l_vt, l_wt;

    l_vt = l_v + drift(l_x, l_v, l_w, l_gam, l_amp, l_force)*l_dt + diffusion(l_gam, l_Dg, l_dt, l_2ndorder, l_state) 
               + regular_jump(l_lambda, l_Dp, l_comp, l_dt, l_state);
    l_xt = l_x + l_v*l_dt;
    l_wt = l_w + l_omega*l_dt;

    nl_v = l_vt;
    nl_x = l_xt;
    nl_w = l_wt;
}

void fold(float &nx, float x, float y, float &nfc, float fc)
//reduce periodic variable to the base domain
{
    nx = x - floor(x/y)*y;
    nfc = fc + floor(x/y)*y;
}

void run_moments(float *d_x, float *d_v, float *d_w, float *d_sv, float *d_sv2, float *d_dx, curandState *d_states)
//actual moments kernel
{
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    float l_x, l_v, l_w, l_sv, l_sv2, l_dx; 
    curandState l_state;

    //cache path and model parameters in local variables
    l_x = d_x[idx];
    l_v = d_v[idx];
    l_w = d_w[idx];
    l_sv = d_sv[idx];
    l_sv2 = d_sv2[idx];
    l_state = d_states[idx];

    float l_amp, l_omega, l_force, l_gam, l_Dg, l_Dp, l_lambda;
    int l_comp;

    l_amp = d_amp;
    l_omega = d_omega;
    l_force = d_force;
    l_gam = d_gam;
    l_Dg = d_Dg;
    l_Dp = d_Dp;
    l_lambda = d_lambda;
    l_comp = d_comp;

    //run simulation for multiple values of the system parameters
    long ridx = (idx/d_paths) % d_points;
    l_dx = d_dx[ridx];

    switch(d_domainx) {
        case 'a':
            l_amp = l_dx;
            break;
        case 'w':
            l_omega = l_dx;
            break;
        case 'f':
            l_force = l_dx;
            break;
        case 'g':
            l_gam = l_dx;
            break;
        case 'D':
            l_Dg = l_dx;
            break;
        case 'p':
            l_Dp = l_dx;
            break;
        case 'l':
            l_lambda = l_dx;
            break;
    }

    //step size & number of steps
    float l_dt;
    long l_steps, l_trigger, i;

    l_dt = 2.0f*PI/l_omega/d_spp; 
    l_steps = d_steps;
    l_trigger = d_trigger;

    //counters for folding
    float xfc, wfc;
    
    xfc = 0.0f;
    wfc = 0.0f;

    int l_2ndorder, pcd;

    l_2ndorder = d_2ndorder;

    if (l_2ndorder) {
        //jump countdown
        pcd = (int) floor( -logf( curand_uniform(&l_state) )/l_lambda/l_dt + 0.5f );
    }
    
    for (i = 0; i < l_steps; i++) {

        //algorithm
        if (l_2ndorder) {
            predcorr(l_x, l_x, l_v, l_v, l_w, l_w, pcd, pcd, &l_state, l_amp, l_omega, l_force, l_gam, l_Dg, l_2ndorder, l_Dp, l_lambda, l_comp, l_dt);
        } else {
            eulermaruyama(l_x, l_x, l_v, l_v, l_w, l_w, &l_state, l_amp, l_omega, l_force, l_gam, l_Dg, l_2ndorder, l_Dp, l_lambda, l_comp, l_dt);
        }
        
        //fold path parameters
        if ( fabs(l_x) > 1.0f ) {
            fold(l_x, l_x, 1.0f, xfc, xfc);
        }

        if ( l_w > (2.0f*PI) ) {
            fold(l_w, l_w, (2.0f*PI), wfc, wfc);
        }

        if (i >= l_trigger) {
            l_sv += l_v;
            l_sv2 += l_v*l_v;
        }

    }

    //write back path parameters to the global memory
    d_x[idx] = l_x + xfc;
    d_v[idx] = l_v;
    d_w[idx] = l_w;
    d_sv[idx] = l_sv;
    d_sv2[idx] = l_sv2;
    d_states[idx] = l_state;
}

void prepare()
//prepare simulation
{
    //number of steps
    h_steps = h_periods*h_spp;
     
    //host memory allocation
    size_f = h_paths*sizeof(float);
    h_x = (float*)malloc(size_f);
    h_v = (float*)malloc(size_f);
    h_w = (float*)malloc(size_f);

    //moments specific requirements
    h_trigger = h_steps*h_trans;

    //moments vec
    h_sv = (float*)malloc(size_f);
    h_sv2 = (float*)malloc(size_f);
}

void initial_conditions()
//set initial conditions for path parameters
{
    long i;

    for (i = 0; i < h_paths; ++i) {
        h_x[i] = (float)rand()/RAND_MAX; //x in (0,1]
        h_v[i] = 4.0f*(float)rand()/RAND_MAX - 2.0f; //v in (-2,2]
        h_w[i] = 2.0f*PI*(float)rand()/RAND_MAX; //w in (0,2\pi]
    }

    memset(h_sv, 0, size_f);
    memset(h_sv2, 0, size_f);
}

void moments(float *sv, float *sv2, float *dc)
//calculate the first two moments of <v> and diffusion coefficient
{
    float sx = 0.0f, 
	  sx2 = 0.0f;
    int i;

    for (i = 0; i < h_paths; i++) {
      sv += h_sv[i];
      sv2 += h_sv2[i];
      sx += h_x[i];
      sx2 += h_x[i]*h_x[i];
    }

    sv /= (h_steps - h_trigger)/h_paths;
    sv2 /= (h_steps - h_trigger)/h_paths;
    sx /= h_paths;
    sx2 /= h_paths;
    dc = (sx2 - sx*sx)/(2.0f*h_periods*2.0f*PI/h_omega);
}

void finish()
//free memory
{
    free(h_x);
    free(h_v);
    free(h_w);
    free(h_sv);
    free(h_sv2);
    free(h_dx);
}

int main(int argc, char **argv)
{
    //weak RNG init
    srand(time(NULL));

    prepare();
    initial_conditions();
    
    //asymptotic long time average velocity <<v>>, <<v^2>> and diffusion coefficient
    float av  = 0.0f, 
	  av2 = 0.0f, 
	  dc  = 0.0f;
    int i;
    run_moments(h_x, h_v, h_w);
    moments(&av,&av2,&dc);
    
    printf("#<<v>> <<v^2>> D_x\n");
    printf("%e %e %e\n", av, av2, dc);

    finish();
    return 0;
}
