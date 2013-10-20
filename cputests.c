/*
 * Massive Brownian Particle
 * $\ddot{x} + \gamma\dot{x} = -V'(x) + a\cos(\omega t) + f + \xi(t) + \eta(t)
 * see J. Spiechowicz, J. Luczka and P. Hanggi, J. Stat. Mech. (2013) P02044
 *
 * cpu version with GSL
 *
 * this is a testing version, 
 * with prints after 10^N periods (N=1,2,3,4,5...)
 * 
 */
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define PI 3.14159265358979f
#define PI2 6.28318530718f
#define true 1
#define false 0

clock_t tic;

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
    {"help", no_argument, NULL, 'p'},
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
    printf("                            0: biased (default); 1: unbiased\n");
    printf("Simulation params:\n");
    printf("    -k, --paths=LONG        set the number of paths to LONG\n");
    printf("    -l, --periods=LONG      set the number of periods to LONG\n");
    printf("    -m, --trans=LONG        set the number of periods which stands for transients to LONG\n");
    printf("    -n, --spp=INT           specify how many integration steps should be calculated\n");
    printf("                            for a single period of the driving force\n\n");
    printf("    -o, --algorithm=STRING  sets the algorithm. STRING can be one of:\n");
    printf("                            predcorr: simplified weak order 2.0 adapted predictor-corrector\n");
    printf("                            euler: simplified weak order 1.0 regular euler-maruyama\n");
    printf("    -p, --help		prints this help and exits\n");
    printf("\n");
}

typedef struct {
  float amp;			//	AC driving amplitude
  float omega;			//	AC driving frequency
  float force;			//	external bias
  float gam;			//	viscosity
  float Dg;			//	Gaussian noise intensity
  float Dp;			//	Poissonian noise intensity
  float lambda;			//	Poissonian kicks frequency
  int biased;			//	choose between biased and unbiased Poissonian noise
  long paths;			//	number of paths
  long periods;			//	number of periods
  long trans;			//	fraction of periods which stands for transients
  long spp;			//	integration step
  int _2ndorder;		//	sets the algorithm
} params;

void set_parameters(int argc, char **argv, params * p){
  int c;
  while( (c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h:k:l:m:n:o:p", options, NULL)) != EOF) {
    switch (c) {
      case 'a':
	sscanf(optarg, "%f", &(p->amp));
	break;
      case 'b':
	sscanf(optarg, "%f", &(p->omega));
	break;
      case 'c':
	sscanf(optarg, "%f", &(p->force));
	break;
      case 'd':
	sscanf(optarg, "%f", &(p->gam));
	break;
      case 'e':
	sscanf(optarg, "%f", &(p->Dg));
	break;
      case 'f':
	sscanf(optarg, "%f", &(p->Dp));
	break;
      case 'g':
	sscanf(optarg, "%f", &(p->lambda));
	break;
      case 'h':
	sscanf(optarg, "%d", &(p->biased));
	break;
      case 'k':
	sscanf(optarg, "%ld", &(p->paths));
	break;
      case 'l':
	sscanf(optarg, "%ld", &(p->periods));
	break;
      case 'm':
	sscanf(optarg, "%ld", &(p->trans));
	break;
      case 'n':
	sscanf(optarg, "%ld", &(p->spp));
	break;
      case 'o':
	p->_2ndorder = !strcmp(optarg, "predcorr") ? true : false;
	break;
      case 'p':
	usage(argv);
	exit(0);
	break;
    }
  }
}

void dump_params(params * p)
{
  fprintf(stdout,"#\n");
  fprintf(stdout,"#amp:%f\n",p->amp);
  fprintf(stdout,"#omega:%f\n",p->omega);
  fprintf(stdout,"#force:%f\n",p->force);
  fprintf(stdout,"#gamma:%f\n",p->gam);
  fprintf(stdout,"#Dg:%f\n",p->Dg);
  fprintf(stdout,"#Dp:%f\n",p->Dp);
  fprintf(stdout,"#lambda:%f\n",p->lambda);
  fprintf(stdout,"#%sbiased posisson\n",p->biased?"un":"");
  fprintf(stdout,"#\n");
  fprintf(stdout,"#paths:%ld\n",p->paths);
  fprintf(stdout,"#periods:%ld\n",p->periods);
  fprintf(stdout,"#trans:%ld\n",p->trans);
  fprintf(stdout,"#spp:%ld\n",p->spp);
  fprintf(stdout,"#%s order algorithm\n",p->_2ndorder?"2nd":"1st");
  fprintf(stdout,"#\n");
  fflush(stdout);
}

float drift(float x, float v, float w, params *p)
{
    return -(p->gam)*v - PI2*cosf(PI2*x) + (p->amp)*cosf(w) + p->force;
}

float diffusion(float dt, params *p, gsl_rng *rg)
{
  float ret = 0.0f;

  if (p->Dg != 0.0f){
    float r = gsl_ran_flat(rg,0,1);
    if (p->_2ndorder){
      ret = sqrtf(6.0f*(p->gam)*(p->Dg)*dt);
      if (r <= 1.0f/6) 
	ret *= -1;
      else if (r > 2.0f/6) 
	ret = 0.0f;
    }
    else{
      ret = sqrtf(2.0f*(p->gam)*(p->Dg)*dt);
      if (r <= 0.5f)
	ret *= -1;
    }
  } 

  return ret;
}

float adapted_jump(int *pcd, float dt, params *p, gsl_rng *rg)
{
  float ret = 0.0f;
  
  if (p->Dp != 0.0f) {
      float comp = sqrtf((p->Dp)*(p->lambda))*dt;
      if (*pcd <= 0) {
	float ampmean = sqrtf((p->lambda)/(p->Dp));
        float r = gsl_ran_flat(rg,0,1);

	//magic...
	//npcd = (int) floor( -logf( curand_uniform(l_state) )/lambda/dt + 0.5f );
	*pcd = (int)roundf(-logf(r)/(p->lambda)/dt);
	
	if (p->biased)
	  ret = -logf(r)/ampmean - comp;
	else
	  ret = -logf(r)/ampmean;
      } 
      else{
	*pcd--;
	if (p->biased) 
	  ret = -comp;
	else 
	  ret = 0.0f;
      }
    } 

  return ret;
}

float regular_jump(float dt, params *p, gsl_rng *rg)
{
  float ret = 0.0f;

  if (p->Dp != 0.0f) {
    float mu = (p->lambda)*dt;
    float ampmean = sqrtf((p->lambda)/(p->Dp));
    float comp = sqrtf((p->Dp)*(p->lambda))*dt;
    unsigned int n = gsl_ran_poisson(rg, mu);
    
    float s = 0.0f;
    int i;
    for (i = 0; i < n; ++i) 
      s -= logf(gsl_ran_flat(rg,0,1))/ampmean;
    
    if (p->biased) 
      s -= comp;
    
    ret = s;
  } 

  return ret;
}

void fold(float *x, float base)
//reduce periodic variable to the base domain
{
  *x -= floor(*x/base)*base;
}

void predcorr(float *x, float *v, float *w, int *ix, int *iw, int *pcd, float dt, params *p, gsl_rng *rg)
/* simplified weak order 2.0 adapted predictor-corrector scheme
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 503, p. 532 )
*/
{
  int i;
  float l_xt, l_xtt, l_vt, l_vtt, l_wt, l_wtt, predl_x, predl_v, predl_w, l_x, l_v, l_w;
  float basex = 1.0f,
	basew = PI2;

  for (i = 0; i < p->paths; ++i){
    l_x = x[i];
    l_v = v[i];
    l_w = w[i];

    l_xt = l_v;
    l_vt = drift(l_x, l_v, l_w, p);
    l_wt = p->omega;

    predl_x = l_x + l_xt*dt;
    predl_v = l_v + l_vt*dt + diffusion(dt, p, rg);
    predl_w = l_w + l_wt*dt;

    l_xtt = predl_v;
    l_vtt = drift(predl_x, predl_v, predl_w, p);
    l_wtt = p->omega;

    predl_x = l_x + 0.5f*(l_xt + l_xtt)*dt;
    predl_v = l_v + 0.5f*(l_vt + l_vtt)*dt + diffusion(dt, p, rg);
    predl_w = l_w + 0.5f*(l_wt + l_wtt)*dt;

    l_xtt = predl_v;
    l_vtt = drift(predl_x, predl_v, predl_w, p);
    //l_wtt = p->omega;

    l_x += 0.5f*(l_xt + l_xtt)*dt;
    l_v += 0.5f*(l_vt + l_vtt)*dt + diffusion(dt, p, rg) + adapted_jump(pcd, dt, p, rg);
    l_w += 0.5f*(l_wt + l_wtt)*dt;

    //fold path parameters
    if (fabs(l_x) > basex){
      l_x -= floor(l_x/basex)*basex;
      ix[i]++;
    }
    
    if (l_w > basew){
      l_w -= floor(l_w/basew)*basew;
      iw[i]++;
    }

    x[i] = l_x;
    v[i] = l_v;
    w[i] = l_w;
  }

}

void eulermaruyama(float *x, float *v, float *w, int *ix, int *iw, float dt, params *p, gsl_rng *rg)
/* simplified weak order 1.0 regular euler-maruyama scheme 
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 508, 
  C. Kim, E. Lee, P. Talkner, and P.Hanggi; Phys. Rev. E 76; 011109; 2007 ) 
*/ 
{
  int i;
  float l_xt, l_vt, l_wt, l_x, l_v, l_w;
  float basex = 1.0f,
	basew = PI2;

  for (i = 0; i < p->paths; ++i){
    l_x = x[i];
    l_v = v[i];
    l_w = w[i];

    l_xt = l_x + l_v*dt;
    l_vt = l_v + drift(l_x, l_v, l_w, p)*dt 
      + diffusion(dt, p, rg) 
      + regular_jump(dt, p, rg);
    l_wt = l_w + (p->omega)*dt;

    //fold path parameters
    if (fabs(l_xt) > basex){
      l_xt -= floor(l_xt/basex)*basex;
      ix[i]++;
    }
     
    if (l_wt > basew){
      l_wt -= floor(l_wt/basew)*basew;
      iw[i]++;
    }
    
    x[i] = l_xt;
    v[i] = l_vt;
    w[i] = l_wt;
  }
}

void ic(float *x, float a, float b, params *p, gsl_rng *rg)
//initial conditions
{
  int i;
  for (i = 0; i < p->paths; ++i)
    x[i] = gsl_ran_flat(rg, a, b);
}

int simulate(params *p, gsl_rng *rg)
//actual moments kernel
{
  //silent counters
  //step size & number of steps
  long i, j, period, steps,
       total_periods = (p->periods) + (p->trans);
  float dt = PI2/(p->omega)/(p->spp); 

  //counters for folding
  size_t sizei = p->paths*sizeof(int);
  int *ix, *iw;
  ix = (int*)malloc(sizei);
  memset(ix,0,sizei);
  iw = (int*)malloc(sizei);
  memset(iw,0,sizei);

  // vectors allocation and initial conditions
  size_t sizef = p->paths*sizeof(float);
  float *x, *v, *w;

  x = (float*)malloc(sizef);
  ic(x,0,1,p,rg);

  v = (float*)malloc(sizef);
  ic(v,-2,2,p,rg);

  w = (float*)malloc(sizef);
  ic(w,0,PI2,p,rg);

  int pcd;
  if (p->_2ndorder) //jump countdown
    pcd = floor(-logf(gsl_ran_flat(rg,0,1))/(p->lambda)/dt + 0.5f);
    
  float sv = 0.0f, 
	sv2 = 0.0f;
  int factor10 = 10;
  fprintf(stdout,"#period <<v>> <<v^2>> tic-tac\n");
  fflush(stdout);
  for (period = 0; period < total_periods; ++period) {
    for (i = 0; i < p->spp; ++i){
      //algorithm
      if (p->_2ndorder)
	predcorr(x, v, w, ix, iw, &pcd, dt, p, rg);
      else
	eulermaruyama(x, v, w, ix, iw, dt, p, rg);
      
      if (period > p->trans){
	for (j = 0; j < p->paths; ++j){
	  sv += v[j];
	  sv2 += v[j]*v[j];
	}
      }
    }

    if (period == factor10){
      steps = (p->spp) * period;
      factor10 *=10;
      sv  /= steps*(p->paths);
      sv2 /= steps*(p->paths);
      fprintf(stdout,"%ld %e %e %f\n",period,sv,sv2,(clock()-tic)/(double)CLOCKS_PER_SEC);
      fflush(stdout);
    }
  }
  
  return 1;//EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
  int ret;
  params p = {	 	/*<<<	___default___	*/
    4.2f,		/*	amp		*/
    4.9f, 		/*	omega		*/
    0.0f, 		/*	force		*/
    0.9f,  		/*	gamma		*/
    0.001f, 		/*	Dg		*/
    0.0f,		/*	Dp		*/
    0.0f,		/*	lambda		*/
    false,		/*	biased		*/
    10,			/*	paths		*/
    10000,		/*	periods		*/
    100,		/*	trans		*/
    200,		/*	spp		*/
    true,		/*	2nd order	*/
  };/*>>>*/	

  tic = clock();
  set_parameters(argc,argv,&p);

  /*<<<      GSL                     */
  const gsl_rng_type * R;
  gsl_rng * rg;
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED")){
    gsl_rng_default_seed = time(NULL);
    fprintf(stderr,"#seed:%ld\n",gsl_rng_default_seed);
  }
  R = gsl_rng_default;
  rg = gsl_rng_alloc (R);
  /*>>>GSL*/

  dump_params(&p);

  ret = simulate(&p, rg);

  gsl_rng_free(rg);
  return ret;
}

