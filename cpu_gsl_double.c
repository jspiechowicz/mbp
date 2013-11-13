/*
 * Massive Brownian Particle
 * $\ddot{x} + \gamma\dot{x} = -V'(x) + a\cos(\omega t) + f + \xi(t) + \eta(t)
 * see J. Spiechowicz, J. Luczka and P. Hanggi, J. Stat. Mech. (2013) P02044
 *
 * cpu version with GSL
 *
 * This file is subject to the terms and conditions of the X11 (MIT) License.
 * LM, lukasz.machura@us.edu.pl
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

#define PI 3.14159265358979
#define PI2 6.28318530718
#define true 1
#define false 0

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
    printf("    -a, --amp=DOUBLE         set the AC driving amplitude 'amp' to DOUBLE\n");
    printf("    -b, --omega=DOUBLE       set the AC driving frequency '\\omega' to DOUBLE\n");
    printf("    -c, --force=DOUBLE       set the external bias 'force' to DOUBLE\n");
    printf("    -d, --gam=DOUBLE         set the viscosity '\\gamma' to DOUBLE\n");
    printf("    -e, --Dg=DOUBLE          set the Gaussian noise intensity 'Dg' to DOUBLE\n");
    printf("    -f, --Dp=DOUBLE          set the Poissonian noise intensity 'Dp' to DOUBLE\n");
    printf("    -g, --lambda=DOUBLE      set the Poissonian kicks frequency '\\lambda' to DOUBLE\n\n");
    printf("    -h, --comp=INT          choose between biased and unbiased Poissonian noise. INT can be one of:\n");
    printf("                            0: biased; 1: unbiased\n");
    printf("Simulation params:\n");
    printf("    -k, --paths=LONG        set the number of paths to LONG\n");
    printf("    -l, --periods=LONG      set the number of periods to LONG\n");
    printf("    -m, --trans=DOUBLE      fraction for  transients\n");
    printf("    -n, --spp=INT           specify how many integration steps should be calculated\n");
    printf("                            for a single period of the driving force\n\n");
    printf("    -o, --algorithm=STRING  sets the algorithm. STRING can be one of:\n");
    printf("                            predcorr: simplified weak order 2.0 adapted predictor-corrector\n");
    printf("                            euler: simplified weak order 1.0 regular euler-maruyama\n");
    printf("\n");
}

typedef struct {
  double amp;			//	AC driving amplitude
  double omega;			//	AC driving frequency
  double force;			//	external bias
  double gam;			//	viscosity
  double Dg;			//	Gaussian noise intensity
  double Dp;			//	Poissonian noise intensity
  double lambda;			//	Poissonian kicks frequency
  int biased;			//	choose between biased and unbiased Poissonian noise
  long paths;			//	number of paths
  long periods;			//	number of periods
  double trans;			//	fraction of periods which stands for transients
  long spp;			//	integration step
  int _2ndorder;		//	sets the algorithm
} params;

void set_parameters(int argc, char **argv, params * p){
  int c;
  while( (c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h:k:l:m:n:o", options, NULL)) != EOF) {
    switch (c) {
      case 'a':
	sscanf(optarg, "%lf", &(p->amp));
	break;
      case 'b':
	sscanf(optarg, "%lf", &(p->omega));
	break;
      case 'c':
	sscanf(optarg, "%lf", &(p->force));
	break;
      case 'd':
	sscanf(optarg, "%lf", &(p->gam));
	break;
      case 'e':
	sscanf(optarg, "%lf", &(p->Dg));
	break;
      case 'f':
	sscanf(optarg, "%lf", &(p->Dp));
	break;
      case 'g':
	sscanf(optarg, "%lf", &(p->lambda));
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
	sscanf(optarg, "%lf", &(p->trans));
	break;
      case 'n':
	sscanf(optarg, "%ld", &(p->spp));
	break;
      case 'o':
	p->_2ndorder = !strcmp(optarg, "predcorr") ? true : false;
	break;
    }
  }
}

void dump_params(params * p)
{
  fprintf(stdout,"#\n");
  fprintf(stdout,"#amp:%lf\n",p->amp);
  fprintf(stdout,"#omega:%lf\n",p->omega);
  fprintf(stdout,"#force:%lf\n",p->force);
  fprintf(stdout,"#gamma:%lf\n",p->gam);
  fprintf(stdout,"#Dg:%lf\n",p->Dg);
  fprintf(stdout,"#Dp:%lf\n",p->Dp);
  fprintf(stdout,"#lambda:%lf\n",p->lambda);
  fprintf(stdout,"#%sbiased posisson\n",p->biased?"un":"");
  fprintf(stdout,"#\n");
  fprintf(stdout,"#paths:%ld\n",p->paths);
  fprintf(stdout,"#periods:%ld\n",p->periods);
  fprintf(stdout,"#trans:%lf\n",p->trans);
  fprintf(stdout,"#spp:%ld\n",p->spp);
  fprintf(stdout,"#%s order algorithm\n",p->_2ndorder?"2nd":"1st");
  fprintf(stdout,"#\n");
  fflush(stdout);
}

double drift(double x, double v, double w, params *p)
{
    return -(p->gam)*v - PI2*cos(PI2*x) + (p->amp)*cos(w) + p->force;
}

double diffusion(double dt, params *p, gsl_rng *rg)
{
  double ret = 0.0;

  if (p->Dg != 0.0){
    double r = gsl_ran_flat(rg,0,1);
    if (p->_2ndorder){
      ret = sqrt(6.0*(p->gam)*(p->Dg)*dt);
      if (r <= 1.0/6) 
	ret *= -1;
      else if (r > 2.0/6) 
	ret = 0.0;
    }
    else{
      ret = sqrt(2.0*(p->gam)*(p->Dg)*dt);
      if (r <= 0.5)
	ret *= -1;
    }
  } 

  return ret;
}

double adapted_jump(int *pcd, double dt, params *p, gsl_rng *rg)
{
  double ret = 0.0;
  
  if (p->Dp != 0.0) {
      double comp = sqrt((p->Dp)*(p->lambda))*dt;
      if (*pcd <= 0) {
	double ampmean = sqrt((p->lambda)/(p->Dp));
        double r = gsl_ran_flat(rg,0,1);

	//magic...
	//npcd = (int) floor( -log( curand_uniform(l_state) )/lambda/dt + 0.5 );
	*pcd = (int)roundf(-log(r)/(p->lambda)/dt);
	
	if (p->biased)
	  ret = -log(r)/ampmean - comp;
	else
	  ret = -log(r)/ampmean;
      } 
      else{
	*pcd--;
	if (p->biased) 
	  ret = -comp;
	else 
	  ret = 0.0;
      }
    } 

  return ret;
}

double regular_jump(double dt, params *p, gsl_rng *rg)
{
  double ret = 0.0;

  if (p->Dp != 0.0) {
    double mu = (p->lambda)*dt;
    double ampmean = sqrt((p->lambda)/(p->Dp));
    double comp = sqrt((p->Dp)*(p->lambda))*dt;
    unsigned int n = gsl_ran_poisson(rg, mu);
    
    double s = 0.0;
    int i;
    for (i = 0; i < n; ++i) 
      s -= log(gsl_ran_flat(rg,0,1))/ampmean;
    
    if (p->biased) 
      s -= comp;
    
    ret = s;
  } 

  return ret;
}

void fold(double *x, double base)
//reduce periodic variable to the base domain
{
  *x -= floor(*x/base)*base;
}

void predcorr(double *x, double *v, double *w, int *ix, int *iw, int *pcd, double dt, params *p, gsl_rng *rg)
/* simplified weak order 2.0 adapted predictor-corrector scheme
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 503, p. 532 )
*/
{
  int i;
  double l_xt, l_xtt, l_vt, l_vtt, l_wt, l_wtt, predl_x, predl_v, predl_w, l_x, l_v, l_w;
  double basex = 1.0,
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

    predl_x = l_x + 0.5*(l_xt + l_xtt)*dt;
    predl_v = l_v + 0.5*(l_vt + l_vtt)*dt + diffusion(dt, p, rg);
    predl_w = l_w + 0.5*(l_wt + l_wtt)*dt;

    l_xtt = predl_v;
    l_vtt = drift(predl_x, predl_v, predl_w, p);
    //l_wtt = p->omega;

    l_x += 0.5*(l_xt + l_xtt)*dt;
    l_v += 0.5*(l_vt + l_vtt)*dt + diffusion(dt, p, rg) + adapted_jump(pcd, dt, p, rg);
    l_w += 0.5*(l_wt + l_wtt)*dt;

    //fold path parameters
    //if (fabs(l_x) > basex){
    //  l_x -= floor(l_x/basex)*basex;
    //  ix[i]++;
   // }
    
    if (l_w > basew){
      l_w -= floor(l_w/basew)*basew;
      iw[i]++;
    }

    x[i] = l_x;
    v[i] = l_v;
    w[i] = l_w;
  }

}

void eulermaruyama(double *x, double *v, double *w, int *ix, int *iw, double dt, params *p, gsl_rng *rg)
/* simplified weak order 1.0 regular euler-maruyama scheme 
( see E. Platen, N. Bruti-Liberati; Numerical Solution of Stochastic Differential Equations with Jumps in Finance; Springer 2010; p. 508, 
  C. Kim, E. Lee, P. Talkner, and P.Hanggi; Phys. Rev. E 76; 011109; 2007 ) 
*/ 
{
  int i;
  double l_xt, l_vt, l_wt, l_x, l_v, l_w;
  double basex = 1.0,
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
    /*if (fabs(l_xt) > basex){
      l_xt -= floor(l_xt/basex)*basex;
      ix[i]++;
    }*/
     
    if (l_wt > basew){
      l_wt -= floor(l_wt/basew)*basew;
      iw[i]++;
    }
    
    x[i] = l_xt;
    v[i] = l_vt;
    w[i] = l_wt;
  }
}

void ic(double *x, double a, double b, params *p, gsl_rng *rg)
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
  long i, j,
       steps = (p->spp) * (p->periods),
       trigger = (p->trans) * steps;
  double dt = PI2/(p->omega)/(p->spp); 

  //counters for folding
  size_t sizei = p->paths*sizeof(int);
  int *ix, *iw;
  ix = (int*)malloc(sizei);
  memset(ix,0,sizei);
  iw = (int*)malloc(sizei);
  memset(iw,0,sizei);

  // vectors allocation and initial conditions
  size_t sized = p->paths*sizeof(double);
  double *x, *v, *w;

  x = (double*)malloc(sized);
  ic(x,0,1,p,rg);

  v = (double*)malloc(sized);
  ic(v,-2,2,p,rg);

  w = (double*)malloc(sized);
  ic(w,0,PI2,p,rg);

  int pcd;
  if (p->_2ndorder) //jump countdown
    pcd = floor(-log(gsl_ran_flat(rg,0,1))/(p->lambda)/dt + 0.5);
    
  double sv = 0.0, 
	sv2 = 0.0;
  for (i = 0; i < steps; i++) {

    //algorithm
    if (p->_2ndorder)
      predcorr(x, v, w, ix, iw, &pcd, dt, p, rg);
    else
      eulermaruyama(x, v, w, ix, iw, dt, p, rg);
    
    if (i > trigger) {
      for (j = 0; j < p->paths; ++j){
	sv += v[j];
	sv2 += v[j]*v[j];
      }
    }
  }

  sv  /= (steps-trigger)*(p->paths);
  sv2 /= (steps-trigger)*(p->paths);
  fprintf(stdout,"%e %e\n",sv,sv2);
  fflush(stdout);

  return 1;//EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
  int ret;
  params p = {	 	/*<<<	___default___	*/
    4.2,		/*	amp		*/
    4.9, 		/*	omega		*/
    0.0, 		/*	force		*/
    0.9,  		/*	gamma		*/
    0.001, 		/*	Dg		*/
    0.0,		/*	Dp		*/
    0.0,		/*	lambda		*/
    false,		/*	biased		*/
    33,			/*	paths		*/
    5000,		/*	periods		*/
    0.1,		/*	trans		*/
    500,		/*	spp		*/
    true,		/*	2nd order	*/
  };/*>>>*/	

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

