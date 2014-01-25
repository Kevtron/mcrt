#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define MUBINS 10
#define NPHOTONS 1000000
#define ALBEDO 1.0
#define TAUMAX 10.0

/*
 This program calculates the emergent intensity from a plane parallel slab of optical depth TAUMAX, using NPHOTONS photon packets.
 */
/*==============================================================*/
/*---------------------Struct definitions-----------------------*/
/*==============================================================*/
typedef struct {float cost, phi, sint, cosp, sinp, x, y, z;}photon;
typedef struct {int bins[MUBINS];float theta[MUBINS];}hist;
/*==============================================================*/
/*------------------Function declarations-----------------------*/
/*==============================================================*/
hist *binphotons(float mu[]);
photon *newphoton();
void *scatter(photon *p);
void *intensity(hist *binned, float *p);
/*==============================================================*/
/*------------Main loop, loop over stellar photons--------------*/
/*==============================================================*/
int main(int argc, char *argv[]){
	int i;
	static float mu[NPHOTONS];
	static float intens[MUBINS];
	float *intensp;
	hist workspace;
	//workspace.bins[MUBINS]={0};
	hist *histp;
	static photon test;
	photon *working;
	working = &test;
	histp = &workspace;
	intensp = &intens;
	FILE *file;
	srand(time(NULL));
	/* Initializing works.*/
	for(i=0;i<=NPHOTONS;i++){
		*working = *newphoton();
		*scatter(working);
		mu[i]=working->cost;
		if(i%10000==0){
			printf("%d photons completed\n",i);
		}
	}
	*histp = *binphotons(mu);
	intensity(histp, intensp);
	file=fopen("intensity-c.dat","a+");
	fprintf(file,"%s","#intensity, theta\n");
	for(i=0;i<10;i++){
		fprintf(file,"%f\t%f\n",intens[i],histp->theta[i]);
	}
	fclose(file); /*done!*/ 
	return 0;
}
/*==============================================================*/
/*------------------Function definitions------------------------*/
/*==============================================================*/

/*-----------------------New photons----------------------------*/
photon *newphoton(){
	photon *new;
	photon test;
	new = &test;
	float tmp;
	new->cost = sqrt((float)rand()/RAND_MAX); 				//cos theta
	new->phi = 2.0*M_PI*(float)rand()/RAND_MAX;			// phi
	tmp = new->cost;
	new->sint =  sqrt(1-tmp*tmp); 				//sin theta
	new->cosp = cos(new->phi);						//cos PHI
	new->sinp = sin(new->phi);						//sin phi
	new->x = 0.0;							// x
	new->y = 0.0;							// y
	new->z = 0.0;							// z
	return new;
}

/*-------------------Scatter photons----------------------------*/
void *scatter(photon *scatterptr){
	float tau,s;
	while (scatterptr->z>=0.0 && scatterptr->z<=1.0){
		tau = -1.0*log((float)rand()/RAND_MAX);
		s=tau/TAUMAX;
		scatterptr->x +=s*scatterptr->sint*scatterptr->cosp;			//x
		scatterptr->y +=s*scatterptr->sint*scatterptr->sinp;			//y
		scatterptr->z +=s*scatterptr->cost;					//z
		if (scatterptr->z<0.0){
		*scatterptr = *newphoton();
		}
		if ((float)rand()/RAND_MAX < ALBEDO && scatterptr->z < 1.0){
			scatterptr->cost = 2.0*((float)rand()/RAND_MAX) - 1.0;			//cos theta
			scatterptr->phi = 2.0*M_PI*((float)rand()/RAND_MAX);			//phi
			scatterptr->sint = sqrt(1-(scatterptr->cost)*(scatterptr->cost));			//sin theta
			scatterptr->cosp = cos(scatterptr->phi);					//cos phi
			scatterptr->sinp = sin(scatterptr->phi);					//sin phi
		}
	}	
	return 0;
}

/*-----------------------Bin photons----------------------------*/
hist *binphotons(float mu[]){
	float dtheta,width;
	hist binspace;
	hist *binned;
	binned = &binspace;
	int i,j;
	dtheta = 1.0/MUBINS;
	width = 0.5*dtheta;
	for(i=0;i<MUBINS;i++){
	binned->bins[i]=0;
	binned->theta[i]=acos(i*dtheta+width)*(180.0/M_PI);
	}
	for(i=0;i<NPHOTONS;i++){
		j=abs((int)(mu[i]*MUBINS));
		binned->bins[j]+=1;
	}
	return binned;
}
/*---------------------Calculate intensity-----------------------*/
void *intensity(hist *binned, float *p){
	int i;
	for(i=0;i<=(int)MUBINS;i++){
		p[i] = ((float)binned->bins[i]*(float)MUBINS)/(2.0*(float)NPHOTONS*cos(binned->theta[i]*M_PI/180.0));
	}
	return 0;
}
