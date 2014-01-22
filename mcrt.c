#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi M_PI
#define max RAND_MAX
#define mubins 10.0
#define nphotons 1000000
/*
This program calculates the emergent intensity from
*/

/*==============================================================*/
/*------------------Structure definitions-----------------------*/
/*==============================================================*/
typedef struct{float cost,sint,phi,cosp,sinp,x,y,z;}photon;
typedef struct{int bins[(int)mubins]; float theta[(int)mubins];}hist;

/*==============================================================*/
/*------------------Function declarations-----------------------*/
/*==============================================================*/
float gen_rand(void);
hist binphotons(float mu[nphotons]);
photon newphoton();
photon scatter(float albedo, float taumax);
float intensity();
/*==============================================================*/
/*------------Main loop, loop over stellar photons--------------*/
/*==============================================================*/
int main(int argc, char *argv[]){
	float albedo, taumax, mu[nphotons];
	int a[10],i;
	//pi = M_PI; //pi defined in math.h 
	//max = RAND_MAX; //need rands between 0 and 1
	//photon p={1,1,1,1,1};
	for(i=0;i<10;i++){
	printf("%f\n",(float)rand()/max);
	}
	//printf("%f,%f\n",p.cost,p.x);
	/*for(i=0;i<10;i++){
		a[i]=i*i;
		printf("%d,%d\n",i,a[i]);
		}*/
}
/*==============================================================*/
/*------------------Function definitions------------------------*/
/*==============================================================*/

/*-----------------------New photons----------------------------*/
photon newphoton(){
	photon new;
	new.cost = sqrt((float)rand()/max);
	new.phi = 2.0*pi*(float)rand()/max;
	new.sint =  sqrt(1-new.cost*new.cost);
	new.cosp = cos(new.phi);
	new.sinp = sin(new.phi);
	new.x = 0.0;
	new.y = 0.0;
	new.z = 0.0;
	return new;
}

/*-------------------Scatter photons----------------------------*/
photon scatter(float albedo, float taumax){
	photon scatter;
	float tau,s;
	scatter.cost = 2.0*((float)rand()/max) - 1.0;
	scatter.sint = sqrt(1-scatter.cost*scatter.cost);
	scatter.phi = 2.0*pi*((float)rand()/max);
	scatter.cosp = cos(scatter.phi);
	scatter.sinp = sin(scatter.phi);
	if ((float)rand()/max < albedo && scatter.z > 0.0 && scatter.z < 1.0){
		tau = -1.0*log((float)rand()/max);
		s=tau/taumax;
		scatter.x +=s*scatter.sint*scatter.cosp;
		scatter.y +=s*scatter.sint*scatter.sinp;
		scatter.z +=s*scatter.cost;
	}
	return scatter;
}

/*-----------------------Bin photons----------------------------*/
hist binphotons(float mu[nphotons]){
	float dtheta,width;
	hist binned;
	int i;
	dtheta = 1.0/mubins;
	width = 0.5*dtheta;
	for(i=0;i<(int)mubins;i++){
		binned.theta[i]=acos(i*dtheta+width)*(180.0/pi);
	}
	for(i=0;i<(int)mubins;i++){
		int j;
		j=abs((int)(mu[i]*mubins));
		binned.bins[j]+=1;
	}
	return binned;
}
/*---------------------Calculate intensity-----------------------*/
float intensity(hist binned){
	int i;
	float intens[(int)mubins];
	for(i=0;i<(int)mubins;i++){
	intens[i] = ((float)binned.bins[i]*(float)mubins)/(2.0*(float)nphotons*cos(binned.theta[i]*(float)pi/180.0));
	}
	return intens;
}
