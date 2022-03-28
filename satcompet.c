/*
   	Adriana Sanz 22/05/2017
   satcompet.c is a simple model of viral competition.
	
	Implements Runge Kutta integration step. 
	Fist order system of equations:
		h' = g - (d +p_x*x+p_y*y + p_sy*s)*h 
		x' = (p_x*h - (d + d_x))*x
		y' = (p_y*h - (d + d_y) - p_s*s)*y 
		s' = (p_s*y -(d + d_s) + p_sy*h)*s
	Variables
		h heathy hosts
		x infected hosts by virusx
		y infected hosts by virusy	
	Parameters
		g linear growth rate	
		p_i prob of get infected by i E {x,y}
		d_i prob of increased death being infected by i E {x,y}

compilation routine: gcc -std=c99 sat.c -o sat 

INCLUDES change of parameters
	
*/
////////////////////////////////////////
// HEADER
////////////////////////////////////////

# include <stdio.h>
# include <stdlib.h>
# include <time.h>

////////////////////////////////////////
// DEFINITIONS DECLARATIONS
////////////////////////////////////////

#define TIME 500000
//Uncomment the scenario you want to run 	
/*Commensalism*/
#define PATH4 "Data/e_gdx_h.dat"
#define PATH1 "Data/e_gdx_x.dat"
#define PATH2 "Data/e_gdx_y.dat"
#define PATH3 "Data/e_gdx_s.dat"  
/*Coexistence 
#define PATH4 "Data/c_gdx_h.dat"
#define PATH1 "Data/c_gdx_x.dat"
#define PATH2 "Data/c_gdx_y.dat"
#define PATH3 "Data/c_gdx_s.dat"*/
/*Bistability
#define PATH4 "Data/b_gdx_h.dat"
#define PATH1 "Data/b_gdx_x.dat"
#define PATH2 "Data/b_gdx_y.dat"
#define PATH3 "Data/b_gdx_s.dat"*/ 

#define U 500 // matrix size
#define V 500

//Declare functions

double dh(double *, double ,double ,double, double);
double dx(double *, double ,double ,double);
double dy(double *, double ,double ,double, double);
double ds(double *, double ,double ,double, double);
void rk(double *, double *, double ,double ,double, double );
void integration_rk(double [][4], double *, double *, double *, double);
void linspace(double *,double, double, int);

////////////////////////////////////////
// START
////////////////////////////////////////
int main(){

//parameters 
	double initial[]={1,0.2,0.2,1};	//vector for initial conditions
	double final[4];			//vector for final conditions
	double k[4][4];			//vector to calculate runge kutta
	double p[U];
	double q[V];
	linspace(p,0.0,0.5,V);
	linspace(q,0.0,0.5,U);
	double param[10];
	
	/*Common parameters*/
	param[0]=1;	// g
	param[1]=0.05;	// d
	param[2]=0.1;	// p_x
	param[3]=0.1;	// p_y
	param[4]=0.25;	// d_x
	param[5]=0.25;	// d_y
	param[6]=0.8;	// p_s
	param[7]=0.04;	// p_sy
	param[8]=0.3;	// d_s
	
	/*Commensalism*/
	param[9]=0.07666;//p_Y
	
	/*Bistability	
	param[9]=0.2;	// p_Y*/
	
	/*Coexistence
	param[9]=0.016;//p_Y*/	

	double tstep=0.01; // time step

////////////////////////////////////////	
// File IO
////////////////////////////////////////
	FILE *fp1, *fp2, *fp3, *fp4; 	// create a pointer "fp" to a type named FILE-->creates an empty file
	fp1=fopen(PATH1,"w"); 		// indicate where to save/name  the file and the option "w" 
	fp2=fopen(PATH2,"w");
	fp3=fopen(PATH3,"w");		
	fp4=fopen(PATH4,"w");
//SET TIME
	time_t start, stop;
	time(&start); 
	
	for (int i=0; i<V; i++){	// initial loops for change of parameters
		param[4]=q[i];
		for (int u=0;u<U;u++){
			initial[0]=1;
			initial[1]=1;// switch to .00001 for bistability
			initial[2]=1;	
			initial[3]=1;	
			param[0]=p[u];
		
			for (int t=0;t<TIME;t++){ 	//temporal loop
		
				//integrate runge kutta for each time
				integration_rk(k,param,initial,final,tstep);		


				//change the intial value
				initial[0]=final[0];
				initial[1]=final[1];
				initial[2]=final[2];
				initial[3]=final[3];
				if (t==TIME-1){				
					fprintf(fp1,"%.10f\t",final[1]);
					fprintf(fp2,"%.10f\t",final[2]);
					fprintf(fp3,"%.10f\t",final[3]);
					fprintf(fp4,"%.10f\t",final[0]);
				}		
			}
		}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
		fprintf(fp3,"\n");
		fprintf(fp4,"\n");
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	time(&stop);

	printf("\nElapsed time: %.0f seconds\n",difftime(stop, start));
	return 0;
}
//END

////////////////////////////////////////
//   FUNCTIONS    
////////////////////////////////////////
void integration_rk(double k[][4], double *variables, double *initial,double *final,double tstep){
//runge-kutta integrator
	
	double x_0,y_0,h_0,s_0;
	h_0=initial[0];
	x_0=initial[1];	
	y_0=initial[2];
	s_0=initial[3];

	rk(k[0],variables,h_0,x_0,y_0,s_0);

	rk(k[1],variables, h_0 + (tstep/2)*k[0][0], x_0 + (tstep/2)*k[0][1], y_0 + (tstep/2)*k[0][2] ,s_0 + (tstep/2)*k[0][3]);

	rk(k[2],variables, h_0 + (tstep/2)*k[1][0], x_0 + (tstep/2)*k[1][1], y_0 + (tstep/2)*k[1][2],s_0 + (tstep/2)*k[1][3]);
		
	rk(k[3],variables, h_0 + tstep*k[2][0], x_0 + tstep*k[2][1], y_0 + tstep*k[2][2],s_0 + tstep*k[2][3]);
	
	final[0]= h_0 + tstep*(k[0][0]+2*k[1][0]+2*k[2][0]+k[3][0])/6;
	final[1]= x_0 + tstep*(k[0][1]+2*k[1][1]+2*k[2][1]+k[3][1])/6;
	final[2]= y_0 + tstep*(k[0][2]+2*k[1][2]+2*k[2][2]+k[3][2])/6;
	final[3]= s_0 + tstep*(k[0][3]+2*k[1][3]+2*k[2][3]+k[3][3])/6;

}

void rk(double v[4], double p[8], double h,double x,double y,double s){
// Call the runge-kutta for each derivative in the system	
	v[0]=dh(p,h,x,y,s);
	v[1]=dx(p,h,x,s);
	v[2]=dy(p,h,x,y,s);
	v[3]=ds(p,h,x,y,s);

}
// Calculate the derivatives of the system
double dh(double *v,double h,double x,double y,double s){
	double g=v[0];
	double d=v[1];
	double p_x=v[2];
	double p_y=v[3];
	double p_sy=v[7];
	double p_Y=v[9];
	return g - (d + p_x*x + p_y*y + (p_sy+p_Y)*s)*h;
}
double dx(double *v,double h,double x,double y){
	double d=v[1];
	double p_x=v[2];
	double d_x=v[4];
	return (p_x*h-(d+d_x))*x;
}
double dy(double *v,double h,double x,double y,double s){
	double d=v[1];
	double p_y=v[3];
	double d_y=v[5];
	double p_s=v[6];
	double p_Y=v[9];
	return p_y*y*h - (d+d_y)*y - p_s*s*y + p_Y*s*h;
}
double ds(double *v,double h,double x,double y, double s){
	double d=v[1];
	double p_s=v[6];
	double p_sy=v[7];
	double d_s=v[8];
	return p_s*y*s -(d + d_s)*s + p_sy*h*s;
}


void linspace(double *v,double MIn, double MAx, int LENGTH){
//creates a vector of length LENGTH that spands from min to max 
	double s=(MAx-MIn)/(LENGTH-1);
	v[0]=MIn;
	for ( int i =1; i< LENGTH; i++ ){v[i]=s+v[i-1];}
}

