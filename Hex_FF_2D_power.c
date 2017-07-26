//
//  Hex_FF_2D_power.c
//  
//
//  Created by Emmanuel Malikides on 27/09/12.
//  Generates the far field output of a hexagonal-lattice optical phased array
//  with dark fibres. calculates the coupling factor in grace-like conditions.
//

//

// UNITS ARE IN MICRONS!
#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>

const double pi=3.14159265358979323846264338327950288419716939937510;

#define df 0
#define lf 1

// function that returns the value of a gaussian mode with source at xo yo on a 
// plane a distance z away, and at position x,y, phase shifted such that the peak 
// is at px,py.
complex double gmode(double x, double y, double z, double w0, double po, 
                     double lambda, double xo, double yo, double px, double py){
    double R=0.0;
    complex double mag;
    // calculate the intensity
    double E0=sqrt(2.0*po/pi/w0/w0);
    // calculate the rayleigh range
    double z0=pi*w0*w0/lambda;
    // calculate the waist size
    double wz=w0*sqrt(1.0+(z*z)/(z0*z0));
    // calculate the wavefront radius
    if (z!=0.0) R=z*(1.0+(z0*z0)/(z*z)); 
    // calculate radial vector squared over the grid
    double rhosq=(x-xo)*(x-xo)+(y-yo)*(y-yo);
    double ph =(px-xo)*(px-xo)+(py-yo)*(py-yo);
    // wavenumber
    double k=2.0*pi/lambda;
    // deal with zero case
    if (z==0.0) {
        mag=(complex double)E0*(w0/wz)*exp(-rhosq/wz/wz);
    } else {
        mag=((complex double)E0*(w0/wz)*exp(-rhosq/wz/wz))*
        cexp(-I*k*(rhosq-ph)/(2.0*R));
    }
    return mag;
}

// calculate the variance of a list
double var(double ns[], int length, double mean){
    double sum=0.0;
    int i=0;
    for (i=0;i<length;i++) {
        sum+=(ns[i]-mean)*(ns[i]-mean);
    }
    return sum/(double)length;
}

// calculates the mean of a list
double mean(double ns[], int length) {
    double sum=0.0;int i=0;
    for (i=0;i<length;i++) sum+=ns[i];
    return sum/(double)length;
}

// calculates a normally distributed random number.
double rand_normal(double mean, double stddev) {
    static double n2 = 0.0;
    static int n2_cached = 0;
    double x,y,r;
    
    if (!n2_cached) {
        do {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x+y*y;
        } while (r == 0.0 || r > 1.0);
        
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    } else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

// rotates a point counter clockwise about the origin by a given angle
void rotate(double source[3], double output[3], double angle) {
    output[0]=source[0]*cos(angle)+source[1]*sin(angle);
    output[1]=source[0]*sin(angle)-source[1]*cos(angle);
    output[2]=source[2];
    return;
}

// generates list of positions for hexagonal array, the centre of the pattern is
// always a dark fibre.
void hex_grid(double array[][3], int nRings, double pos[2],double rad, double error) {
    
    int tmp=0,tmp1=0,i,j;
    array[tmp][0]=pos[0];
    array[tmp][1]=pos[1];
    array[tmp][2]=df;
    tmp++;
    //Array stores the coorinates as well as the type of fibre (dark = 0,
    //light=1)
    int nTri=nRings*(nRings+1)/2;
    double triangle[nTri][3]; 
    double dx=rad/2,dy=rad*sqrt(3)/2;
    
    //generate a triangle
    for (i=1;i<=nRings;i++) {
        for (j=0;j<=i-1;j++) {
            triangle[tmp1][0]=(2*j-i)*dx+rand_normal(0,error);
            triangle[tmp1][1]=i*dy+rand_normal(0,error);
            triangle[tmp1][2]=((i+j)%3==0?df:lf);
            tmp1++;
        }
    }
    
    //add rotated versions of the triangle to the list to complete the hexagon
    for (i=0;i<6;i++) {
        for (j=0;j<nTri;j++) {
            rotate(triangle[j],array[tmp],i*pi/3);
            tmp++;
        }
    }
    return;
}


int main()

{
    
    //wavelength of light
    const double LL=1.064;
    //fibre aperture, 3.3 is Mode field Radius for thorlabs PM980-XP PM fibre
    const double WW=70;
    // calculate the rayleigh range
    double z0=pi*WW*WW/LL;
    //distance to far field
    //    const double ZZ=1000*z0;
    //220km
    const double ZZ=220000000000;
    printf("distance to far field: %e\n",ZZ);
    const double sep=150;
    //number of rings in hexagon
    const int Rings=27;
    //resolution of simulation grid
    const int  m=100;
    //receiving mode width in microns
    const double WR=1000;
    //length and width of simulation grid in microns
    const double l=2*WR;
    const double w=2*WR;
    
    //pointing position
    const double PX=0.0;
    const double PY=0.0;
    double pos[2]={PX,PY};
    const double dy=(double)l/(double)(m-1);
    const double dx=(double)w/(double)(m-1);
    int i=0;
    // number of sources
    const int NSources=6*(Rings*(Rings+1)/2)+1;
    double sources[NSources][3];
    hex_grid(sources,Rings,pos,sep,0);
    int NBright=0;
    for (i=0;i<NSources;i++) {
        NBright+=sources[i][2];
    }
    printf("The number of bright fibres was: %d\n",NBright);
    printf("The number of dark fibres was: %d\n",NSources-NBright);
    // intensity of each element.
    const double P0=1.0/NBright;
    double amplitudes[NSources];
    double badamplitudes[NSources];
    double phaseerrors[NSources];
    
    //stdev in position as percentage of spacing 
    const double vars[2]={0.05,0.15};
    const int origin=m/2;
    FILE *fp; 
    // output file 
    fp=fopen("out.dat","wt");
    
    // temporary variables to store scores, iterators and such 
    int j=0,k=0,ii=0,jj=0,kk=0;
    double x=0.0;
    double y,tmp,intensity, dphase;
    double z[m];
    double dphi=0;
    double meanpower=0;
    
    int NRuns=5;
    double powers[NRuns];
    // standard deviation of position error (in microns)
    double poserror=44;
    complex double tmode,rmode,cftemp;
    double recpower,tpower;
    srand((int)time(NULL));
    
    for (i=0;i<NRuns;i++) {
        
        // generate amplitude and phase errors.
        for (ii=0;ii<NSources; ii++) {
            // polarisation error, with mean 0 and standard deviation of 
            // 1 degree, or about 0.02 radians
            dphi=rand_normal(0,0.02);
            // intensity, a normally distributed random variable with mean 1 
            // and standard deviation of 5%
            amplitudes[ii]=cos(dphi)*abs(rand_normal(1.0,0.05));
            // amplitudes of signals with orthogonal polarisations
            badamplitudes[ii]=sin(dphi);
            phaseerrors[ii]=rand_normal(0,0.001);
        }
        
        //fill sources with relevant positions.
        hex_grid(sources,Rings,pos,sep,poserror);
        cftemp=recpower=tpower=0.0;
        
        //calculate the coupling factor over the area
        for (ii=0;ii<m;ii++) {
            for (jj=0;jj<m;jj++) {
                // generate coordinate which we are looking at.
                y=dy*jj-l/2;
                x=dx*ii-w/2;
                tmode=rmode=0.0;
                // if we are within the aperture, add to coupling factor.
                // add the contributions from each source at this position with 
                // some phase error.
                if(x*x+y*y<WR*WR) {
                    for (k=0;k<NSources;k++) {
                        dphase=phaseerrors[k];
                        tmode += 
                        cexp(I*dphase)*
                        (amplitudes[k]+badamplitudes[k])*
                        (sources[k][2]?
                         gmode(x, y, ZZ, WW, P0, LL, sources[k][0], sources[k][1], PX, PY)
                         :0.0);
                    }
//                    tmode=gmode(x, y, ZZ, WW, 1.0, LL, 0, 0, 0, 0);
                    rmode=gmode(x, y, 0.0, WR, 1.0, LL, 0, 0, 0, 0);
                    cftemp+=tmode*conj(rmode)*dx*dy;
                    recpower+=rmode*conj(rmode)*dx*dy;
                    tpower+=tmode*conj(tmode)*dx*dy;
                } else {                
                }
                //                fprintf(fp," %e",cabs(tmode)*cabs(tmode));
            }
            fprintf(fp,"\n");    
        }
        // normalise.
        printf("tpower %e\nrpower %e\n",tpower,recpower);
        powers[i]=cabs(cftemp)*cabs(cftemp)/recpower;
        printf("power coupled: %e\n",powers[i]);
    }
    printf("dx %e\ndy %e\n",dx,dy);
    meanpower=mean(powers,NRuns);
    printf("Power coupled: %e+- %e\n",meanpower,sqrt(var(powers,NRuns,meanpower)));
    
    fclose(fp);
    return 0;
}





