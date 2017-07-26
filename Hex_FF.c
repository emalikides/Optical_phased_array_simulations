//
//  FarField.c
//  
//
//  Created by Emmanuel Malikides on 18/04/12.  
//  Generates the far field output of a hexagonal-lattice optical phased array
//  with dark fibres.
//

// UNITS ARE IN MICRONS!
#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>


#define nvar 10
//randomisation resolution
#define res 10
//spacing resolution-number of increments between maximum and minimum element spacing
#define nsp 100
//spacing between the fibres
#define SPMAX 2000
#define SPMIN 200

//length and width of simulation grid in microns
#define l 20
#define w 20000
//resolution of simulation grid
#define n 10
#define m 200
//fibre aperture
const double WW=6.2;
//number of rings in hexagon
#define Rings 27

const double pi=3.14159265358979323846264338327950288419716939937510;

#define df 0
#define lf 1

// function that returns the value of a gaussian mode with source at xo yo on a plane a distance z away, and at position x,y, phase shifted such that the peak is at px,py.
complex double gmode(double x, double y, double z, double w0, double po, double lambda, double xo, double yo, double px, double py){
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
void hex_grid(double array[][3], int nRings, double pos[2],double rad) {
    
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
            triangle[tmp1][0]=(2*j-i)*dx;
            triangle[tmp1][1]=i*dy;
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

// Helper functions to find the radius of the peak.
int investigatex(int x, int y, double power, double intensity[n][m],int dir) {
    int r=x;
    do {
        r+=dir;
    } while ((intensity[r][y]>power)&&(r<n));
    return r;
}

int investigatey(int x, int y, double power, double intensity[n][m],int dir) {
    int r=y;
    do {
        r+=dir;
    } while ((intensity[x][r]>power)&&(r<m));
    return r;
}

// find the radius of the peak in grid units by finding the distance before the intensity
// drops to "power" in the x and y directions.
int findrad(int x, int y, double power, double intensity[n][m]){
//    int xp=investigatex(x,y,power,intensity,1);
//    int xn=investigatex(x,y,power,intensity,-1);
    int yp=investigatey(x,y,power,intensity,1);
    int yn=investigatey(x,y,power,intensity,-1);
    //printf("xp=%d,xn=%d,yp=%d, yn=%d, av=%d\n",xp,xn,yp,yn,(xp-xn+yp-yn)/4);
    return (yp-yn)/4;
}

//double findpeak(int x, int y, double power, double intensity[n][m],int width){    
//    int r=width;
//    printf("r: %d\n",r);
//    do {
//        r+=1;
//    } while ((intensity[x][y+r]<power)&&(r<n));
//    return r;
//}

//Find the distance between x and the nearest peak with power greater than that given by searching through list intensity.
int findpeak(int x,int y, double power, double intensity[n][m]){
    int r=y;   
    int left=0;
    do {
        r+=1;
        if ((intensity[x][y+r]<power)&&(!left)) left=1;
    } while (((!left)||(intensity[x][r]<power)||(intensity[x][r+1]-intensity[x][r]>0.0))&&(r+1<n));
    int result=r-y;
    if (result>n) printf("bad news,distance too large...\n");
    //    printf("power at peak=%e\n",intensity[r]);
    return result;
}
int main()

{
    //distance to far field
    const double ZZ=2000000;
    const int NSources=6*(Rings*(Rings+1)/2)+1;
    const int NBright=6*Rings;
    double amplitudes[NSources];
    double badamplitudes[NSources];
    
    const double LL=1.064;
    // power of each element.
    const double P0=1.0/NBright;
    //pointing position
    const double PX=0.0;
    const double PY=0.0;
    double pos[2]={PX,PY};
    const double dx=(double)w/(double)(n-1);
    const double dy=(double)l/(double)(m-1);

    const double dsp=(SPMAX-SPMIN)/(nsp-1);

    //stdev in position as percentage of spacing 
    const double vars[2]={0.05,0.15};
    const int oo[2]={w/2/dx,l/2/dy};
    FILE *fp; 
    // output file 
    fp=fopen("out.dat","wt");

    // temporary variables to store scores, iterators and such 
    int i=0,j=0,k=0,ii=0,jj=0,kk=0;
    double x,y,tmp,tmps,power,sources[NSources][3];
    double z[n][m];
    double peaktmp[res];
    double dphi=0;
    double dphase=0;
    double meanwidth=0;
    int NRuns=4;
    double rad[NRuns];
    double dpeak[NRuns];
    double powers[NRuns];
    double meandpeak=0;
    complex double etemp;
    srand((int)time(NULL));
    
    printf("dx=%e, dy=%e, expected width of peak=%e\n",dx,dy,ZZ*LL/pi/WW>WW?ZZ*LL/pi/WW:WW);
    
    // generate amplitudes with random numbers. 
    for (i=0;i<NSources; i++) {
        // polarisation error, with mean 0 and standard deviation of 1 degree, or about 0.02 radians
        dphi=rand_normal(0,0.02);
        // power, a normally distributed random variable with mean 1 and standard deviation of 5%
        amplitudes[i]=cos(dphi)*abs(rand_normal(1.0,0.05));
        badamplitudes[i]=sin(dphi);
    }
    
    for (i=0;i<NRuns;i++) {
        tmps=300;
        //fill sources with relevant positions.
        hex_grid(sources,Rings, pos,tmps);
        //calculate the peak power for this hexagon
        for (k=0;k<NSources;k++) {
            etemp += (amplitudes[k]+badamplitudes[k])*(sources[k][2]?gmode(PX, PY, ZZ, WW, P0, LL, sources[k][0], sources[k][1], PX, PY):0.0);
        }
        power=cabs(etemp)*cabs(etemp);
        //calculate the fringe pattern over the peak.
        for (ii=0;ii<n;ii++) {
            for (jj=0;jj<m;jj++) {
                // generate coordinate at which we are looking at.
                x=dx*ii-w/2+PX;
                y=dy*jj-l/2+PY;
                etemp=0.0;
                // add the contributions from each source at this position
                for (k=0;k<NSources;k++) {
                    //                        dphase=rand_normal(0,0.001);
                    //                        etemp += cexp(I*dphase)*amplitudes[k]*(sources[k][2]?gmode(x, y, ZZ, WW, P0, LL, sources[k][0], sources[k][1], PX, PY):0.0);
                    etemp += amplitudes[k]*(sources[k][2]?gmode(x, y, ZZ, WW, P0, LL, sources[k][0], sources[k][1], PX, PY):0.0);
                }
                z[ii][jj]=cabs(etemp)*cabs(etemp);
                if (i==0) {fprintf(fp,"%e %e %e\n",x/ZZ,y/ZZ,z[ii][jj]);}
            }
        }
        //store the current peak power in an array
        powers[i]=power;
        printf("power: %e\n",power);
        //store the current peak width in radians in an array
        rad[i]=findrad(oo[0],oo[1], power/2, z);
        //store the distance to the nearest peak of 1/4 the peak height in radians.
        dpeak[i]=findpeak(oo[0],oo[1],power/2,z)*dx/ZZ;
        printf("steering range: %e\n",dpeak[i]);
        //convert rad to radians
        rad[i]=rad[i]*dx/ZZ;
        printf("width of peak in radians: %e\n",rad[i]);
    }
    power=mean(powers,NRuns);
    meanwidth=mean(rad,NRuns);
    meandpeak=mean(dpeak,NRuns);
    printf("mean peak width (FWHM): %e +- %e\n", meanwidth, var(rad,NRuns,meanwidth));
    printf("mean peak Intensity, or efficiency: %e +- %e \n", power,var(powers,NRuns,power));
    printf("mean distance to nearest peak: %e +- %e \n", meandpeak, var(dpeak, NRuns,meandpeak));
    printf("width of grid in radians: %e by %e\n",l/ZZ,w/ZZ);
    fclose(fp);
    return 0;
}





