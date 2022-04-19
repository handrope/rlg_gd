#include<iostream>
#include<fstream>
#include<sstream>
#include"casuale.h"
#include<time.h>

using namespace std;

int Poisson(double) ;

int main(int argc, char *argv[]) {
    
    /* system initial parameters */
    int d=atoi(argv[1]);
    double phi=atof(argv[2]);   //scaled packing fraction of obstacles
    double A=atof(argv[3]);
    //long seed=atol(argv[4]);
    int nrun=10000 ;
    //int nrun=atoi(argv[4]);
    double R=1.+A/d; //interaction cutoff
    double deltaF_th=.99, dtadapt=.2 ;  // adaptive integration parameters
    
    /* run parameters */
    double tfin=1000.;
    //double tprint=0.01;
    //double tav = 1. ;
    double vstop=1e-8;
    double esat=vstop*vstop*1e3 ;
    //int freqprint=(int)(tprint/dt) ;
    //int freqav=(int)(tav/dt) ;
    long seed=-time(NULL);
    //long seed=-4378654;
    
    /* other constants */
    //int N=int(d*phi*pow(R,d));
    //cout<<"# of obstacles: "<<N<<endl;
    
    int NP=1000;
    double *t_av = (double*) calloc(NP,sizeof(double)) ;
    double dt0=0.001 ;
    
    double db=(1.+log(tfin/dt0))/((double)NP) ;
    double N1=1.+1./db+log(db*tfin/dt0)/log(1.+db) ;
    while(fabs(1.-N1/NP)>1e-4) {
        db *=NP/(2.*NP-N1) ;
        N1=1.+1./db+log(db*tfin/dt0)/log(1.+db) ;
    }
    double b=1.+db ;
        
    int np=0;
    t_av[np]=0. ; np++ ;
    t_av[np]=dt0 ;
    while(t_av[np]*b<t_av[np]+dt0&&np<NP-1) {
        np++ ;
        t_av[np] = t_av[np-1]+dt0 ;
    }
    double t1=t_av[np] ;
    while(np<NP-1) {
        np++ ;
        t1 *= b ;
        t_av[np] = floor(t1/dt0)*dt0 ;
    }

    double *msd_av = (double*) calloc(NP,sizeof(double)) ;
    double *v_av = (double*) calloc(NP,sizeof(double)) ;
    double *e_av = (double*) calloc(NP,sizeof(double)) ;
    double *c_av = (double*) calloc(NP,sizeof(double)) ;
    double *c1_av = (double*) calloc(NP,sizeof(double)) ;
    
    double N0= d*phi*pow(R,d);
    for (int ru=0;ru<nrun;ru++) {
    
        int N=Poisson(N0);
        
        /* initialize the obstacles */
        double **xi = (double**) calloc(N,sizeof(double*));
        for(int i=0;i<N;i++) xi[i] = (double*) calloc(d,sizeof(double)) ;
        for (int i=0;i<N;i++) {
            double rtemp=0;
            for (int mu=0;mu<d;mu++) {
                xi[i][mu]=gasdev(&seed);
                rtemp+=xi[i][mu]*xi[i][mu];
            }
            double rfin=R*pow(ran1(&seed),1./d)/sqrt(rtemp);
            for (int mu=0;mu<d;mu++) {
                xi[i][mu]*=rfin;
            }
        }
        
        double *X = (double*) calloc(d,sizeof(double)) ;
        double *X1 = (double*) calloc(d,sizeof(double)) ;
        
        double *h = (double*) calloc(N,sizeof(double)) ;
        double **r = (double**) calloc(N,sizeof(double*)) ;
        for(int i=0;i<N;i++) r[i] = (double*) calloc(d,sizeof(double)) ;
        double *Fc = (double*) calloc(d,sizeof(double)) ;
        double *Fa = (double*) calloc(d,sizeof(double)) ;
        double msd, v, e, p, c, c1 ;
        double msd0, v0, e0, p0, c0, c10, t0 ;
        
        /* Initialize the system */
        
        /* compute the hi and ri */
        for (int i=0;i<N;i++) {
            h[i]=0;
            for (int mu=0;mu<d;mu++) {
                r[i][mu]=X[mu]-xi[i][mu];
                h[i]+=r[i][mu]*r[i][mu];
            }
            h[i]=sqrt(h[i]);
            for (int mu=0;mu<d;mu++) {
                r[i][mu]/=h[i];
            }
            h[i]=d*(h[i]-1);
        }
        /* compute the conservative force */
    
        for (int mu=0;mu<d;mu++) Fc[mu]=0;
        
        for (int i=0;i<N;i++) {
            if (h[i]<0) {
                for (int mu=0;mu<d;mu++)     Fc[mu]-=h[i]*r[i][mu]/d;
            }
        }
        
        msd=0;
        v=0;
        double FF=0. ;
        for (int mu=0;mu<d;mu++) {
            FF += Fc[mu]*Fc[mu] ;
            msd+=X[mu]*X[mu];
        }
        v=sqrt(FF);
        
        /* Time loop */
        int np=0, npp=0 ;
        double t=0., dt ;
        do {
            
            /* output */
            e=0;
            p=0;
            c=0;
	    c1=0;
            for (int i=0;i<N;i++) {
                if (h[i]<-0.0) {
                    e+=h[i]*h[i]/2.;
                    p-=h[i];
                    c++;
                }
		if (h[i]<-vstop) c1++ ;
            }
	   
            
            //if(t>=tprint*npp) {
	    //cout<<"t"<<ru<<"t "<<t<<" "<<c<<" "<<p<<" "<<" "<<e<<" "<<v<<" "<<msd<<endl;
	    //npp++ ;
	    //}
            

            if(np<NP) {
                if(t>=t_av[np]) {
                    if(t!=t_av[np]) {
                        double t0=t-dt ;
                        msd=(msd0*(t-t_av[np])+msd*(t_av[np]-t0))/(t-t0);
                        c=(c0*(t-t_av[np])+c*(t_av[np]-t0))/(t-t0);
                        v=(v0*(t-t_av[np])+v*(t_av[np]-t0))/(t-t0);
                        e=(e0*(t-t_av[np])+e*(t_av[np]-t0))/(t-t0);
			c1=(c10*(t-t_av[np])+c1*(t_av[np]-t0))/(t-t0);
                    }
                    msd_av[np] = (ru*msd_av[np]+msd)/(ru+1) ;
                    c_av[np] = (ru*c_av[np]+c)/(ru+1) ;
                    v_av[np] = (ru*v_av[np]+v)/(ru+1) ;
                    e_av[np] = (ru*e_av[np]+e)/(ru+1) ;
		    c1_av[np] = (ru*c_av[np]+c1)/(ru+1) ;
                    np++ ;
                }
            }
            msd0=msd ; c0 = c ; v0 = v ; e0 = e ; c10 = c1 ;
            
            // Adaptive algorithm
            dt=dt0;  // bare integration time step
            double deltaF, Fpost;
            do {
                /* compute the attempted move */
                for (int mu=0;mu<d;mu++) X1[mu] = X[mu] + Fc[mu]*dt/d;
                
                /* compute the new gaps */
                for (int i=0;i<N;i++) {
                    h[i]=0;
                    for (int mu=0;mu<d;mu++) {
                        r[i][mu]=X1[mu]-xi[i][mu];
                        h[i]+=r[i][mu]*r[i][mu];
                    }
                    h[i]=sqrt(h[i]);
                    for (int mu=0;mu<d;mu++) {
                        r[i][mu]/=h[i];
                    }
                    h[i]=d*(h[i]-1);
                }
                /* compute the new conservative force */
                for (int mu=0;mu<d;mu++) Fa[mu]=0;
                
                for (int i=0;i<N;i++) {
                    if (h[i]<0) {
                        for (int mu=0;mu<d;mu++) Fa[mu]-=h[i]*r[i][mu];
                    }
                }
                /* compute the force square modulus */
                deltaF=0. ; Fpost=0. ;
                for(int mu=0;mu<d;mu++) {
                    deltaF += Fc[mu]*Fa[mu] ;
                    Fpost += Fa[mu]*Fa[mu] ;
                }
                deltaF /= sqrt(FF*Fpost);
                dt *= dtadapt ;
                //cout<<scientific<<"# "<<deltaF<<" "<<dt<<endl;
            }while(deltaF<deltaF_th);
            
            /* update the configuration */
	    msd=0 ;
            for(int mu=0;mu<d;mu++) {
                X[mu] = X1[mu] ;
                Fc[mu] = Fa[mu] ;
                msd += X[mu]*X[mu] ;
            }
            dt /= dtadapt ;
            FF=Fpost ;
            v=sqrt(FF) ;
            //msd=sqrt(msd) ;
            t+=dt ;
            
        }while(v>vstop&&t<tfin) ;
        
        while(np<NP) {
            msd_av[np] = (ru*msd_av[np]+msd)/(ru+1) ;
            c_av[np] = (ru*c_av[np]+c)/(ru+1) ;
            v_av[np] = (ru*v_av[np]+v)/(ru+1) ;
            e_av[np] = (ru*e_av[np]+e)/(ru+1) ;
	    c1_av[np] = (ru*c1_av[np]+c1)/(ru+1) ;
            np++ ;
        }
        cout<<"t"<<ru<<"t "<<t<<" "<<c<<" "<<p<<" "<<" "<<e<<" "<<v<<" "<<msd<<" "<<N<<" "<<c1<<endl;
        //cout<<"\n\n";
        for(int i=0;i<N;i++) { free(xi[i]) ; free(r[i]) ; }
        free(xi) ; free(X) ; free(h) ; free(r) ; free(Fc) ;
        free(X1) ; free(Fa) ;
        
        ofstream fV ;
        fV.open("ave.dat");
        for(int np=0;np<NP;np++) fV<<t_av[np]<<" "<<c_av[np]<<" "<<v_av[np]<<" "<<e_av[np]<<" "<<msd_av[np]<<" "<<c1_av[np]<<endl ;
        fV.close() ;
    }
    free(t_av) ;
    free(msd_av) ;
    free(c_av) ;
    free(v_av) ;
    free(e_av) ;
    free(c1_av) ;
    
    return 0;
    
}


int Poisson(double lambda) {
    
    long seed=-time(NULL) ;
    double STEP=500. ;
    
    double p=1. ;
    int k=0 ;
    do {
        k++ ;
        p *= ran1(&seed);
        while(p<1&&lambda>0) {
            if(lambda>STEP) {
                p *= exp(STEP) ;
                lambda -= STEP ;
            }
            else {
                p *= exp(lambda) ;
                lambda=0. ;
            }
        }
    }while(p>1.);

    return (k-1) ;
    
}
