#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <complex.h>
#include <string.h>





void HypergeompFqMat(int dima,double* (a),int dimb, double* (b),int coldimx,double complex (x)[2],int coldimy,double complex (y)[2],double alpha,int MAX,double(lam)[2]){
	

	int n=coldimx;
	double sumlam;	
	int lambda[n];

	//floor(MAX./(1:n)
	for(int i=0; i<n ; i++){
		 lambda[i] = MAX/(i+1);

	}
	
	//if ~empty (lam)
	if(lam[1] != 0 && lam[2] != 0){
		for(int i=0; i<n; i++){
			if(lambda[i]>lam[i]){
				lambda[i]=lam[i];
			}
			sumlam+=lam[i];
		}
		if((int)sumlam<MAX){
			MAX=(int)sumlam;
		}
	}


	int Lp=n;	
	while(lambda[Lp-1] == 0){
		Lp--;
	}	
	int lambda2[Lp];
	memcpy(&lambda2[0],&lambda[0],2*sizeof(int));


	
	int farr[MAX+1];
	int fval;
	for(int k=1; k<(MAX+2); k++){
		farr[k-1]=k;
	}

	for(int i=2; i<(n-1);i++){
		for(int j=(i+1); j<(MAX+1); j++){
			farr[j]= farr[j] + farr[j-1];
		}
	}

	fval=farr[MAX];	
	
	int D[fval];
	double complex Sx[fval][n];
	double complex Sy[fval][n];

	for(int i=0; i<fval; i++){
		D[i]=0;
	}
	for(int i=0; i<fval; i++){
		for(int j=0;j<n;j++){
			Sx[i][j]=0+0*I;
			if(i==0){
				Sx[i][j]=1+0*I;
			}
		}
	}
	
	/*    for X     */
	double complex xn[n][MAX+1];
	for(int k=0; k<n; k++){
		for(int j=0;j<(MAX+1); j++){
			xn[k][j]=1+0*I;	
			
		}
	}
	
	
	for(int i=1; i<(MAX+1);i++){
		for(int j=0; j<n;j++){
			xn[j][i]=xn[j][i-1]*x[j];
		}
		
	}	

	double complex prodx[coldimx];
	double complex prody[coldimx];
	prodx[0]=x[0];
	for(int i=1;i<coldimx;i++){
		prodx[i]=prodx[i-1]*x[i];
	}

	/*************************************/
	


	/*    for Y     */

	complex double yn[n][MAX+1];
	if(coldimy!=0){
		for(int i=0;i<fval;i++){
			memcpy(&Sy[i][0],&Sx[i][0],n*sizeof(double complex));
		}
		for(int k=0; k<coldimy; k++){
			for(int j=0;j<(MAX+1); j++){
				yn[k][j]=1+0*I;	
				
			}
		}
		
		
		for(int i=1; i<(MAX+1);i++){
			for(int j=0; j<n;j++){
				yn[j][i]=yn[j][i-1]*y[j];

			}
			
		}
		prody[0]=y[0];
		for(int i=1;i<coldimx;i++){
			prody[i]=prody[i-1]*y[i];
		}

	}
	/*************************************/

	/* DEFINITIONS (CHECK FOR INPUT ARRAY DIMENSIONS) */
	int rowdima=1;   // Change if a is  a matrix!
	double z[rowdima][Lp];
	//double (*z)[Lp] = malloc(sizeof(double[rowdima][Lp]));
	//memset(z,1,sizeof(double)*rowdima*Lp);
	for(int i=0;i<2;i++){ 
		memset(z[i], 1.0, Lp*sizeof(z[i][0]));
		for(int j=0;j<Lp;j++){
			z[i][j]=1;
		}
	}
	int* l = malloc(Lp * sizeof(int));
	int* ww = malloc(Lp * sizeof(int));
	int* kt = malloc(Lp * sizeof(int));
	int* d = malloc(Lp * sizeof(int));
	int* g = malloc(Lp * sizeof(int));
	int* mu = malloc(Lp * sizeof(int));
	int* mt = malloc(Lp * sizeof(int));
	memset(l, 0, Lp*sizeof(l[0]));
	memset(d, 0, Lp*sizeof(d[0]));
	memset(g, -1, Lp*sizeof(g[0]));
	for(int i=0;i<Lp;i++){
		ww[i]=1;
		kt[i]=-i-1;
	}
	int sl=1,h=0,heap=(lambda[0]+2);
	double complex  ss[rowdima][MAX+1];
	for(int i=0;i<rowdima;i++){
		ss[rowdima-1][0]=1+0*I;
		for(int j=1;j<MAX+1;j++){
			ss[rowdima-1][j]=0+0*I;
		}
	}

	
	int m;
	int pp;
	int w;
	double c;
	double complex cc1=0+0*I;
	double zn;
	double dn;
	
	double cc;
	int nhstrip;
	


	/*********** ALGORITHIM ***********/
	
	while(h>-1){
		if(l[h]<lambda[h] && (h==0||l[h]<l[h-1]) && MAX>=sl && z[h]){
			l[h]=l[h]+1;
			if(l[h]==1 && h>0 && h<(n-1)){
				D[ww[h]]=heap;
				ww[h] = heap;
				if(lambda[h]<(MAX-sl+l[h])){
					m = lambda[h];	
				}
				else{
				m=(MAX-sl+l[h]);
				}
				if(m<l[h-1]){
					heap+=m;
				}
				else{
					heap+=l[h-1];
				}
			}
			else{
				ww[h]=ww[h]+1;
			}


			w=ww[h];
			c=l[h]-1-(h)/alpha;
			zn=1;
			for(int i=0;i<dima;i++){
				zn=zn*(a[i]+c);
			}
			zn=zn*alpha;
			dn=1;
			for(int i=0;i<dimb;i++){
				dn=dn*(b[i]+c);
			}
			dn=dn*(kt[h]+h+2);


			int delta;
			if(coldimy!=0){ //if argy (or XY) exists
				zn=zn*alpha*l[h];
				dn=dn*(n+(alpha*c));
				if(h>0){
					for(int j=0;j<h;j++){
						delta=kt[j]-kt[h];
						zn=zn*delta;
						dn=dn*(delta-1);
					}
				}
				}

			kt[h]=kt[h]+alpha;

			
			if(h>0){
				for(int j=0;j<h;j++){
					delta=kt[j]-kt[h];
					zn=zn*delta;
					dn=dn*(delta+1);
				}
			}
			
			
			
			for(int i=0;i<rowdima;i++){
				z[i][h]=z[i][h]*(zn/dn);
			}
			sl=sl+1;



			if(h<n-1){
				if(h>0){
					d[h-1]=d[h-1]-1;
				}
				
				
				double val1=1,val2=1;
				for(int i=0;i<h+1;i++){
					val1=val1*(h+2-alpha+kt[i]);
					val2=val2*(h+1+kt[i]);
				}
				
				cc=val1/val2;
				d[h]=l[h];
				pp=l[0];
				int k=2;


				while(k<h+1 && l[k]>1){//unchecked
					pp=D[pp-1]+l[k]-2;
					k++;
				}

				Sx[w-1][h]=cc*prodx[h]*Sx[pp-1][h];
				if(coldimy!=0){ //if argy (or XY) exists
					Sy[w-1][h]=cc*prody[h]*Sy[pp-1][h];
				}



				int lg=0;
				for(int i=0;i<h+1;i++){
					if(d[i]>0){
						g[lg]=i;
						lg++;
					}
				}
				
				int slm=0;
				nhstrip=1;
				for(int i=0;i<lg;i++){
					nhstrip=nhstrip*(d[g[i]]);	
				}

				memcpy(&mu[0],&l[0],Lp*sizeof(int));
				memcpy(&mt[0],&kt[0],Lp*sizeof(int));
				

			
				double* blm = malloc( lg * sizeof(double));
				int* lmd = malloc( lg * sizeof(int));
				//INEFFICIENT	
				for(int i=0;i<lg;i++){
					blm[i]=1;
					lmd[i]=l[g[i]]-d[g[i]];
				}


				for(int i=0;i<nhstrip;i++){	

					int lgvar=lg-1;
					int gz=g[lg-1];

					while( mu[gz] == lmd[lgvar] ){
						mu[gz]=l[gz];
						mt[gz]=kt[gz];
						slm=slm-d[gz];
						lgvar=lgvar-1;
						gz=g[lgvar];
					}

					int t=kt[gz]-mt[gz];
					blm[lgvar]=blm[lgvar]*(t+1);
					dn=t+alpha;

					int q1,q2;//CHECK
					for (int r=0;r<gz-1;r++){
						q1=mt[r]-mt[gz];
						q2=kt[r]-mt[gz];
						blm[lgvar]=blm[lgvar]*(q1+alpha-1)*(1+q2);
						dn=dn*q1*(alpha+q2);
					}
					
					blm[lgvar]=blm[lgvar]/dn;
					mu[gz]=mu[gz]-1;
					mt[gz]=mt[gz]-alpha;
					slm=slm+1;

					if(lgvar<lg){//CHECK
						for(int j=lgvar;j<lg;j++){
							blm[j+1]=blm[j];
						}
					}
					int hvar=0;
					int nmu=mu[0];
					if(mu[h]==0){
						hvar=1;
					}
					for(int f=1;f<h-hvar;f++){
						nmu=D[nmu]+mu[f]-1;
					}

					for(int f=h+1;f<n;f++){
						Sx[w-1][f]=Sx[w-1][f]+(blm[lgvar]*Sx[nmu][f-1]*xn[f][slm]);
					}

					if(coldimy!=0){ //if argy (or XY) exists
						for(int f=h+1;f<n;f++){
							Sy[w-1][f]=Sy[w-1][f]+(blm[lgvar]*Sy[nmu][f-1]*yn[f][slm]);
						}
						
					}


				}


				for(int f=h+1;f<n;f++){
					Sx[w-1][f]=Sx[w-1][f]+Sx[w-1][f-1];
				}
				if(coldimy!=0){ //if argy (or XY) exists
					for(k=h+1;k<n;k++){
						Sy[w-1][k]=Sy[w-1][k]+Sy[w-1][k-1];
					}

					for(int i=0;i<rowdima;i++){
						ss[i][sl-1]=ss[i][sl-1]+(z[i][h]*Sx[w-1][n-1]*Sy[w-1][n-1]);
					}
				}
				else{
					for(int i=0;i<rowdima;i++){
						ss[i][sl-1]=ss[i][sl-1]+(z[i][h]*Sx[w-1][n-1]);
					}
				}
			}
			else{
				pp=l[0]+1-l[n-1];
				int k=1;
				while(k<n-1 && l[k]>l[n-1]){
					pp=D[pp]+l[k]-1-l[n-1];
					k++;
				}
				if(l[n-1]==1){
					k=1;
				}
				else{
					k=0;
				}

				int g1=1,g2=1;
				for(int i=0;i<n-1;i++){
					g1=g1*(1+kt[i]-kt[n-1]);
					g2=g2*(alpha+kt[i]-kt[n-1]);
				}

				cc1=k+(1-k)*cc1;
				double fn;
				if(coldimy!=0){ //if argy (or XY) exists
					fn=(g1*((1/alpha) +l[n-1]-1))/(g2*l[n-1]);
					cc1=cc1*prodx[n-1]*prody[n-1]*fn*fn;

					
					for(int i=0;i<rowdima;i++){
						ss[i][sl-1]=ss[i][sl-1]+(z[i][n-1]*Sx[pp-1][n-1]*Sy[pp-1][n-1]*cc1);
					}


				}
				else{
					fn=(g1*((1/alpha) +l[n-1]-1))/(g2*l[n-1]);
					cc1=cc1*prodx[n-1]*fn;

					for(int i=0;i<rowdima;i++){
						ss[i][sl-1]=ss[i][sl-1]+(z[i][n-1]*Sx[pp-1][n-1]*cc1);
					}
				}

				

			}
			if(h<Lp-1){
				for(int i=0;i<2;i++){
					z[i][h+1]=z[i][h];
				}
				ww[h+1]=w;
				h=h+1;
			}
		}
		else{
			sl=sl-l[h];
			l[h]=0;
			kt[h]=-(h+1);
			h=h-1;
		}
	//FREEING MEMORY
//	free(d);

	}


	double complex ans[rowdima],sum=0;

	for(int i=0;i<rowdima;i++){
		for(int j=0;j<MAX+1;j++){
			sum+=ss[i][j];
		}
		ans[i]=sum;
		sum=0;
	}

	
	for(int i=0;i<rowdima;i++){
		printf("%f + i%f\n",creal(ans[i]),cimag(ans[i]));
	}
}



int main(int argc,char *argv[]){

	double para[2] = {3,4};
	double parb[3] = {5,6,7};
	double complex argx[2] = {1+2*I,2+3*I};
	double complex argy[2] = {8+0*I,9+0*I};
	double hyperalpha = 9;
	int Partitions =6;
	double parlam[2] = {4,3};


	
	HypergeompFqMat(2,para,3,parb,2,argx,2,argy,hyperalpha,Partitions,parlam);

	return 0;
}
