/*
### FctqLDA.cpp  (2011-04-23)
###    
###
### Copyright 2011 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA
*/

#include <vector>
#include "FctqLDA.h"

#define TRUE 1
#define FALSE 0

inline double dmax(double a,double b) { return( a > b ? a : b ); }
inline int imin(int a,int b) { return( a < b ? a : b ); }
inline int imax(int a,int b) { return( a > b ? a : b ); }
inline int imin3(int a,int b,int c) { return( a < b ? (a < c ? a : c) : (b < c ? b : c ) ); }
inline int iabs(int a) { return( a > 0 ? a : -a); }

double FnDist(double *B,int p,int q,double *Sigma,vector<int>& ind,double k0,double penF)
{
	double dist,halfsqdist=0.,viol=0.;

    { for (int i=0;i<p;i++) ind[i] = i*p-i*(i-1)/2; }	
	for (int i=0;i<p;i++) for (int j=0;j<=i;j++) {
		dist = Sigma[ind[j]+i-j] - B[i]*B[j];
		for (int a=1;a<=imin3(i,j,q-1);a++) dist -= B[ind[a]+i-a]*B[ind[a]+j-a];
		if (j<i) halfsqdist += dist*dist;
		else viol += dmax(k0-dist,0.);
	}
    return(2*halfsqdist + penF*viol*viol);
} 

int fgrad(double *B,vector<double>& Bij,vector<int>& iviol,vector<double>& fg,int p,int q,double *Sigma,vector<int>& ind,double k0,double penF)
{
	int indij;
	double tmp,dist,viol=0.;

	{ for (int i=0;i<p;i++)  ind[i] = i*p-i*(i-1)/2; }	
	for (int j=0;j<p;j++) for (int i=j;i<p;i++)  {
		Bij[indij=ind[j]+i-j] = B[j]*B[i];
		for (int b=1;b<=imin3(i,j,q-1);b++) Bij[indij] += B[ind[b]+j-b]*B[ind[b]+i-b];
	}
	for (int i=0;i<p;i++) {
		dist = Sigma[ind[i]] - Bij[ind[i]];
		if (dist < k0)  {
			iviol[i] = TRUE;
			viol += k0-dist;
		}
		else iviol[i] = FALSE;
	}
    for (int a=0;a<q;a++) for (int j=a;j<p;j++)  {
		tmp = 0.;
		for (int i=a;i<p;i++) if(i!=j) {
			indij = ind[imin(i,j)] + iabs(j-i);
			tmp += B[ind[a]+i-a] * (Sigma[indij] - Bij[indij]);
		}
		fg[ind[a]+j-a] = -4*tmp;
		if (iviol[j]) fg[ind[a]+j-a] += 4*penF*viol*B[ind[a]+j-a];
    }
    return(0);
}

int fhess(double *B,vector<double>& Baij,vector<double>& Basq,vector<double>& Bij,vector<double>& Bjl,vector<double>& Bab,vector<double>& Babj,
		vector<int>& iviol,vector<double>& fh,int p,int q,double *Sigma,
		vector<int>& ind,vector<int>& ind1,vector<int>& ind2,vector<int>& ind3,vector<int>& ind4,vector<int>& ind5,double k0,double penF)
{
	int npar,combij,indaj,indbj,indbl,ind4ij,indab,indabj;
	double dist,viol=0.;

    { for (int i=0;i<p;i++) ind[i] = i*p-i*(i-1)/2; }	
	npar = q*(q+1)/2 + (p-q)*q; 
	{ for (int i=0;i<npar;i++) ind1[i] = i*npar-i*(i-1)/2; }
	combij = p*(p+1)/2 ; 
	{ for (int i=0;i<p;i++) ind2[i] = i*combij-i*(i-1)/2; }
    { for (int i=0;i<q;i++) ind3[i] = i*q-i*(i-1)/2; }	
	{ for (int i=0;i<=p;i++) ind4[i] = i*(p-1)-i*(i-1)/2; }
	{ int ind0=0; for (int i=0;i<q-1;i++) ind0 = ind5[i] = ind0+ind4[q-i]; }

	{ for (int a=0;a<q;a++) {
		Basq[a] = 0.;
		for (int j=a;j<p;j++)  {
			for (int i=j;i<p;i++) Baij[ind2[a]+ind[j-a]+i-j] = B[ind[a]+j-a]*B[ind[a]+i-a];
			Basq[a] += Baij[ind2[a]+ind[j-a]];
		}
	} }
	for (int j=0;j<p;j++) for (int i=j+1;i<p;i++)  {
		Bij[ind4ij=ind4[j]+i-j-1] = Baij[ind[j]+i-j];		
		for (int a=1;a<=imin3(i,j,q-1);a++) Bij[ind4ij] += Baij[ind2[a]+ind[j-a]+i-j];
// Note : Replacing Bij by Bi (summing later) may be more efficient
	}
	{ for (int a=0;a<q;a++) for (int b=a+1;b<q;b++)  {
 		Bab[indab=ind3[a]+b-a-1] = Babj[indabj=ind5[a]+ind4[b-a-1]] = B[ind[a]+b-a]*B[ind[b]];	
		for (int i=b+1;i<p;i++) Bab[indab] += (Babj[indabj+i-b]=B[ind[a]+i-a]*B[ind[b]+i-b]);
	} }

	for (int a=0;a<q;a++) for (int j=a;j<p;j++)  {
		indaj = ind[a]+j-a;
		fh[ind1[indaj]] = 4*(Basq[a]-Baij[ind2[a]+ind[j-a]]); 
		for (int l=j+1;l<p;l++)  fh[ind1[indaj]+ind[a]+l-a-indaj] = 
				-4 * ( Sigma[ind[j]+l-j] - Baij[ind2[a]+ind[j-a]+l-j] - Bij[ind4[j]+l-j-1] );
		for (int b=a+1;b<q;b++)  {
			indbj = ind[b]+j-b;
			if (j>=b) fh[ind1[indaj]+indbj-indaj] = 
				4 * ( Bab[ind3[a]+b-a-1] - Babj[ind5[a]+ind4[b-a-1]+j-b] );
			for (int l=b;l<p;l++)  {
				indbl = ind[b]+l-b;
				if (indbl>indaj && l!=j) 
					if (j>=b) fh[ind1[indaj]+indbl-indaj] = 4 * B[ind[a]+l-a]*B[ind[b]+j-b];
					else fh[ind1[indaj]+indbl-indaj] = 0.;
			}
		}
	}

	for (int i=0;i<p;i++)  {
		dist = Sigma[ind[i]] - Baij[ind[i]];
		for (int a=1;a<imin(q,i+1);a++) dist -= Baij[ind2[a]+ind[i-a]];
		if (dist < k0)  {
			iviol[i] = TRUE;
			viol += k0-dist;
		}
	}
	{ for (int j=0;j<p;j++) if (iviol[j]) {
		for (int a=0;a<imin(q,j+1);a++) {
			indaj = ind[a]+j-a;
			fh[ind1[indaj]] += 4*penF*(viol+2*Baij[ind2[a]+ind[j-a]]);
			for (int l=j+1;l<p;l++) if (iviol[l]) 
				fh[ind1[indaj]+ind[a]+l-a-indaj] += 8*penF*Baij[ind2[a]+ind[j-a]+l-j];
			for (int b=a+1;b<q;b++) {
				indbj = ind[b]+j-b;
				if (b<=j) fh[ind1[indaj]+indbj-indaj] += 8*penF*Babj[ind5[a]+ind4[b-a-1]+j-b];
				for (int l=b;l<p;l++)  if (iviol[l]) {
					indbl = ind[b]+l-b;
					if (indbl>indaj && l!=j) 
						fh[ind1[indaj]+indbl-indaj] += 8*penF*B[ind[a]+j-a]*B[ind[b]+l-b];
				}
			}
		}
	} }
	return(0);
}	

double FnDist1(double *B,int p,int q,double *SigmSr,int Srank,vector<int>& ind,double k0,double penF)
{
	double Sigmaij,dist,halfsqdist=0.,viol=0.;
	int indi,indj;

	{ for (int i=0;i<p;i++) ind[i] = i*p-i*(i-1)/2; }	
	for (int i=0;i<p;i++) for (int j=0;j<=i;j++) {
		Sigmaij = SigmSr[indi=i*Srank] * SigmSr[indj=j*Srank];
		for (int k=1;k<Srank;k++)  Sigmaij += SigmSr[indi+k] * SigmSr[indj+k]; 		
		dist = Sigmaij - B[i]*B[j];
		for (int a=1;a<=imin3(i,j,q-1);a++) dist -= B[ind[a]+i-a]*B[ind[a]+j-a];
		if (j<i) halfsqdist += dist*dist;
		else viol += dmax(k0-dist,0.);
	}
    return(2*halfsqdist + penF*viol*viol);
} 

int fgrad1(double *B,vector<double>& Bij,vector<int>& iviol,vector<double>& fg,int p,int q,double *SigmSr,int Srank,vector<int>& ind,double k0,double penF)
{
	int indi,indj,indij;
	double Sigmaij,tmp,dist,viol=0.;

	{ for (int i=0;i<p;i++)  ind[i] = i*p-i*(i-1)/2; }	
	for (int j=0;j<p;j++) for (int i=j;i<p;i++)  {
		Bij[indij=ind[j]+i-j] = B[j]*B[i];
		for (int b=1;b<=imin3(i,j,q-1);b++) Bij[indij] += B[ind[b]+j-b]*B[ind[b]+i-b];
	}
	for (int i=0;i<p;i++) {
		Sigmaij = (tmp=SigmSr[indi=i*Srank])*tmp;
		for (int k=1;k<Srank;k++)  Sigmaij += (tmp=SigmSr[indi+k])*tmp; 		
		dist = Sigmaij - Bij[ind[i]];
		if (dist < k0)  {
			iviol[i] = TRUE;
			viol += k0-dist;
		}
		else iviol[i] = FALSE;
	}
    for (int a=0;a<q;a++) for (int j=a;j<p;j++)  {
		tmp = 0.;
		for (int i=a;i<p;i++) if(i!=j) {
			indi = imin(i,j)*Srank;
			indj = imax(i,j)*Srank;
			Sigmaij = SigmSr[indi]*SigmSr[indj];
			for (int k=1;k<Srank;k++)  Sigmaij += SigmSr[indi+k]*SigmSr[indj+k];; 		
			indij = ind[imin(i,j)] + iabs(j-i);
			tmp += B[ind[a]+i-a] * (Sigmaij - Bij[indij]);
		}
		fg[ind[a]+j-a] = -4*tmp;
		if (iviol[j]) fg[ind[a]+j-a] += 4*penF*viol*B[ind[a]+j-a];
    }
    return(0);
}

int fhess1(double *B,vector<double>& Baij,vector<double>& Basq,vector<double>& Bij,vector<double>& Bjl,vector<double>& Bab,vector<double>& Babj,
		vector<int>& iviol,vector<double>& fh,int p,int q,double *SigmSr,int Srank,
		vector<int>& ind,vector<int>& ind1,vector<int>& ind2,vector<int>& ind3,vector<int>& ind4,vector<int>& ind5,double k0,double penF)
{
    int npar,combij,indi,indj,indl,indaj,indbj,indbl,ind4ij,indab,indabj;
    double Sigmaij,tmp,dist,viol=0.;

    { for (int i=0;i<p;i++) ind[i] = i*p-i*(i-1)/2; }	
	npar = q*(q+1)/2 + (p-q)*q; 
	{ for (int i=0;i<npar;i++) ind1[i] = i*npar-i*(i-1)/2; }
	combij = p*(p+1)/2 ; 
	{ for (int i=0;i<p;i++) ind2[i] = i*combij-i*(i-1)/2; }
    { for (int i=0;i<q;i++) ind3[i] = i*q-i*(i-1)/2; }	
	{ for (int i=0;i<=p;i++) ind4[i] = i*(p-1)-i*(i-1)/2; }
	{ int ind0=0; for (int i=0;i<q-1;i++) ind0 = ind5[i] = ind0+ind4[q-i]; }

	{ for (int a=0;a<q;a++) {
		Basq[a] = 0.;
		for (int j=a;j<p;j++)  {
			for (int i=j;i<p;i++) Baij[ind2[a]+ind[j-a]+i-j] = B[ind[a]+j-a]*B[ind[a]+i-a];
			Basq[a] += Baij[ind2[a]+ind[j-a]];
		}
	} }
	for (int j=0;j<p;j++) for (int i=j+1;i<p;i++)  {
		Bij[ind4ij=ind4[j]+i-j-1] = Baij[ind[j]+i-j];		
		for (int a=1;a<=imin3(i,j,q-1);a++) Bij[ind4ij] += Baij[ind2[a]+ind[j-a]+i-j];
// Note : Replacing Bij by Bi (summing later) may be more efficient
	}
	{ for (int a=0;a<q;a++) for (int b=a+1;b<q;b++)  {
 		Bab[indab=ind3[a]+b-a-1] = Babj[indabj=ind5[a]+ind4[b-a-1]] = B[ind[a]+b-a]*B[ind[b]];	
		for (int i=b+1;i<p;i++) Bab[indab] += (Babj[indabj+i-b]=B[ind[a]+i-a]*B[ind[b]+i-b]);
	} }

	for (int a=0;a<q;a++) for (int j=a;j<p;j++)  {
		indaj = ind[a]+j-a;
		fh[ind1[indaj]] = 4*(Basq[a]-Baij[ind2[a]+ind[j-a]]); 
		{ for (int l=j+1;l<p;l++)  { 
			indj = j*Srank;
			indl = l*Srank;
			Sigmaij = SigmSr[indj]*SigmSr[indl];
			for (int k=1;k<Srank;k++)  Sigmaij += SigmSr[indj+k]*SigmSr[indl+k];; 		
			fh[ind1[indaj]+ind[a]+l-a-indaj] = 
				-4 * ( Sigmaij - Baij[ind2[a]+ind[j-a]+l-j] - Bij[ind4[j]+l-j-1] ); 
		} }
		for (int b=a+1;b<q;b++)  {
			indbj = ind[b]+j-b;
			if (j>=b) fh[ind1[indaj]+indbj-indaj] = 
				4 * ( Bab[ind3[a]+b-a-1] - Babj[ind5[a]+ind4[b-a-1]+j-b] );
			for (int l=b;l<p;l++)  {
				indbl = ind[b]+l-b;
				if (indbl>indaj && l!=j) 
					if (j>=b) fh[ind1[indaj]+indbl-indaj] = 4 * B[ind[a]+l-a]*B[ind[b]+j-b];
					else fh[ind1[indaj]+indbl-indaj] = 0.;
			}
		}
	}

	for (int i=0;i<p;i++)  {
		Sigmaij = (tmp=SigmSr[indi=i*Srank])*tmp;
		for (int k=1;k<Srank;k++)  Sigmaij += (tmp=SigmSr[indi+k])*tmp; 		
		dist = Sigmaij - Baij[ind[i]];
		for (int a=1;a<imin(q,i+1);a++) dist -= Baij[ind2[a]+ind[i-a]];
		if (dist < k0)  {
			iviol[i] = TRUE;
			viol += k0-dist;
		}
	}
	{ for (int j=0;j<p;j++) if (iviol[j]) {
		for (int a=0;a<imin(q,j+1);a++) {
			indaj = ind[a]+j-a;
			fh[ind1[indaj]] += 4*penF*(viol+2*Baij[ind2[a]+ind[j-a]]);
			for (int l=j+1;l<p;l++) if (iviol[l]) 
				fh[ind1[indaj]+ind[a]+l-a-indaj] += 8*penF*Baij[ind2[a]+ind[j-a]+l-j];
			for (int b=a+1;b<q;b++) {
				indbj = ind[b]+j-b;
				if (b<=j) fh[ind1[indaj]+indbj-indaj] += 8*penF*Babj[ind5[a]+ind4[b-a-1]+j-b];
				for (int l=b;l<p;l++)  if (iviol[l]) {
					indbl = ind[b]+l-b;
					if (indbl>indaj && l!=j) 
						fh[ind1[indaj]+indbl-indaj] += 8*penF*B[ind[a]+j-a]*B[ind[b]+l-b];
				}
			}
		}
	} }
	return(0);
}	

