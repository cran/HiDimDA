/*
### RInterface.cpp  (2011-04-20)
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
#include "RInterface.h"
#include "FctqLDA.h"

extern "C"
SEXP RFnDist(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    vector<int> ind;
    double *cres;;
    ind.resize(nvar); 

    PROTECT(res = allocVector(REALSXP,1));
    cres = REAL(res);
    
    cres[0] = FnDist(REAL(B),INTEGER(p),INTEGER(q),REAL(Sigma),ind,REAL(k0),REAL(penF));
 
    PROTECT(dim = allocVector(INTSXP,1));
    INTEGER(dim)[0] = 1;
    SET_DIM(res,dim); 
    UNPROTECT(2);

    return(res);
}

extern "C"
SEXP Rfgrad(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    int nfact = INTEGER(q)[0];	
    int npar = nfact*(nfact+1)/2 + (nvar-nfact)*nfact;
    vector<int> iviol,ind;
    vector<double> Bij,fg;
    double *cres;

    iviol.resize(nvar); 
    ind.resize(npar); 
    fg.resize(npar); 
    Bij.resize(nvar*(nvar+1)/2); 
    
    fgrad(REAL(B),Bij,iviol,fg,INTEGER(p),INTEGER(q),REAL(Sigma),ind,REAL(k0),REAL(penF));
 
    PROTECT(res = allocVector(REALSXP,npar));
    cres = REAL(res);
    for(int i=0;i<npar;i++) cres[i]= fg[i];
    PROTECT(dim = allocVector(INTSXP,1));
    INTEGER(dim)[0] = npar;
    SET_DIM(res,dim); 
    UNPROTECT(2);

    return(res);
}

SEXP Rfhess(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    int nfact = INTEGER(q)[0];	
    int npar = nfact*(nfact+1)/2 + (nvar-nfact)*nfact;
    double *cres;
    vector<int> iviol,ind,ind1,ind2,ind3,ind4,ind5;
    vector<double> fh,Baij,Basq,Bij,Bjl,Bab,Babj;

    fh.resize(npar*(npar+1)/2); 
    Baij.resize(nfact*nvar*(nvar+1)/2); 
    Basq.resize(nfact); 
    Bij.resize(nvar*(nvar+1)/2); 
    Bjl.resize(nfact*(nfact+1)/2); 
    Bab.resize(nfact*(nfact+1)/2); 
    Babj.resize(nfact*nvar*nvar); 
    iviol.resize(nvar); 
    ind.resize(nvar); 
    ind1.resize(npar); 
    ind2.resize(nvar); 
    ind3.resize(nfact); 
    ind4.resize(nvar+1); 
    ind5.resize(nfact); 
     
    fhess(REAL(B),Baij,Basq,Bij,Bjl,Bab,Babj,iviol,fh,INTEGER(p),INTEGER(q),REAL(Sigma),ind,ind1,ind2,ind3,ind4,ind5,REAL(k0),REAL(penF));
 
   PROTECT(res = allocVector(REALSXP,npar*npar));
   cres = REAL(res);
   for(int i=0,fhind=0,r=0,c=0;i<npar;i++,fhind+=i)
	for(int j=0;j<=i;j++)  {
		cres[c*npar+r] = fh[fhind+j];
		if (r < npar-1) r++;
		else r = ++c;   
	}

   PROTECT(dim = allocVector(INTSXP,2));
   INTEGER(dim)[0] = npar;
   INTEGER(dim)[1] = npar;
   SET_DIM(res,dim); 
   UNPROTECT(2);

   return(res);
}

extern "C" 
SEXP RFnDist1(SEXP B,SEXP p,SEXP q,SEXP SigmSr,SEXP Srank,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    vector<int> ind;
    double *cres;;
    PROTECT(res = allocVector(REALSXP,1));
    cres = REAL(res);
    
    cres[0] = FnDist1(REAL(B),INTEGER(p),INTEGER(q),REAL(SigmSr),INTEGER(Srank),ind,REAL(k0),REAL(penF));
 
    PROTECT(dim = allocVector(INTSXP,1));
    INTEGER(dim)[0] = 1;
    SET_DIM(res,dim); 
    UNPROTECT(2);

    return(res);
}

extern "C"
SEXP Rfgrad1(SEXP B,SEXP p,SEXP q,SEXP SigmSr,SEXP Srank,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    int nfact = INTEGER(q)[0];	
    int npar = nfact*(nfact+1)/2 + (nvar-nfact)*nfact;
    double *cres;
    vector<int> iviol,ind;
    vector<double> Bij,fg;

    iviol.resize(nvar); 
    ind.resize(npar); 
    fg.resize(npar); 
    Bij.resize(nvar*(nvar+1)/2); 
    
    fgrad1(REAL(B),Bij,iviol,fg,INTEGER(p),INTEGER(q),REAL(SigmSr),INTEGER(Srank),ind,REAL(k0),REAL(penF));
 
    PROTECT(res = allocVector(REALSXP,npar));
    cres = REAL(res);
    for(int i=0;i<npar;i++) cres[i]= fg[i];
    PROTECT(dim = allocVector(INTSXP,1));
    INTEGER(dim)[0] = npar;
    SET_DIM(res,dim); 
    UNPROTECT(2);

    return(res);
}

extern "C"  
SEXP Rfhess1(SEXP B,SEXP p,SEXP q,SEXP SigmaSr,SEXP Srank,SEXP k0,SEXP penF)
{
    SEXP res,dim;
    int nvar = INTEGER(p)[0];	
    int nfact = INTEGER(q)[0];	
    int npar = nfact*(nfact+1)/2 + (nvar-nfact)*nfact;
    double *cres;
    vector<int> iviol,ind,ind1,ind2,ind3,ind4,ind5;
    vector<double> fh,Baij,Basq,Bij,Bjl,Bab,Babj;

    fh.resize(npar*(npar+1)/2); 
    Baij.resize(nfact*nvar*(nvar+1)/2); 
    Basq.resize(nfact); 
    Bij.resize(nvar*(nvar+1)/2); 
    Bjl.resize(nfact*(nfact+1)/2); 
    Bab.resize(nfact*(nfact+1)/2); 
    Babj.resize(nfact*nvar*nvar); 
    iviol.resize(nvar); 
    ind.resize(nvar); 
    ind1.resize(npar); 
    ind2.resize(nvar); 
    ind3.resize(nfact); 
    ind4.resize(nvar+1); 
    ind5.resize(nfact); 
     
    fhess1(REAL(B),Baij,Basq,Bij,Bjl,Bab,Babj,iviol,fh,INTEGER(p),INTEGER(q),REAL(SigmaSr),INTEGER(Srank),ind,
	ind1,ind2,ind3,ind4,ind5,REAL(k0),REAL(penF));
 
   PROTECT(res = allocVector(REALSXP,npar*npar));
   cres = REAL(res);
   for(int i=0,fhind=0,r=0,c=0;i<nvar;i++,fhind+=i)
	for(int j=0;j<=i;j++)  {
		cres[c*nvar+r] = fh[fhind+j];
		if (r < nvar-1) r++;
		else r = ++c;   
	}

   PROTECT(dim = allocVector(INTSXP,2));
   INTEGER(dim)[0] = npar;
   INTEGER(dim)[1] = npar;
   SET_DIM(res,dim); 
   UNPROTECT(2);

   return(res);
}
