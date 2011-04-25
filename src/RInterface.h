/*
### RInterface.h  (2011-04-20)
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

#include <R.h>
#include <Rdefines.h>

extern "C" 
SEXP RFnDist(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF);
extern "C" 
SEXP Rfgrad(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF);
extern "C" 
SEXP Rfhess(SEXP B,SEXP p,SEXP q,SEXP Sigma,SEXP k0,SEXP penF);
extern "C" 
SEXP RFnDist1(SEXP B,SEXP p,SEXP q,SEXP SigmSr,SEXP Srank,SEXP k0,SEXP penF);
extern "C" 
SEXP Rfgrad1(SEXP B,SEXP p,SEXP q,SEXP SigmSr,SEXP Srank,SEXP k0,SEXP penF);
extern "C" 
SEXP Rfhess1(SEXP B,SEXP p,SEXP q,SEXP SigmaSr,SEXP Srank,SEXP k0,SEXP penF);

