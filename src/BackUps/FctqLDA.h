/*
### FctqLDA.h  (2011-04-20)
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

using namespace std;

double FnDist(double *B,int *p,int *q,double *Sigma,vector<int>& ind,double *k0,double *penF);
int fgrad(double *B,vector<double>& Bij,vector<int>& iviol,vector<double>& fg,int *p,int *q,double *Sigma,vector<int>& ind,double *k0,double *penF);
int fhess(double *B,vector<double>& Baij,vector<double>& Basq,vector<double>& Bij,vector<double>& Bjl,vector<double>& Bab,vector<double>& Babj,
		vector<int>& iviol,vector<double>& fh,int *p,int *q,double *Sigma,
		vector<int>& ind,vector<int>& ind1,vector<int>& ind2,vector<int>& ind3,vector<int>& ind4,vector<int>& ind5,double *k0,double *penF);
double FnDist1(double *B,int *p,int *q,double *SigmSr,int *Srank,vector<int>& ind,double *k0,double *penF);
int fgrad1(double *B,vector<double>& Bij,vector<int>& iviol,
		vector<double>& fg,int *p,int *q,double *SigmSr,int *Srank,vector<int>& ind,double *k0,double *penF);
int fhess1(double *B,vector<double>& Baij,vector<double>& Basq,vector<double>& Bij,vector<double>& Bjl,vector<double>& Bab,vector<double>& Babj,
		vector<int>& iviol,vector<double>& fh,int *p,int *q,double *SigmSr,int *Srank,
		vector<int>& ind,vector<int>& ind1,vector<int>& ind2,vector<int>& ind3,vector<int>& ind4,vector<int>& ind5,double *k0,double *penF);

