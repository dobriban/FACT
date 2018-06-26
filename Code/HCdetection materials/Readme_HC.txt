HCdetection

This archive contains a Matlab implementation of HC-detection algorithm as described in the paper 
D. Donoho and J. Jin (2004) Higher criticism for detecting sparse heterogeneous mixtures



Current version 

V1.1


SHORT DOCUMENTATION

For more information, type help functionname in the Matlab command line.


Usage:
     [H, stat] = HCdetection(p, alpha, pvalcut)

Inputs:
  p        n-by-1 array of p-values from data

Inputs(Optional):
  alpha    0 < alpha < 1, indicates the smallest alpha*n p-values will be
            used to calculate the HC statistic
           default 1/2 
  pvalcut  0 < pvalcut < 1, indicates those small p-values (smaller than
            pvalcut) will be taken away to avoid heavy tails of test
            statistic
           default 1/n

Outputs
  H        "0" or "1" scalar, indicates whether H_0 is accepted ("0") or
           rejected ("1").
  stats    3-by-1 struct, including the test statistic
     stats.HCT        Higher Criticism test statistic
     stats.numselect  Number of variables at which HCT is achieved
     stats.HC         HC score for each feature


 
 
LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

If you use this code for your publication, please include a reference to the paper "Higher criticism for detecting sparse heterogeneous mixtures".
 
 
CONTACT

For any problem about the code or the paper, please contact

Jiashun Jin
Department of Statistics 
Carnegie Mellon University 
Email: jiashun@stat.cmu.edu 

Wanjie Wang 
Department of Statistics 
University of Pennsylvania 
Email: wanjiew@wharton.upenn.edu