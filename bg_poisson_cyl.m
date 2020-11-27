## Copyright (C) 2016 Guillermo Regodón
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} poisson_cyl_bg (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-07-27









%
% Solves Poisson equation for this system, by means of the ion density obtained
%  using bg_NnormEuler
%










function [ydot C] = bg_poisson_cyl (y, x)   
  
  global beta
  global gamma
  global alpha0
  
  global Ip
  
  ydot(1) = y(2);
  
  C = bg_NnormEuler(Ip, x, y(1), beta, alpha0) - e^(-y(1)) - alpha0*e^(-gamma*y(1));
  ydot(2) = -y(2)/x + C;
  
endfunction
