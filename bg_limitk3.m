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
## @deftypefn {Function File} {@var{retval} =} limitk3 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@rabfis30.uco.es>
## Created: 2016-09-20










%
% Derivative of x with respect to y in the quasineutral solution
% valid for kappa == 3
%










function [f_y] = bg_limitk3 (y)
% prior to use, kappa must be ensured to be kappa == 3, also beta!=0 and y!=0

  global alpha0
  global gamma
  global beta
  
  f_y = (e.^(2*y)*(1+alpha0)^2)./(( 1 + alpha0*e.^(-y*(gamma-1)) ).^2) - 3*beta./(y + 1.5*beta);

endfunction
