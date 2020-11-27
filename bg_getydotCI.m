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
## @deftypefn {Function File} {@var{retval} =} getydotCI (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-08-30









%
% Returns the value of ydot, the derivative of y with respect to x
%










function [ydot0] = bg_getydotCI (x0, y0)

  global Ip
  global gamma
  global alpha0
  global beta
  global kappa
  
  if (kappa == 2)
    ydot0_1 = -x0*(bg_fNp(y0)/bg_fN(y0) + 0.5*(1-2*beta*bg_fNp(y0)/(1+alpha0))/(y0 - 2*beta*(bg_fN(y0)/(1+alpha0)-1)));
    ydot0 = 1/ydot0_1;
    
  elseif (kappa == 3)
    ydot0_1 = -x0*(bg_fNp(y0)/bg_fN(y0) + 0.5*(1 - 1.5*beta*(2*bg_fN(y0)/(1+alpha0)^2)*bg_fNp(y0))/(y0 - 1.5*beta*((bg_fN(y0))^2/(1+alpha0)^2 - 1)));
    ydot0 = 1/ydot0_1;
    
  endif
  
  
endfunction
