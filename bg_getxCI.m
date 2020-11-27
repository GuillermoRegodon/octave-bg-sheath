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
## @deftypefn {Function File} {@var{retval} =} getxCI (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-08-30









%
% Returns the x of the quasineutral solution
%










function [x] = bg_getxCI (y)

  global Ip
  global gamma
  global alpha0
  global beta
  global kappa

  c = e^(-y)+ alpha0*e^(-gamma*y);
  
  if (kappa == 2)
    x = Ip/(c*sqrt(y-2*beta*(c/(1+alpha0)-1)));
    
  elseif (kappa == 3)
    x = Ip/(c*sqrt(y-1.5*beta*(c^2/(1+alpha0)^2 - 1)));;
    
  endif

endfunction
