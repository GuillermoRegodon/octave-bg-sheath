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
## @deftypefn {Function File} {@var{retval} =} findquarticroots (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@rabfis30.uco.es>
## Created: 2016-11-10

function [result] = findquarticroots (a, b, c, d, e)
  
  p = (8*a*c - 3*b^2)/(8*a^2);
  q = (b^3 - 4*a*b*c + 8*(a^2)*d)/(8*a^3);
  
  D0 = c^2 - 3*b*d + 12*a*e;
  D1 = 2*c^3 - 9*b*c*d + 27*(b^2)*e + 27*a*d^2 - 72*a*c*e;
  
  Q = (0.5*(D1 + sqrt(D1^2 - 4*D0^3)))^(1/3);
  S = 0.5*sqrt(-2*p/3 + (Q + D0/Q)/(3*a));
  
  result(1) = abs(-b/(4*a) + S + 0.5*sqrt( -q/S - 4*S^2 - 2*p ));
  result(2) = abs(-b/(4*a) + S - 0.5*sqrt( -q/S - 4*S^2 - 2*p ));
  result(3) = abs(-b/(4*a) - S + 0.5*sqrt(  q/S - 4*S^2 - 2*p ));
  result(4) = abs(-b/(4*a) - S - 0.5*sqrt(  q/S - 4*S^2 - 2*p ));
  
endfunction
