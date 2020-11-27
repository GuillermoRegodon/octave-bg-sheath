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
## @deftypefn {Function File} {@var{retval} =} findrootv (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-10-09

function [result] = findcubicroots (a, b, c, d)
  
  % Trigonometric method
  p = (3*a*c - b^2)/(3*a^2);
  q = (2*b^3 - 9*a*b*c + 27*a^2*d)/(27*a^3);
  
  m=3*q*sqrt(-3/p)/(2*p);
  
  rm = real((1/3)*acos(m));
  im = imag((1/3)*acos(m));
  
  %if (abs(im) >= 0.001)
  %  warning(['Imaginary solution. imag(3*q*sqrt(-3/p)/(2*p)) = ' num2str(im, '%f')]);
  %endif
  
  n1=rm;
  n2=rm - 2*pi/3;
  n3=rm + 2*pi/3;
  
  t1 = 2*sqrt(-p/3)*cos( n1 );
  t2 = 2*sqrt(-p/3)*cos( n2 );
  t3 = 2*sqrt(-p/3)*cos( n3 );
  
  result = [t1 t2 t3] - b/(3*a);
  
endfunction
