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
## @deftypefn {Function File} {@var{retval} =} findvaluev (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-10-09

%
% Returns the image of x1 in function y(x) sampled in x and y
%

function [y1] = findvaluev (x, y, x1)
  
  for i = 2:length(x)
    if (x(i)>x1)||(i==length(x))
      y1 = y(i-1) + (x1-x(i-1))*(y(i)-y(i-1))/(x(i) - x(i-1));
      break
    endif
  endfor
  
endfunction
