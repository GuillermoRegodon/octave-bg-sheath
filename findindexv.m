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
## @deftypefn {Function File} {@var{retval} =} findindexv (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-10-23

%
% IMPORTANT: this function is only applicable if y is increasing
% i is the first index for which y(i) is greater or equal to yl
%
% i_frac returns the fractional index of yl, even if i is not in [1:end]
%

function [i, i_frac] = findindexv (y, yl)
  
  test_incr = y(2:end)-y(1:end-1); 
  test_incr = (test_incr>0);
  incr = sum(test_incr==0)==0;
  decr = sum(test_incr==1)==0;
  
  for i = 2:length(y)
    if (y(i)>=yl)||(i==length(y))
      break
    endif
  endfor
  
  if incr
    i_int = length(y)-sum(y>yl);
  elseif decr
    i_int = sum(y>yl);
  endif
  
  if (incr||decr)
    if i_int<=1
      i_frac=1+(yl-y(1))/(y(2)-y(1));
    elseif i_int<length(y)
      i_frac=i_int+(yl-y(i_int))/(y(i_int+1)-y(i_int));
    elseif i_int==length(y)
      i_frac = i_int-1+(yl-y(end-1))/(y(end)-y(end-1));
    endif
  endif
  
endfunction
