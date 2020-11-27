## Copyright (C) 2019 Guillermo Regodón
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
## @deftypefn {Function File} {@var{retval} =} findvaluei (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2019-05-22

%
% Returns the image of index i in vector y, even if i is not in [1:end]
%

function [yi] = findvaluei (y, i)
  
  if i<1
    yi = y(1) + (i-1)*(y(2)-y(1));
  elseif i>=length(y)
    yi = y(end) + (i-length(y))*(y(end)-y(end-1));
  else
    i_f = floor(i);
    yi = y(i_f)+(i-i_f)*(y(i_f+1)-y(i_f));
  endif
  
endfunction
