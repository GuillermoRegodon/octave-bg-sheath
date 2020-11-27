## Copyright (C) 2017 Guillermo Regodón
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
## @deftypefn {Function File} {@var{retval} =} find2graderoots (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2017-04-28

% coef = [a b c], en ax^2`bx+c=0
function [retval] = find2graderoots (coef)
  delta2 = coef(2)^2-4*coef(1)*coef(3);
  delta = sqrt(delta2);
  retval(1) = (-coef(2) + delta)/(2*coef(1));
  retval(2) = (-coef(2) - delta)/(2*coef(1));
endfunction
