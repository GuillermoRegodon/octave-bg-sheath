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
## @deftypefn {Function File} {@var{retval} =} getylimit (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@rabfis30.uco.es>
## Created: 2016-09-20










%
% Returns the turning point of the quasineutral solution, for each kappa value
%










function [retval] = bg_getylimit (a)

  global kappa
  if nargin >=1
    kappa = a;
  endif
  
  if ( kappa == 2 )
  
    retval = fsolve(@bg_limitk2, 0);
  
  elseif ( kappa == 3 )
  
    retval = fsolve(@bg_limitk3, 0);
  
  endif

endfunction
