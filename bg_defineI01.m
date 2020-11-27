## Copyright (C) 2018 Guillermo Regodón
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
## @deftypefn {Function File} {@var{retval} =} bg_defineI01 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2018-03-21









%
% I02 is the coeficient for the current formula for the electrons
% I01 is the one for negative ions, that is calculated using I02
%










function [I01] = bg_defineI01 ()
  
  global alpha0
  global gamma
  global I02
  global m_e
  global B
  
  I01 = I02*alpha0*sqrt(m_e/(B*gamma));
  
endfunction
