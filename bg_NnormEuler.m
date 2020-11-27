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
## @deftypefn {Function File} {@var{retval} =} NnormEuler (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-07-27










%
% Return the solution to the Euler equation when applied to this system
%  the branch is selected using global sol_Euler
%  kappa allows the selection of the right expressions from theory
%  show == 1 to show intermediate results
%










% Note that Ip, beta, alpha0 are not globals, so it can be used more easily in 
%  plotting of representative
function result = bg_NnormEuler (Ip, x, y, beta, alpha0, show)

  global sol_Euler
  global kappa
  global __config__cylindrical
  
  if nargin>=6
    if show == 1
      a = (kappa/(kappa-1))*beta/(1+alpha0)^(kappa-1)
      b = - y - (kappa/(kappa-1))*beta
      c = 0
      %d = -4*b^3/(27*a^2)
      if __config__cylindrical
        d = Ip^2/x^2
      else
        d = Ip^2
      endif
    endif
  else
      a = (kappa/(kappa-1))*beta/(1+alpha0)^(kappa-1);
      b = -y - (kappa/(kappa-1))*beta;
      c = 0;
      %d = -4*b^3/(27*a^2) 
      if __config__cylindrical
        d = Ip^2/x^2;
      else
        d = Ip^2;
      endif
  endif
  
  if (beta != 0)
    
    if kappa==2 %a*N^3 + b*N^2 ∫+ d = 0
      
      % trigonometric solution for cubic
      p = - b^2/(3*a^2);
      q = (2*b^3 + 27*a^2*d)/(27*a^3);
      
      m=3*q*sqrt(-3/p)/(2*p);
      
      rm = real((1/3)*acos(m));
      im = imag((1/3)*acos(m));
      
      if nargin>=6
        if show==1
          rm
          im
        endif
      endif
      
      if sol_Euler == -1
        n=rm- 2*pi/3; %Smaller positive solution
      elseif sol_Euler == +1
        n=rm; %Higher positive solution
      endif
      
      t0 = 2*sqrt(-p/3)*cos( n );
      
      result = t0 - b/(3*a);
      
    elseif kappa == 3 %a*N^4 + b*N^2 + d = 0
      
      if sol_Euler == -1
        result = real(sqrt((-b - sqrt(b^2-4*a*d))/(2*a))); %Smaller positive solution
      elseif sol_Euler == +1
        result = real(sqrt((-b + sqrt(b^2-4*a*d))/(2*a))); %Higher positive solution
      endif
      
    endif
      
  else %return the ABR quasineutral solution
  
    result = Ip/(x*sqrt(y));
  
  endif






endfunction
