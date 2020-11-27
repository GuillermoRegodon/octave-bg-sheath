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
## @deftypefn {Function File} {@var{retval} =} getInitFloat (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-08-01










%
% Returns the values of x, y and ydot (the derivative of y with respect to x)
%  for a given x0a in the singularity curve. It return xl as well, the maximum
%  value for x in the singularity curve.
%
% It is usually called only to get xl
% If x0a == 0, we will get three divisions over 0 that do not affect
%  the whole calculation
%










function [x0 y0 ydot0 xl] = bg_getInitFloat (x0a)

  global Ip
  global gamma
  global alpha0
  global beta
  global xp
  global xl
  global kappa
  global __config__cylindrical
  
  if 1 %beta!=0 has to be ensured before calling
    
    if nargin <1 
      option = 1;
    else 
      option = 4;
    endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if option == 4
    
      if __config__cylindrical
        
        if beta!=0
          
          if kappa==2
            x0 = x0a;
            y0 = 3*(beta*Ip/(x0*(1+alpha0)))^(2/3) - 2*beta;
            xl = (1.5)^(1.5)*Ip/(sqrt(beta)*(1+alpha0));
            ydot0 = -2*(y0 + 2*beta)/(3*x0);
          elseif kappa==3
            x0 = x0a;
            y0 = sqrt(6*beta)*(Ip/(x0*(1+alpha0))) - 1.5*beta;
            xl = sqrt(8/(3*beta))*Ip/(1+alpha0);
            ydot0 = -sqrt(6*beta)*Ip/(x0^2*(1+alpha0));
          endif
          N0 = bg_NnormEuler(Ip, x0, y0, beta, alpha0);
        else %beta == 0
          y0 = Ip*0.0003;
          N0 = e^(-y0) + alpha0*e^(-gamma*y0);
          x0 = Ip/(sqrt(y0)*N0);
          ydot0 = 1/(x0*((e^(-y0) + alpha0*gamma*e^(-gamma*y0))/N0 - 0.5/y0));
          xl = x0;
          if x0>40
            warning(['x0 very high, it may take time, x0 = ' num2str(x0, '%g')]);
          endif
        endif
        
      else
        
        warning('Using rutines for cylindrical with __config__cylindrical != 1');
        x0 = 0;
        y0 = 0;
        ydot0 = 0;
        xl = 0;
        
      endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif option == 1
      if beta!=0
        if kappa==2
          xl = (1.5)^(1.5)*Ip/(sqrt(beta)*(1+alpha0));
          
        elseif kappa==3
          xl = sqrt(8/(3*beta))*Ip/(1+alpha0);
          
        endif
        
        y0asin = 1.1*bg_getylimit(kappa);
  
        x0 = bg_getxCI(y0asin);
        
        y0 = y0asin;
        ydot0 = bg_getydotCI(x0, y0);
        
      else %beta == 0
        y0 = Ip*0.0003/(1+alpha0);
        N0 = e^(-y0) + alpha0*e^(-gamma*y0);
        x0 = Ip/(sqrt(y0)*N0);
        ydot0 = 1/(x0*((e^(-y0) + alpha0*gamma*e^(-gamma*y0))/N0 - 0.5/y0));
        xl = x0;
      endif

    endif
    
  endif  
  
endfunction
