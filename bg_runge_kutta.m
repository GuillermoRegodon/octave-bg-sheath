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
## @deftypefn {Function File} {@var{retval} =} runge_kutta (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-07-31










%
% Runge-Kutta algorithm with an exit option if the positive values diverge
%  or if the values become negative
%










function [y , t] = bg_runge_kutta (f, y_0, tlim, h)
  
  y(1,:) = y_0;
  
  t = [tlim(1):h:tlim(2)]';
  
  for i=2:length(t)
    
    k1 = f(y(i-1,:), t(i-1));
    k2 = f(y(i-1,:) + h/2*k1, t(i-1) + h/2);
    k3 = f(y(i-1,:) + h/2*k2, t(i-1) + h/2);
    k4 = f(y(i-1,:) + h*k3, t(i-1) + h);
    y(i,:) = y(i-1,:) + h/6*( k1 + 2*k2 + 2*k3 + k4 );
    %warning('vuelta en bucle dentro de RK') %optional loop warning
    if h>0
      %exit condition
      if ((i<length(t))&&(t(i+1)<=0))||(y(end,1)*y(end,2) > 0)
        break
      endif
    else
      global beta
      global Ip
      global alpha0
      global kappa
      if beta!=0
        N = (2*Ip^2*((1+alpha0)^(kappa-1))./(kappa*beta*t(i).^2)).^(1/(kappa+1));
        if (y(i,1) < (Ip./(t(i).*N)).^2 + (kappa/(kappa-1))*beta*((N/(1+alpha0)).^(kappa-1)-1))
          y = y(1:i-1,:);
          break
        endif
      endif
    endif
    
  endfor
  t = t(1:length(y(:,1)));
  
endfunction
