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
## @deftypefn {Function File} {@var{retval} =} plotlimitcond (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-09-24










%
% PLots the singularity curves for the potential, electricl field, ion density
%  and ion velocity
% Plots the asimptotic curve, with its chacracteristic return to the plasma
% After initializing (see bg_sheath.m), use:
%
%    S = bg_sheath;
%    x_array = S(:,1);
%    y_array = S(:,2);
%    z_array = S(:,3);
%    N_array = S(:,4);
%    plot(x_array, y_array, x_array, z_array, x_array, N_array)
%
% to check that the potential solution is tangent to the
%   singularity potential curve
%
% new argument for a new plot in the actual figure
% hold is set to on at the end of the rutine
%










function [retval] = bg_plotlimitcond (new)

  global Ip
  global beta
  global alpha0
  global kappa
  global x0a
  
%  figure 1
  
  if nargin >=1
    if new
      hold off
    else
      hold on
    endif
  endif
  
  y0asin = bg_getylimit(kappa);
  xc = bg_getxCI(y0asin);
  
  global xl
  
  if beta != 0
    if kappa==2
      xl = (1.5)^(1.5)*Ip/(sqrt(beta)*(1+alpha0));
    elseif kappa==3
      xl = sqrt(8/(3*beta))*Ip/(1+alpha0);
    endif
  else
    xl = 300;
  endif
  
  x = [xl/1000:xl/1000:xl];
  
  N = (2*Ip^2*((1+alpha0)^(kappa-1))./(kappa*beta*x.^2)).^(1/(kappa+1));

% plot ion density singularity curve
  plot(x, N, 'color', [0.5 1.0 0.5])
% velocity singularity curve, if required
%  plot(x, Ip./(x.*N), 'color', [0.5 0.5 0.5] );
  hold on
  

  z = -2*(Ip^2)./(x.^3.*N.^2);
  plot(x, z, 'color', [0.5 0.5 1.0])
  
  y = (Ip./(x.*N)).^2 + (kappa/(kappa-1))*beta*((N/(1+alpha0)).^(kappa-1)-1);
  plot(x, y, 'color', [0.5 0.5 0.5] , [xl xl], [-2 y(1)], 'color', [0.5 0.5 0.5], [xc xc], [-2 y(1)], 'color', [0.5 0.5 0.5]);
  
  sol_Euler = +1;
  for i=1:length(x)
    Nmas(i) = bg_NnormEuler(Ip, x(i), y(i), beta, alpha0);
  endfor
  
  sol_Euler = -1;
  for i=1:length(x)
    Nmenos(i) = bg_NnormEuler(Ip, x(i), y(i), beta, alpha0);
  endfor
  
  % asymptotic solution. It has a turn in x, so it is better to set y first
  yasin_array = [0.001:0.0001: 2];
  xasin_array = zeros(1, length(yasin_array));
  for i=1:length(yasin_array)
    xasin_array(i) = bg_getxCI(yasin_array(i));
  endfor
  plot(xasin_array, yasin_array, 'color', [1 0.6 0.8]);
  
  if y0asin != 0
    axis([0 xl*1.1 -20*y0asin 40*y0asin])
  else % y0asin == 0 %it happens if beta equals 0, see bg_limitk2
    axis([0 2*x0a -1 (alpha0+1.5)])
  end
  grid
    
endfunction
