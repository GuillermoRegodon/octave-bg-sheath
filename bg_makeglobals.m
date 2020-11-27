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
## @deftypefn {Function File} {@var{retval} =} makeglobals (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-07-27










%
% Definition of global for the rutines in this solution
%
% globals had to be used, as Runge-Kutta and other rutines used call to other rutines
%  with a series of arguments to obtain values for which the expresions depended on
%  the same parameters Ip, beta, alpha0, gamma...
% Some struct could have been used, but simplicity was chosen above
%  other considerations.
%
% bg stands for beta_gamma, the main global constants of this solution
%










% Sheath
global kB %Boltzmann constant
kB = 1.380648e-23; %SI

global k %k_B/e
k = 8.61732739e-5; %SI

global m_e %electron mass
m_e = 5.4858e-4; %umas

global uma
uma = 1.66054e-27; %SI

global A %positive ion mass
global B %negative ion mass
global beta %positive ion temperature to electron temperature ratio
global gamma %electron temperature to negative ion temperature ratio
global alpha0 %negative ion density to electron density ratio
global kappa %thermodynamic adiabatic index
global vmas0 %Bohm velocity
global xp %probe raduis to Debye length ratio

global Ip %adimensional ion current
global I02 %electron current constant
global I01 %negative ion current constant
global yp %adimensional probe potential
global yf %adimensional floating potential


global xl %maximum x for singularity, according to theory
global x0a %singularity x coordinate for  bounded quasineutral solution
global xMAX %Maximum integration limit if singularity x coordinate is small

global sol_Euler %index to select Euler's equation branch for the solution


global __config__cylindrical
global __config__numprec
global __config__graph
global __config__save
global __config__warn
global __config__time_calculation_reduction %reduces compute time, 3 is ok

global q_e
q_e = 1.602176e-19; %C
global Ap
Ap = pi*6e-3*2e-4;

warning('off', 'backtrace')
warning('Please initialize: Ip, beta, alpha0, gamma, kappa, xp')
warning('I02, __config__numprec, __config__cylindrical')

page_output_immediately(1)
page_screen_output(0)

Ip = 5
beta = 0.1
alpha0 = 0.5
gamma = 10
kappa = 2

xp = 1
xMAX = 50

% Argon and electrons
A = 39.948
B = m_e

%Neaon and electrons
%A = 20.18
%B = m_e

%Oxigen O+ and O-
%A = 15.9989
%B = 15.9999

I02 = bg_defineI02
I01 = bg_defineI01
sol_Euler = 1

__config__cylindrical = 1
__config__numprec = 1e-9
__config__graph = 1
__config__save = 0
__config__warn = 1
__config__time_calculation_reduction = 3
