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
## @deftypefn {Function File} {@var{retval} =} bg_sheath (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Guillermo Regodón <Billy@MacBook-Pro-de-Guillermo-3.local>
## Created: 2016-10-07










%
% Calculates the sheath in a cylindrical Langmuir probe according to radial
%  theory that includes the positive ion temperature and takes into account
%  the singularity of the speed of sound for the ion fluid
%
% First, the singularity is found by using that the system, when integrated
%  towards the plasma, has to remain bounded
% Then, the solution towards the probe is found
% If necessary, the cquaineutral solution is used to add more values
%  towards the plasma
%
% IMPORTANT: makeglobals has to be used prior to this rutine
% 
% This rutine uses bg_getInitFloat to obtain the points where to start RK
% It uses global variables, because the system function depends on parameter
%  and these parameters cannot be passed as arguments in the RK rutine
%
% It uses a RK modified for exit conditions particular for this system
%
% It uses bg_NnormEuler, that is the solution of Euler equation in this problem
%  bg_NnormEuler gives the higher or lower solution depending on global solEuler
%
% It uses bg_getxCI and bg_getydotCI, which are the quasineutral variables as
% a  function of the potential y
%
%  plotting and saving by means of globals __config__warn and __config__save
%
%
%
%
% xl is the maximum x position of the singularity, according to theory
% It is written by bg_getInitFloat
%
% x0a is where the x position of the singularity is stored
%
% __config__numprec es la precisión numérica de x0a
%
% xp is the probe radius to Debye length ratio
%










function S = bg_sheath ()
  
  global xl
  global kappa
  global xp
  global sol_Euler
  global beta
  global alpha0
  global gamma
  global Ip
  global x0a
  global vmas0
  
  global __config__numprec
  
  global __config__graph;
  global __config__save;
  global __config__warn;
  global __config__time_calculation_reduction
  global xMAX;
  _xMAX = xMAX + 20; % local of xMAX
  frac_tail = 0.5;
  
  xd = 0;
  xu = 0;
  xd_prev = 0;
  xu_prev = 0;
  
  _xp = xp; %local of xp
  
% Nrk is the number of rk steps for each x unit
% These values are very high (many steps)
  Nrk = 1000;
  if beta < 0.1 
    if Ip<=5 Nrk = 1000;
    elseif Ip<=14 Nrk = 4000;
    elseif Ip<=25 Nrk = 8000;
    %else Nrk = 12000;
    endif
  else
    if Ip <=60 Nrk = 250;
    elseif Ip<100 Nrk = 1000;
    else Nrk = 2000;
    endif
  endif

  % to reduce calculation time if such precission is unneeded
  Nrk = Nrk/__config__time_calculation_reduction;

  if __config__warn
    Nrk
    Ip
    beta
  endif
  
  xinf = 20;    % x units above x0 that be assume very high in plots
  
  if (kappa == 2)||(kappa == 3)
    
% for beta not equals 0, the singularity must be found
    if beta!=0
      sol_Euler = +1;
      
%Looks for the singularity
      if 1
        [x0 y0 ydot0 xl] = bg_getInitFloat;
        xd = 0;
        xu = xl;
        do
          
          xtest = 0.5*(xu+xd);     %search test variable
          
          [x0, y0, ydot0] = bg_getInitFloat(xtest);
          if !(__config__graph)
            _xp = x0;
          endif
          rkincrement = min([1 xl])/Nrk; % key line for calculation time
          
% bg_runge_kutta stops if non validity of the solution, valid only for this problem
          [y, x] = bg_runge_kutta(@bg_poisson_cyl, [y0 ydot0], [x0 xinf+x0], rkincrement );
          
          if y(end,1)>0
            xd_prev = xd;
            xu_prev = [];
            xd = xtest;
          else
            xd_prev = [];
            xu_prev = xu;
            xu = xtest;
          endif
          
          if __config__warn
            warning(['bg_sheath looking for singularity loop, xu - xd = ' num2str(xu-xd, '%g') '; xtest = ' num2str(xtest, '%f')]);
          endif
          
          if (xtest<0.001)
            warning(['Ip too low for sheath to form, Ip = ' num2str(Ip, '%g')])
            xu = xd;
            return ;
          endif
          
        until (abs(xu-xd)<max(__config__numprec*0.001, 1e-12))
        
      endif
      
      x0a = xtest;
      
      if __config__graph||__config__save
        _xp = max(_xMAX,4*x0a);
      else
        _xp = xp;
      endif
      
% if xp is to the right, to the plasma, more values have to be calculated
      if (length(y(:,1))>1) % &&(nargin==2)
        
        if length(xd_prev)==0
          xu = xu_prev;
        elseif length(xu_prev)==0
          xd = xd_prev;
        endif
        
% We have to calculate where the solutions u and d are too distint
        [x0u, y0u, ydot0u] = bg_getInitFloat(xu);
        [y_d, x_d] = bg_runge_kutta(@bg_poisson_cyl, [y0u ydot0u], [x0u xinf+x0u], rkincrement );
        [x0d, y0d, ydot0d] = bg_getInitFloat(xd);
        [y_u, x_u] = bg_runge_kutta(@bg_poisson_cyl, [y0d ydot0d], [x0d xinf+x0d], rkincrement );
        l_y_dif = min(length(y_u(:,1)),length(y_d(:,1)));
        y_dif = y_u(1:l_y_dif,1)-y_d(1:l_y_dif,1);
        for i = 2:length(y_dif)
          if y_dif(i)>1e-8 break endif
        endfor
        
        do
          
          x0w = x(i);
          y0w = y(i,1);
          
          ydotd = y(i,2)-0.01;
          ydotu = y(i,2)+0.01;
          
          x = x(1:(i-1));
          y = y(1:(i-1),:);
          
          do
            ydotwin = 0.5*(ydotu+ydotd);
            [ywin, xwin] = bg_runge_kutta(@bg_poisson_cyl, [y0w ydotwin], [x0w xinf+x0w], rkincrement );
            
            if ywin(end,1)<0
              ydotd = ydotwin;
            else
              ydotu = ydotwin;
            endif
          until abs(ydotu-ydotd)<__config__numprec
          
          if __config__warn
            warning(['quasineutral search 0, xwin = ' num2str(xwin(1)) ' , i = ' num2str(i)])
          endif
          
% Calculate the index where the solutions u and d are too distint
          [y_u, x_u] = bg_runge_kutta(@bg_poisson_cyl, [y0w ydotu], [x0w xinf+x0w], rkincrement );
          [y_d, x_d] = bg_runge_kutta(@bg_poisson_cyl, [y0w ydotd], [x0w xinf+x0w], rkincrement );
          l_y_dif = min(length(y_u(:,1)),length(y_d(:,1)));
          y_dif = y_u(1:l_y_dif,1)-y_d(1:l_y_dif,1);
          for j = 2:length(y_dif)
            if y_dif(j)>1e-7 break endif
          endfor

          %join
          i = i + j;
          x = [x; xwin];
          y = [y; ywin];
          
        until x(i)>_xp
        
      endif
      
      i = floor(findindexv(x, _xp));
      x = x(1:i);
      y = y(1:i,:);
      
      for i = 1:length(x)
        N(i) = bg_NnormEuler(Ip, x(i), y(i), beta, alpha0);
      endfor
      
      
      
      
      
      
% Calculates towards the axis
      _xp = 0;
      
      if (_xp<x0+rkincrement) %calculates to the left if necessary
        sol_Euler = -1;
        
        % this line is key to control de duration of the calculations
        rkincrement = min([1 xl])/Nrk;
        
        if rkincrement<x0-_xp
          
          x_up = _xp-rkincrement;
          if x_up<0
            x_up = 0.1*__config__numprec;
          endif
          
          [y2, x2] = bg_runge_kutta(@bg_poisson_cyl, [y0 ydot0], [x0 x_up], -rkincrement );
          
          for i = 1:length(x2)
            N2(i) = bg_NnormEuler(Ip, x2(i), y2(i), beta, alpha0);
          endfor
          
          
          
          x_array = [x2(end:-1:2); x];
          y_array = [y2(end:-1:2,1); y(:, 1)];
          z_array = [y2(end:-1:2,2); y(:, 2)];
          N_array = [N2(end:-1:2)'; N'];
        else
          y2 = [];
          x2 = [];
          
        endif
        
        if length(x2)!=0
          x_array = [x2(end:-1:2); x];
          y_array = [y2(end:-1:2,1); y(:, 1)];
          z_array = [y2(end:-1:2,2); y(:, 2)];
          N_array = [N2(end:-1:2)'; N'];
        else
          x_array = x;
          y_array = y(:,1);
          z_array = y(:,2);
          N_array = N';
        endif
      
      else
        x_array = x;
        y_array = y(:, 1);
        z_array = y(:, 2);
        N_array = N';
      endif
      
      
      
      
      
      
      
% to calculate the sheath for beta equals 0, we follow the ABR squeme
    else %beta==0
      
      if __config__graph
        _xp = 0;
      else
        _xp = xp;
      endif
      
      beta_0 = 4;
      [x0 y0 ydot0] = bg_getInitFloat(beta_0);
      
% if Ip is high, the calculation is long, we try to accelerate without
%  losing precission
      rkincrement = 0.005*(max(200, x0))/Nrk;
      
      [y2, x2] = bg_runge_kutta(@bg_poisson_cyl, [y0 ydot0], [x0 _xp-rkincrement], -rkincrement );
      
      for i = 1:length(x2)
        N2(i) = bg_NnormEuler(Ip, x2(i), y2(i), beta, alpha0); %el vector N2 queda horizontal
      endfor
      
      
      x_array = x2(end:-1:1);
      y_array = y2(end:-1:1,1);
      z_array = y2(end:-1:1,2);
      N_array = N2(end:-1:1)'; %for loop creates a horizontal array
            
    endif
    
  endif
  
  if __config__graph
    
    figure 1

    bg_plotlimitcond(1);
    %axis set by bg_plot_limit_cond
    if beta != 0
      plot(x_array, y_array, 'k', x_array, z_array, 'b', x_array, N_array, 'g', [x0a x0a], [-1 alpha0+2], 'color', [0.5 0.5 0.5])
    else %beta==0
      plot(x_array, y_array, 'k', x_array, z_array, 'b', x_array, N_array, 'g')
    end

%Other figures, V_profile
    figure 2   % V_profile, it is not saved, but can be computed easily
    hold off
    V_array = Ip./(x_array.*N_array);
    plot(x_array, V_array, 'r', x_array, sqrt((kappa*beta/2)*(N_array/(1+alpha0)).^(kappa-1)), 'color', [0.5 0 0])
    axis([0 2*x0a -0.5 alpha0+1.5])

%Other figures, electron and negative ion densities
    figure 3   % these are neither calculated, but easy to compute
    Ne = e.^(-y_array);
    Nmenos = alpha0*e.^(-gamma*y_array);
    plot(x_array, N_array, 'g', x_array, Ne, 'r', x_array, Nmenos, 'b', [x0a x0a], [-1 alpha0+2], 'color', [0.5 0.5 0.5])
    axis([0 2*x0a -0.5 alpha0+1.5])

    figure 1
    
  endif
  
  if (__config__save)
    data =  [x_array y_array z_array N_array];
    save('-ascii', ['bg_sheath_potential_density_profile_' num2str(beta, 'beta_%g_') num2str(alpha0, 'alpha0_%g_') num2str(gamma, 'gamma_%g_') num2str(kappa, 'kappa_%g_') num2str(Ip,'_Ip_%g_') num2str(x0a, 'x0a_%g.txt')], 'data')
    global I02
    global A
    global B
    global m_e
    I01 = I02*alpha0*sqrt(m_e/(B*gamma));
    data = [x_array Ip*ones(length(x_array), 1) (x_array).*(I01*e.^(-gamma*y_array) + I02*e.^(-y_array))];
    save('-ascii', ['bg_sheath_currents_profile_' num2str(beta, 'beta_%g_') num2str(alpha0, 'alpha0_%g_') num2str(gamma, 'gamma_%g_') num2str(kappa, 'kappa_%g_') num2str(Ip,'_Ip_%g.txt')], 'data')
  endif
  S = [x_array y_array z_array N_array];
  
endfunction
