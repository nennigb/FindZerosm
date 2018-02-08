% ========================================================================%
% Compute the roots of f(z) = cos(z)*P(z)
% where P(z) is a polynomial.
% the objective is to illustrated the adaptive ring splitting
% ========================================================================%

% This file is part of FindZerom, A package to compute the zeros of
% analytic functions Copyright (C) 2018  Benoit Nennig,
% benoit.nennig@supmeca.fr
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

close all
clear all



% definition of f, as an anonymous function
% very simple function, all roots are know !
f = @(z)  cos(z) * (z-1)*(z-2)*(z-3);


% =========================================================================
% Adaptive ring
% =========================================================================
NmaxPerCircle=15;       % number of root per ring before splitting it ~ [15 - 20]
NRootTarget = 50;      % we are looking for 100 zeros along Re direction
RootDistance = pi;      % average distance between 2 zeros


% FindZerom parameters
% You can change this value and see the effect on the error
Nrouche = 6000;                       % Nb point in the intergration path with max Radius
Rrouche = NRootTarget*RootDistance;   % Max Radius of the integration path
Refine = 0;	                          % Refine by making small circle close to each root
R0 = 0;

% compute the number of root inside the whool contour Rrouche
% not sensitive to round off, the main limitation is the exponential growth
% and the number of points Nrouche
[S0, R ] = NumberZerosm(Rrouche,Nrouche,@(z)f(z),0,R0);

% Split linearly the max radius into Ncircle
Ncircle = max([1, round(S0/NmaxPerCircle)]);        % Number of ring
Ri = round(max(R)/Ncircle);                         % ring radius

% compute the zeros inside each ring and collect it
K=[];

for nc=1:Ncircle
    if nc==1
        [Ktemp, Ci] = FindZerosm(Ri,round(Nrouche*(nc*Ri)/Rrouche),@(z)f(z),Refine,R0);
    else
        % Possible to adapt the number of integration points to the path perimeter
%         [Ktemp, Ci] = FindZerosmAnnular(nc*Ri,round(Nrouche*(nc*Ri)/Rrouche),@(z)f(z),Ci,Refine,R0);
        % the same number of integration point for all paths
        [Ktemp, Ci] = FindZerosmAnnular(nc*Ri,Nrouche,@(z)f(z),Ci,0,R0);        
    end
    K = [K ; Ktemp ];
    
    % -----------------------------------------------------------------
    % ploting all the steps
    figure(1)
    hold on
    % plot the roots
    plot(real(Ktemp),imag(Ktemp),'.')
    % plot the path
    plot(Ci.Zi)
    axis equal
    grid on
    
    title('Adaptative ring')
    xlabel('\Re z')
    ylabel('\Im z')
    % -----------------------------------------------------------------
    
    
end

% Error estimation
% as the root are real, the imaginary part leads to rougth an estimation
fprintf(' > Max error from Imaginary part = %g\n', max(abs(imag(K))))
fprintf(' > Average error from Imaginary part = %g\n', mean(abs(imag(K))))

%K*2/pi