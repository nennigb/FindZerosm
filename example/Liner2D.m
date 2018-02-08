% ========================================================================%
%  Compute the guided mode in a 2D lined acoustic duct
%  Z - Air - rigid wall
%  validation with Poernomo:2008
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

%Computation frequency step   ------------------------------------------------
fmin = 100;
fmax = 7/2/pi;
Nfcalc =1;

% FindZerom parameters
Nrouche = 1000;          % Nb point in the intergration path
Rrouche = 10;             % Radius of the integration path
Refine = 0;	          % Refine by making small circle close to each root

%  data from Poernomo:2008
p.h=1
p.Z=(1+1i)*3.5;
c    = 1;

% frequency vector
fcalc = linspace(fmin,fmax,Nfcalc);

% Frequency sweep
for ii = 1:Nfcalc
    % frequency dependant part
    f = fcalc(ii);
    omega = 2*pi*f;
    k0 = (omega/c);
    p.k=k0;    
    
    % compute transverse wavenumber
    % Cicular contour, use external function DetMat2DZ.m
    % the info of the previous contour are put in Ci
    [K, Ci] = FindZerosm(Rrouche,Nrouche,@(alpha)DetMat2DZ(alpha,p),Refine,0);
    % uncomment this line to use several concentic contour
    % Annular
    [Ktemp, Ci] = FindZerosmAnnular(2*Rrouche,Nrouche,@(alpha)DetMat2DZ(alpha,p),Ci,Refine,0);
    K = [K ; Ktemp];
    % Ellipse
    [Ktemp, Ci] = FindZerosmAnnularEllipse([3*Rrouche 2*Rrouche],Nrouche,@(alpha)DetMat2DZ(alpha,p),Ci,Refine,0);
    K = [K ; Ktemp];
    % Ellipse
    [Ktemp, Ci] = FindZerosmAnnularEllipse([4*Rrouche 3*Rrouche],Nrouche,@(alpha)DetMat2DZ(alpha,p),Ci,Refine,0);
    K = [K ; Ktemp];
    
    % Plotting                 ------------------------------------------------
    figure(1)
    plot(real(K),imag(K),'k.')
    axis equal
    grid on
    hold on
    title('\alpha_n')
    xlabel('\Re \alpha')
    xlabel('\Im \alpha')
    % Assymptotique values are n*pi
    % plot((1:10)*pi, zeros(1,10),'ko')
end






