function [S0, R, Ci ] = NumberZerosm(R,N,fhandle,Refine,R0)

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


% ========================================================================%
% Package to compute zeros of analytic function
% ========================================================================%
% > Based on Cauchy Integration Method (CIM) or the Argument Principle Method (APM)
%   see the documentation for more information and references
% > B. Nennig, E. Perrey-Debain, and M. Ben Tahar. A mode matching method 
%   for modelling dissipative silencers lined with poroelastic materials 
%   and containing mean ow. J. Acoust. Soc. Am. , 128(6) :33083320, 2010.
% > C. Chen, P. Berini, D. Feng, S. Tanev, and V. Tzolov. Ecient and 
%   accurate numerical, analysis of multilayer planar optical waveguides 
%   in lossy anisotropic media. Opt. Express, 7(8) :260272, 2000.
% > Poles location is not included

% Mandatory input args :
%  R : integration radius
%  N : number of integration points
%  fhandle : is the anonymous function of which the root are sought
% Optionnal input args :
%  Refine = 1 ou 0 (local refinement with small circle around each root)
%  R0 = 0 si l'origine est en O, coordonée de l'origine sinon
% Mandatory output args :
%  K : roots
% Optionnal output args :
%  Ci : save of contour usefull value for annular computation
%==========================================================================
%                                      benoit.nennig@supmeca.fr     03/2009
%==========================================================================


% ------------------------------------------------------------------------%
% détermination du nombre de zéros et de pôles
% ------------------------------------------------------------------------%

% constant definition
RShift = 1.02;   % percent of increase of the radius when bad convegency

% check input args
if nargin<3
    error('Missing argument\n')
elseif nargin<=3
    R0 = 0;    
    Refine = 0;
elseif nargin<=4
    R0 = 0;
end

tic
fprintf('How many zeros inside... (R=%i)\n',R)
ARRET = 1;

while ARRET~= 0
    % calcul de f'/f
    Z = R*exp(2i*pi*(0:N)/N) + R0 ;
    fp_f = zeros(1,N+1);
    ff = fp_f;

    for ii = 1:(N+1)
        ff(ii) = fhandle(Z(ii));
    end
    fp_f = diffZcircleTheta(ff,Z,9,R0)./ff;


    % ---------------------------------------------------------------------
    % Sauvegarde de f'/f pour accélérer le calcul sur une couronne
    % put all the usefull variable in a struct interal Contour struct Ci
    Ci.fp_fi = fp_f;
    Ci.Zi = Z;
    Ci.Ni = N;
    Ci.Ri = R;
%     % Sauvegarde de f'/f pour acc�l�er le calcul sur une couronne
%     fp_fi = fp_f;
%     Zi = Z;
%     Ni = N;
%     Ri = R;
%     save integrande fp_fi Zi Ni Ri
    % ---------------------------------------------------------------------



    % Int�gartion Nzero (- Npoles)
    S0 = ( 1/(2*1i*pi) )*trapz(Z,fp_f);
    fprintf('  > there is %i zeros \n', S0)


    % crit�re d'arret de la boucle
    S0cut100 = round(10*S0)/10;
    if S0cut100== round(S0cut100)
        ARRET = 0;
    else
        ARRET = ARRET+1;
        R = R*RShift;
        fprintf('    -> Second tour R * %g \n', RShift)
        if ARRET > 10
            disp('Erreur : Pas assez de points.')
            return
        end
    end



end
