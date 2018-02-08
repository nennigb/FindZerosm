function [K, Ci] = FindZerosCouronneEllipse(R,N,fhandle,Ci,Refine,R0)
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
%  Ci : load contour usefull value for annular computation
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


% constant definition
RShift = 1.02;   % percent of increase of the radius when bad convegency
NRefine = 1000;       % number of integration points for the local refiement if Refine=1
RefineFraction = .05; % Refine Radius is RefineFraction*root value

if nargin<=3
    R0 = 0;
    Refine = 0;
end

if length(R)>1
    Contour = 'Ellipse';
end

tic
fprintf('Recherche de zéros par CIM Ellipse...\n')

ARRET = 1;

while ARRET~= 0
    % calcul des Sm
%     load('integrande')

    % contour exterieur
    if strcmp(Contour,'Ellipse')
        Theta = 2*pi*(0:N)/N;
        Z = R(1)*cos(Theta) +1i*R(2)*sin(Theta) + R0 ;
    else
        Theta = 2*pi*(0:N)/N;
        Z = R*exp(1i*Theta)+ R0  ;
    end
    % contour interieur
    Thetai = 2*pi*(0:Ci.Ni)/Ci.Ni;

    % calcul de f'/f
    % -----------------------------------------------------------------
    fp_f = zeros(1,N+1);
    ff = fp_f;

    for ii = 1:(N+1)
        ff(ii) = fhandle(Z(ii));
    end
    fp_f = diffZcircleEllipse(ff,Z,9,R0, Theta,R)./ff;
    %    save integrande fp_f Z
    % Int�gartion Nzero - Npoles
    S0 = ( 1/(2*1i*pi) )*(trapz(Z,fp_f) - trapz(Ci.Zi,Ci.fp_fi) );
    fprintf('  > Le nombres de zéros est %i\n', S0)
    if S0 >20
        fprintf('    -> Attention, plus de 20 zéros, la méthode peut être mal conditionnée\n')
        fprintf('       utiliser plusieurs courrones...\n')
    end

    % crit�re d'arret de la boucle
    S0cut100 = round(10*S0)/10;
    if S0cut100== round(S0cut100)
        ARRET = 0;
    else
        ARRET = ARRET+1;
        R = R*RShift;
        fprintf('    -> Second tour R * %g \n', RShift)
        if ARRET > 10
            disp('Erreur : Pas assez de points!!')
            return
        end
    end

end


% ------------------------------------------------------------------------%
% d�termination ddes coefs et polynome
% ------------------------------------------------------------------------%
S0 = round(S0);
S = zeros(1,S0); 

% calcul des Sm
for ii=1:S0
    S(ii) = ( 1/(2*1i*pi) )*(trapz(Z,(Z.^ii).*fp_f) - trapz(Ci.Zi,(Ci.Zi.^ii).*Ci.fp_fi) ) ;
end
% calcul des coeff de p
C =  CoefC(S0,S);
p = C(end:-1:1);
% racines
K = roots(p);
% display time
fprintf('  > '); toc


% ------------------------------------------------------------------------%
% Rafinement autour des zeros
% ------------------------------------------------------------------------%
if Refine == 1
    for ii=1:S0
        Nr = NRefine;
        Rr = abs(abs(K(ii)))*RefineFraction;
        Kr= FindZerosm(Rr,Nr,fhandle,0,K(ii));
        % s'il n'a rien trouv�
        while (isempty(Kr)==1)
            Rr = 2*Rr; Nr = 2*Nr;
            Kr= FindZerosm(Rr,Nr,fhandle,0,K(ii));
        end
        % attention si plusieurs zeros pris dans la tourmente...
        K(ii) = Kr(min(abs(Kr)) == abs(Kr));
    end
end



% Update and put all the usefull variables in a struct interal Contour struct Ci
Ci.fp_fi = fp_f;
Ci.Zi = Z;
Ci.Ni = N;
Ci.Ri = R;
% save integrande fp_fi Zi Ni Ri
% sauvegarde -------------- Couronne


% ------------------------------------------------------------------------%
% trac� des z�ros
% ------------------------------------------------------------------------%
% plot(real(K),imag(K), '.','MarkerSize',6)
% plot(Z)
% grid on
% axis equal
