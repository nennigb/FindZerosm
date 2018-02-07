function [K, Ci] = FindZerosmCourrone(R,N,fhandle,Ci,Refine,R0)

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
%  fhandle : finds the zero of the anonymous function fhandle
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


% check input args
% Si seulement 4 argugments, pas de rafinement et cercle centré en O
if nargin<4
    error('Missing argument...\n')
elseif nargin<=4
    R0 = 0;    
    Refine = 0;
elseif nargin==5
    R0 = 0;
end


tic
fprintf('Recherche de z�ros par CIM...\n')
ARRET = 1;

while ARRET~= 0
    % calcul des Sm
    %load('integrande')
    
    Theta = 2*pi*(0:N)/N;
    Thetai = 2*pi*(0:Ci.Ni)/Ci.Ni;

    % calcul de f'/f

    Z = R*exp(1i*Theta) + R0 ;
    fp_f = zeros(1,N+1);
    ff = fp_f;

    for ii = 1:(N+1)
        ff(ii) = fhandle(Z(ii));
    end
    fp_f = diffZcircleTheta(ff,Z,9,R0)./ff;



    %    save integrande fp_f Z
    % Int�gartion Nzero (- Npoles)
    S0 = ( 1/(2*1i*pi) )*(trapz(Z,fp_f) -trapz(Ci.Zi,Ci.fp_fi));
    fprintf('  > Le nombres de z�ros est %i\n', S0)
    if S0 >20
        fprintf('    -> Attention, plus de 20 z�ros, la m�thode peut �tre mal conditionn�e\n')
        fprintf('       utiliser plusieurs courrones...\n')
    end

    % crit�re d'arret de la boucle

    S0cut100 =round(10*S0)/10;
    if S0cut100== real(round(S0cut100))
        ARRET = 0;
    else
        ARRET = ARRET+1;
        R = R*Rshift;
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
% s0 = S0
S0 = round(S0);
% if imag(S0) ~= 0
%     imag(S0)
% end
S = zeros(1,S0); %p = zeros(1,S0+1);
% calcul des Sm


if R0==0
    % On factorise R^n (quand R est grand ca am�liore...)
    for ii=1:S0
        S(ii) = ( 1/(2*1i*pi) )*( (1i*R^(ii+1))*trapz(Theta, exp(1i*(ii+1)*Theta) .*fp_f) ...
            - (1i*Ci.Ri^(ii+1))*trapz(Thetai, exp(1i*(ii+1)*Thetai) .*Ci.fp_fi) );
    end
else
    % on ne factorise pas Rn (si R0 =/=0)
    for ii=1:S0
        S(ii) = ( 1/(2*1i*pi) )*( trapz(Z,fp_f) - trapz(Ci.Zi, Ci.fp_fi) );
    end
end

% calcul des coeff de p
C =  CoefC(S0,S);
p = C(end:-1:1);
% racines
K = roots(p);
% on trace les racines dans le plan complexe
fprintf('  > '); toc


hold all
% plot(Z)


% ------------------------------------------------------------------------%
% Rafinement autour des zeros
% ------------------------------------------------------------------------%
if Refine == 1
    for ii=1:S0
        Nr = 5000;
        Rr = abs(abs(K(ii)))*.05;
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



% ---------------------------------------------------------------------
% Update and put all the usefull variables in a struct interal Contour struct Ci
Ci.fp_fi = fp_f;
Ci.Zi = Z;
Ci.Ni = N;
Ci.Ri = R;
% save integrande fp_fi Zi Ni Ri
% ---------------------------------------------------------------------




% ------------------------------------------------------------------------%
% trac� des z�ros
% ------------------------------------------------------------------------%
% plot(real(K),imag(K), 'ob','MarkerSize',6)
% plot(Z)
% grid on