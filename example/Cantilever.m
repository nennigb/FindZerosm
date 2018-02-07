% ========================================================================%
% Compute the roots of a cantilever beam
% f(alphaL) = cos(alphaL) * cosh(alphaL) + 1
% see https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory
% ref value alpha_n.L = 1.875, 4.69, 7.85 and alpha_n.L ~ (2n-1)*pi/2 
% ========================================================================%

close all
clear all

% FindZerom parameter
Nrouche = 500;          % Nb point in the intergration path
Rrouche = 15;             % Radius of the integration path
Refine = 0;	          % Refine by making small circle close to each root

% definition of f, as an anonymous function
cantileverBeam = @(alphaL)  cos(alphaL) * cosh(alphaL) +1;


% =========================================================================
% First Approach 0 centred contour
% =========================================================================
% call of FindZerom function
K = FindZerosm(Rrouche,Nrouche,cantileverBeam,Refine,0);

% Keep only real and positive solution
AlphanL = sort(real(K(find(real(K) > 0 & abs(imag(K))< 1e-2))),'ascend');

% Plotting                 ------------------------------------------------
figure(1)
plot(real(K),imag(K),'.')
plot(real(AlphanL),imag(AlphanL),'o')
axis equal
grid on
hold on
title('\alpha_n L')
xlabel('\Re \alpha L')
ylabel('\Im \alpha L')




% =========================================================================
% Second Approach decentred contour
% more suitable here because only real positive root are wanted
% =========================================================================
% call of FindZerom function
Kd = FindZerosm(Rrouche,Nrouche,cantileverBeam,Refine,Rrouche/2);

% Keep only real and positive solution
AlphanLd = sort(real(Kd(find(real(Kd) > 0 & abs(imag(Kd))< 1e-2))),'ascend');

% Plotting                 ------------------------------------------------
figure(2)
plot(real(Kd),imag(Kd),'.')
plot(real(AlphanLd),imag(AlphanLd),'o')
axis equal
grid on
hold on
title('\alpha_n L')
xlabel('\Re \alpha L')
ylabel('\Im \alpha L')



