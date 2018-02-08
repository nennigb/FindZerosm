function DET = DetMatAPP(alpha,p)
% compute the transverse wavenumber from a lined 2D acoustic duct of heigh
% p.h, with the impedance p.Z with the free field wavenumber p.k
% validation with Poernomo:2008, exp(- i omega t) time convention

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

DET = alpha * sin (alpha*p.h) + (1i* p.k)/p.Z * cos (alpha *p.h);






