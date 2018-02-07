function C = CoefC(S0,S)
% retroune directement le vecteur C avec tout les coef.


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

C = zeros(1, S0);
C(S0+1) = 1;

for n=1:S0
    kp = S0 - n +1;
    C(kp) = (1/(kp-1-S0)) *( S(1:n) * C( (kp+1):(S0+1) ).' ); 
end
