function dfdz = diffZcircleTheta(f, Z, ordre,R0)
% > Différence finie centrée sur un cercle dfdz = dfdtheta * 1/(iZ)
% > utile pout FindZerosm.m

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

if nargin==3
    R=0;
end

Z =Z-R0;

dfdz = zeros(size(f));
h = abs(angle(Z(1)/Z(2)));

% les point f(1) et f(ebd) contiennent la m�me info... il faut donc d�caler

if ordre == 2
    % schema d'ordre 2 centr�

    dfdz(1) = f(2) - f(end-1);
    dfdz(end) = f(2) - f(end-1);
    dfdz(2:(end-1)) = f(3:end) - f(1:(end-2));


    dfdz = dfdz./(2*h * 1i*Z);

    
elseif ordre == 5
    %  schema centr� d'ordre 5
    dfdz(end) = -f(3) + 8*f(2) - 8*f(end-1) + f(end-2);
    dfdz(end-1) = -f(2) + 8*f(end) - 8*f(end-2) + f(end-3);
    dfdz(1) = -f(3) + 8*f(2) - 8*f(end-1) + f(end-2);
    dfdz(2) = -f(4) + 8*f(3) - 8*f(1) + f(end-1);

    dfdz(3:(end-2)) = - f(5:end) + 8*f(4:(end-1)) -8*f(2:(end-3)) + f(1:(end-4));

    dfdz = dfdz./(12*h * 1i*Z);
    
    
elseif ordre == 9
    %  schema centr� d'ordre 9
    a1 = 224/280;
    a2 = -56/280;
    a3 = 32/840;
    a4 = -1/280;
    dfdz(5:(end-4)) = -a4*f(1:(end-8)) -a3*f(2:(end-7)) -a2*f(3:(end-6)) -a1*f(4:(end-5)) + a1*f(6:(end-3)) +a2*f(7:(end-2)) +a3*f(8:(end-1)) + a4*f(9:end);
    dfdz(4) = -a4*f(end-1) -a3*f(1) -a2*f(2) -a1*f(3) + a1*f(5) +a2*f(6) +a3*f(7) + a4*f(8);
    dfdz(3) = -a4*f(end-2) -a3*f(end-1) -a2*f(1) -a1*f(2) + a1*f(4) +a2*f(5) +a3*f(6) + a4*f(7);
    dfdz(2) = -a4*f(end-3) -a3*f(end-2) -a2*f(end-1) -a1*f(1) + a1*f(3) +a2*f(4) +a3*f(5) + a4*f(6);
    dfdz(1) = -a4*f(end-4) -a3*f(end-3) -a2*f(end-2) -a1*f(end-1) + a1*f(2) +a2*f(3) +a3*f(4) + a4*f(5);

    dfdz(end) =     dfdz(1)  ;
    dfdz(end-1) = -a4*f(end-5) -a3*f(end-4) -a2*f(end-3) -a1*f(end-2) + a1*f(end) +a2*f(2) +a3*f(3) + a4*f(4);
    dfdz(end-2) = -a4*f(end-6) -a3*f(end-5) -a2*f(end-4) -a1*f(end-3) + a1*f(end-1) +a2*f(end) +a3*f(2) + a4*f(3);
    dfdz(end-3) = -a4*f(end-7) -a3*f(end-6) -a2*f(end-5) -a1*f(end-4) + a1*f(end-2) +a2*f(end-1) +a3*f(end) + a4*f(2);
    dfdz = dfdz./(h * 1i*Z);
end