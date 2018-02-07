function dfdz = diffZcircleEllipse(f, Z, ordre,R0,Theta,R)
% Diff�rence finie centr�e dfdz = dfdtheta * 1/(-a*sin(theta) +1i*b*cos(theta)


%    pt a  delta/2   delta/2 pt b
%     (i-1) ----- (i) ------ (i+1)
%       -1         0           1     /(2h)

if nargin==3
    R=0;
end

Z =Z-R0;

dfdz = zeros(size(f));
h = Theta(2)-Theta(1);

Den = -R(1)*sin(Theta) + 1i*R(2)*cos(Theta);
% les point f(1) et f(ebd) contiennent la m�me info... il faut donc d�caler

if ordre == 2
    % schema d'ordre 2 centr�

    dfdz(1) = f(2) - f(end-1);
    dfdz(end) = f(2) - f(end-1);
    dfdz(2:(end-1)) = f(3:end) - f(1:(end-2));


   dfdz = dfdz./(2*h * Den);
   

elseif ordre == 5
    %     schema centr� d'ordre 5
    dfdz(end) = -f(3) + 8*f(2) - 8*f(end-1) + f(end-2);
    dfdz(end-1) = -f(2) + 8*f(end) - 8*f(end-2) + f(end-3);
    dfdz(1) = -f(3) + 8*f(2) - 8*f(end-1) + f(end-2);
    dfdz(2) = -f(4) + 8*f(3) - 8*f(1) + f(end-1);

    dfdz(3:(end-2)) = - f(5:end) + 8*f(4:(end-1)) -8*f(2:(end-3)) + f(1:(end-4));

    dfdz = dfdz./(12*h * Den);
elseif ordre == 9
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
    dfdz = dfdz./(h * Den);
end