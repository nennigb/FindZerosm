function DET = DetMatAPP(alpha,p)
% compute the transverse wavenumber from a lined 2D acoustic duct of heigh
% p.h, with the impedance p.Z with the free field wavenumber p.k
% validation with Poernomo:2008, exp(- i omega t) time convention
% B. Nennig 

DET = alpha * sin (alpha*p.h) + (1i* p.k)/p.Z * cos (alpha *p.h);






