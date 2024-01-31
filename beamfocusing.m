function [a] = beamfocusing(r, theta, Nt, d, f)
c = 3e8;
delta = (-(Nt-1)/2 : (Nt-1)/2)' * d;

% distance
r_n = sqrt(r^2 + delta.^2 - 2*r*delta*sin(theta));


% beamfocusing vector
a = exp( -1i * 2 * pi * f/c * r_n );

a_normalize = a/a(1);        
end