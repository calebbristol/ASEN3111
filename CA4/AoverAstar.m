function [Msup, Msub] = AoverAstar(gamma,AoAstar)
% Solves A/A* [Eq. 10.32] Function begins rootsolve at M>1 to ensure
% supersonic flow at exit
options = optimset('Display','off');
star = @(x) ((1/x^2)*((2/(gamma+1))*(1+((gamma-1)/2)*x^2))^((gamma+1)/(gamma-1)))-AoAstar^2;
Msup = fsolve(star,1.01,options);
Msub = fsolve(star,.1,options);
end