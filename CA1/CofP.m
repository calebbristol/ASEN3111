function [C_p] = CofP(theta)
%% Coefficient of Pressure Calculator (CofP): Calculate Coefficient of Pressure
%   Uses equation for coefficient of pressure as defined in the lab
%   document to calculate value of coefficient of pressure at a given theta
%   value. Accepts vector inputs for theta.
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 08/31/21

%% Simplify Equation

% Temporary Variables to Reduce Clutter
v1 = 4 * (sin(theta) .^ 2);
v2 = 4 * sin(theta);
v3 = 1;

%These are simplified by the premise that Γ = 2πRV_inf

%v1 = 4sin^2(θ)
%v2 = 2Γsin(θ)/πRV_inf = 2(2)sin(θ) = 4sin(θ)
%v3 = (Γ/2πRV_inf)^2 = (1)^2 = 1

%% Calculate Coefficient of Pressure

C_p = 1 - (v1 + v2 + v3);

end

