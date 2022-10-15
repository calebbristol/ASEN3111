function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%% Function PLLT: Solves the Fundamental Equation for Prandtl's Lifting Line Theory
%
% Utilizes N discrete points on the wing along with N odd terms from PLLT
% to create a system of equations which solves for the A_n variables. These
% are then used to calculate e, c_L, and c_Di

%% Create Discrete Points
i = 1:N;
theta = (i*pi)./(2*N);

%% Analyze Geometry of Airfoil
%
% Create functions to calculate the chord length and lift slope at a given
% point on the airfoil. Done by converting theta values into values of y
% and then assuming linear variation in the wing geometry.

% Aspect Ratio
%
% Assumes Wing to be trapezoidal
S = (b/2) * (c_t + c_r);
AR = b^2 / S;

% Function for converting to y value
y = @(theta) b/2 * cos(theta);

% Functions for a(y) and c(y), absolute value sign to assure wing symmetry
a = @(y) (a0_t - a0_r) * abs(y / (b/2)) + a0_r;

c = @(y) (c_t - c_r) * abs(y / (b/2)) + c_r;

% Zero Lift Angle of Attack
alpha_0 = @(y) (aero_t - aero_r) * abs(y / (b/2)) + aero_r;

% Geometric Angle of Attack
alpha_g = @(y) (geo_t - geo_r) * abs(y/ (b/2)) + geo_r;

%% Create A Matrix
%
% This will be the matrix of known values through which the linear model
% will be created. When multiplied by the vector of A_n values it produces
% the b vector, which contains known values for geometric angle of attack
% minus zero lift angle of attack

% Preallocate NxN A Matrix
A = zeros(N,N);

% Fill Matrix With Appropriate Values
for i = 1:N %Iterate values of theta
    for j = 1:N %Iterate attached A_n term, only including odd terms
        A(i,j) = ((4*b)/(a(y(i))*c(y(i)))) * sin((2*j-1)*theta(i)) + (2*j-1)*(sin((2*j-1)*theta(i))/theta(i));
    end
end


%% Create b vector
%
% This vector contains all the known values for geometric angle of attack
% minus zero lift angle of attack. This completes the linear system with N
% unknowns
b = zeros(N,1);

for i = 1:N
    b(i) = alpha_g(y(theta(i))) - alpha_0(y(theta(i)));
end

%% Solve linear system for x vector
%
% The x vector contains all of the unknown values for A_n, and can be
% solved using MATLAB's built in linear algebra tools

x = linsolve(A,b);

%% Calculate Return Values
%
% With all of the unknown values calculated, the return values of e, c_L,
% and c_Di can all be calculated

% Coefficient of Lift
c_L = x(1) * pi * AR;

% Calculating δ
%
% All indexing from 2:N is taken for x, because 2j-1 was accounted for when
% solving for the vector, and it only contains the odd A_n values
delta = 0;
for j = 2:N
    delta = delta + (2*j-1) * (x(j)/x(1))^2;
end

% Span Efficiency Factor
e = 1 / (1 + delta);

% Coefficient of Induced Lift
c_Di = c_L^2 / (pi * e * AR);

end

