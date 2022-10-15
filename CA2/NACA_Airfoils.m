function [x,y] = NACA_Airfoils(m,p,t,c,N)
%% NACA_Airfoils: Returns N panels locations for NACA Airfoil
%
% Takes inputs of N and airfoil geometry
%
% Utilizes Cosine Map from KC Handout

%% Define Circle 
%
% Centered at (c/2,0) with radius c/2
%
% Maps to N+1 points, where the trailing edge counts twice

theta = linspace(0,2*pi,N+1);

% Because this is a cosine map, only the x values will be used
x_circ = zeros(1,length(theta));
for i = 1:length(theta)
    x_circ(i) = (c/2) + (c/2) * cos(theta(i));
end


%% Define Airfoil Thickness Function
%
% Airfoil Camber requires if statement so it's defined further below
y_t = @(x) (t/0.2) * c * (0.2969*sqrt(x/c) - 0.1260*(x/c) - ...
    0.3516*(x/c)^2 + 0.2843*(x/c)^3 - 0.1036*(x/c)^4);

%% Define x and y functions
%
%

xu = @(x) x - y_t(x)* sin(atan(dycdx(x)));
xl = @(x) x + y_t(x)* sin(atan(dycdx(x)));
yu = @(x) y_c(x) + y_t(x) * cos(atan(dycdx(x)));
yl = @(x) y_c(x) - y_t(x) * cos(atan(dycdx(x)));

%% Iterate Across Circle to Determine x and y
%
% Note: The leading edge is only included for even N
% For odd N the leading edge is the center of a panel
x = zeros(N+1,1);
y = zeros(N+1,1);

if mod(N,2) %odd panels case
    for i = 2:ceil(N/2)
        x(i) = xu(x_circ(i));
        x(N+2-i) = xl(x_circ(i));
        y(i) = yu(x_circ(i));
        y(N+2-i) = yl(x_circ(i));
    end
    
    x(1) = c;
    y(1) = 0;
    x(N+1) = c;
    y(N+1) = 0;

elseif ~mod(N,2) %even panels case
    for i = 2:(N/2)
        x(i) = xu(x_circ(i));
        x(N+2-i) = xl(x_circ(i));
        y(i) = yu(x_circ(i));
        y(N+2-i) = yl(x_circ(i));
    end
    x(N/2+1) = 0; %leading edge
    y(N/2+1) = 0;
    
    x(N+1) = c; %trailing edge
    y(N+1) = 0;
    x(1) = c; %counted twice
    y(1) = 0;
end


%% Define Airfoil Camber Function
%
% Outlined in Lab Document
    function [y_camb] = y_c(x)
        if (x >= 0) && (x <= p*c)
            y_camb = (m * (x/(p^2))) * ((2*p) - (x/c));
        elseif (x >= p*c) && (x <= c)
            y_camb = (m * (c-x)/((1-p)^2)) * (1 + (x/c) - (2*p)); 
        end
    end

%% Define Camber Derivative Function
%
% This was solved for analytically
    function [dy_camb] = dycdx(x)
        if (x >= 0) && (x < p*c)
            dy_camb = (2*m/p) - (2 * ((m*x)/(p^2*c)));
        elseif (x >= p*c) && (x <= c)
            dy_camb = (m/(1-p)^2) * (-(2*x) + (2*p));
        end
    end
end

