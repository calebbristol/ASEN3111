function [] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
%% Function: Plot_Airfoil_Flow
%
% Plots streamlines, equipotential lines, and pressure contours for flow
% about a thin symmetric airfoil
%
% Utilizes thin airfoil theory and numeric integration of vortex sheet

%% Define Domain
x_step = c/N; %Delta x
x_min = -c;
x_max = 2*c;
y_min = -1.5;
y_max = 1.5;

%% Create mesh over domain with NxN grid points
[x,y] = meshgrid(linspace(x_min,x_max,N),linspace(y_min,y_max,N));


%% Create function for strenth of vortex sheet
gamma = @(x) 2*alpha*V_inf*sqrt((1-(x/c))./(x/c));


%% Create function to solve polar coordinates at point from reference sheet
%
% This will allow the stream and potential function to be in cartesian
%
% The vortex sheet is on the x-axis, where y = 0

theta = @(x,y,x_j) atan2((y),(x-x_j)); %y_j always = 0
radius = @(x,y,x_j) sqrt((x - x_j).^2+y.^2); %Thin airfoil, y_j = 0


%% Solve for strength and circulation at each individual vortex
%
% x_foil describes the x locations evaluated on the airfoil
%
% x = 0 not included because it divides by zero and creates a singularity
% Instead, center of panels is taken, excluding leading and trailing edge
% Because the center of panels is taken, there are N-1 points

x_foil = linspace((x_step/2),c-(x_step/2),N-1);

gamma_x = gamma(x_foil);
Gamma = gamma_x * x_step;


%% Stream Function

    %% Uniform Flow Stream Function
    %
    % Accounts for angle of attack below negative x axis
    psi_UF = V_inf*cos(alpha)*y-V_inf*sin(alpha)*x;


    %% Vortex Stream Function
    %
    % Sum every single vortex up to combine their stream functions together
    % This utilizes the law of superposition with every vortex on the sheet
    %
    % psi_gamma = (Gamma/2pi) * ln(r), r = sqrt(x^2 + y^2)
    %
    % r is measured from the panel in question
    %
    % the sum of each psi_gamma provides the stream function contribution from
    % the entire vortex sheet

    psi_gamma = zeros(N);

    for i = 1:N
        for j = 1:N
            r = radius(x(i,j),y(i,j),x_foil);
            psi_i = Gamma/(2*pi) .* log(r);
            psi_gamma(i,j) = sum(psi_i);
        end
    end


    %% Add Stream Functions Together
    StreamFunction = psi_UF + psi_gamma;
    
    %% Plotting:
    
        %% Determine color levels for contours
        %
        % Code borrowed from Lifting_Cylinder.m

        levmin_psi = min(min(StreamFunction)); % defines the color levels -> trial and error to find a good representation
        levmax_psi = max(max(StreamFunction));
        levels_psi = linspace(levmin_psi,levmax_psi,50)';


        %% Plot streamfunction at levels
        figure()
        contour(x,y,StreamFunction,levels_psi); hold on

        %% Plot airfoil on same graph
        plot(x_foil,zeros(1,length(x_foil)),'k','LineWidth',3);
        title(['Stream Function of Thin Airfoil (N = ',num2str(N),')'])

        %% Adjust axis and label figure
        axis equal
        axis([x_min x_max y_min y_max])
        ylabel('y [m]')
        xlabel('x [m]')
        hold off;
    
    
%% Velocity Potential

    %% Uniform Flow Velocity Potential
    %
    % Accounts for alpha
    phi_UF = V_inf*cos(alpha)*x + V_inf*sin(alpha)*y;
    
    
    %% Vortex Potential Function
    %
    % Sum every single vortex up to combine their potential functions together
    % This utilizes the law of superposition with every vortex on the sheet
    %
    % phi_gamma = -(Gamma/2pi) * theta, theta = arctan(y/(x-x_j))
    %
    % theta is measured from each panel to point in flow being observed
    %
    % the sum of each phi_gamma provides the velocity potential contribution from
    % the entire vortex sheet

    phi_gamma = zeros(N);

    for i = 1:N
        for j = 1:N
            theta_xy = theta(x(i,j),y(i,j),x_foil);
            phi_i = -Gamma/(2*pi) .* theta_xy;
            phi_gamma(i,j) = sum(phi_i);
        end
    end
    
    
    %% Add Velocity Potentials Together
    PotentialFunction = phi_UF + phi_gamma;
    
    
    %% Plotting:
    
        %% Determine color levels for contours
        %
        % Code borrowed from Lifting_Cylinder.m

        levmin_phi = min(min(PotentialFunction)); % defines the color levels -> trial and error to find a good representation
        levmax_phi = max(max(PotentialFunction));
        levels_phi = linspace(levmin_phi,levmax_phi,50)';


        %% Plot PotentialFunction at levels
        figure()
        contour(x,y,PotentialFunction,levels_phi); hold on

        %% Plot airfoil on same graph
        plot(x_foil,zeros(1,length(x_foil)),'k','LineWidth',3);
        title(['Velocity Potential of Thin Airfoil (N = ',num2str(N),')'])

        %% Adjust axis and label figure
        axis equal
        axis([x_min x_max y_min y_max])
        ylabel('y [m]')
        xlabel('x [m]')
        hold off;
    


%% Pressure Contours

    %% Caluclate Velocity from Stream Function
    u = diff(StreamFunction,1,1) ./ diff(y,1,1);
    v = diff(StreamFunction,1,2) ./ diff(x,1,2);
    
    V = sqrt(u(:,1:N-1).^2 + v(1:N-1,:).^2);
    
    
    %% Calculate pressure with Bernoulli's Equation
    p = p_inf - 0.5 * rho_inf * V.^2;
    
    
     %% Plotting:
    
        %% Determine color levels for contours
        %
        % Code borrowed from Lifting_Cylinder.m

        levmin_p = min(min(p)); % defines the color levels -> trial and error to find a good representation
        levmax_p = max(max(p));
        levels_p = linspace(levmin_p,levmax_p,50)';
        
        %% Create Adjusted x and y matrices for plot
        %
        % These have dimensions N-1 becuase the derivative taken above
        x_adj = x(1:N-1,1:N-1);
        y_adj = y(1:N-1,1:N-1);


        %% Plot pressure contour at levels
        figure()
        contour(x_adj,y_adj,p,levels_p); hold on

        %% Plot airfoil on same graph
        plot(x_foil,zeros(1,length(x_foil)),'k','LineWidth',3);
        title(['Pressure Distribution of Thin Airfoil (N = ',num2str(N),')'])

        %% Adjust axis and label figure
        axis equal
        axis([x_min x_max y_min y_max])
        ylabel('y [m]')
        xlabel('x [m]')
        hold off;




end

