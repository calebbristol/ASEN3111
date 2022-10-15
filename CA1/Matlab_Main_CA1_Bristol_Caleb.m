%% ASEN 3111 Computational Assignment 01
%  
% Problem Statement:
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 09/13/21
%      


%% Housekeeping

clc
clear
close all;


%% Problem #1
% 
% The following parts outline an ideal flow over a rotating cylinder


%% Part 1: Analytic Lift and Drag Coefficients
%
% Analytic derivation of sectional lift and drag coefficients
%
% Done using MATLAB built in analytic integration functions

    % Define theta syms variable for MATLAB
    syms theta

    % Get a function for C_p in terms of theta
    C_p_func = CofP(theta);

    % Run integrations for sectional lift and drag coefficients
    c_l_a = double(-0.5 * int(C_p_func * sin(theta),[0 2*pi]));

    c_d_a = double(-0.5 * int(C_p_func * cos(theta),[0 2*pi]));

    % Print values to the command window
    fprintf("Analytic Value for Sectional Coefficient of Lift: %.2f\n",c_l_a);
    fprintf("Analytic Value for Sectional Coefficient of Drag: %.2f\n",c_d_a);


%% Part 2: Composite Trapezoidal Rule
%
% Numeric integration of sectional lift and drag coefficients using the
% composite trapezoidal rule
%
% Done with the help of built in MATLAB function trapz.m which performs a
% composite trapezoidal numeric integration with 

    %% Organize Variables
    
    % Define variable N_cyl
    % N_cyl = vector with N defining number of panels to be used in 
    % integration for cylinder
    N_cyl = [3:2:50];

    % Create cell matrix to store sectional lift and drag coefficients for
    % various values of N, as well as respective theta values
    TrapMat = ones(length(N_cyl),2);
    % Orgainzed as follows
    % Column 1, TrapMat{i,1}: sectional lift coefficient
    % Column 2, TrapMat{i,2}: sectional drag coefficient


    %% Perform Integration 
    
    for i = 1:length(N_cyl)
        %% Create linear vector for theta based off N panels

        theta = linspace(0,2*pi,N_cyl(i));

        %% Calculate Coefficient of Pressure with Respective Theta Values

        c_p = CofP(theta);

        %% Calculate Integrand for Coefficient of Lift and Drag   

        % Lift
        iGrand_l = c_p .* sin(theta);    
        % Drag
        iGrand_d = c_p .* cos(theta);

        %% Perform Numerical Integration with Trapezoidal Rule

        % Coefficient of Lift
        c_l = -0.5 *  traprule(theta,iGrand_l);

        % Coefficient of Drag
        c_d = -0.5 * traprule(theta,iGrand_d);

        %% Write Outputs to TrapMat Matrix for Plotting

        % Column 1, TrapMat(i,1): sectional lift coefficient
        TrapMat(i,1) = c_l;

        % Column 2, TrapMat(i,2): sectional drag coefficient
        TrapMat(i,2) = c_d;

    end

    %% Plotting
    %
    % I commented this section out and added it to the Simpson's Rule
    % section so that the plots could be combined
    
    % Figure 1: Coefficient of Lift vs. N (Trapezoidal Rule)
%     figure(1)
%     
%     plot(N_cyl,TrapMat(:,1)); hold on
%     xlabel("Number of Panels (N)")
%     ylabel("Sectional Coefficient of Lift (c_l)")
%     title("c_l vs. N (Trapezoidal)")
%     xlim([N_cyl(1) N_cyl(end)])
%     set(gca,'fontsize',12)
%     hold off
    
    % Figure 2: Coefficient of Drag vs. N (Trapezoidal Rule)   
%     figure(2)
%     
%     plot(N_cyl,TrapMat(:,2)); hold on
%     xlabel("Number of Panels (N)")
%     ylabel("Sectional Coefficient of Drag (c_d)")
%     title("c_d vs. N (Trapezoidal)")
%     xlim([N_cyl(1) N_cyl(end)])
%     set(gca,'fontsize',12) 
%     hold off
    
    
%% Part 3: Simpson's Rule
%
% Numeric integration of sectional lift and drag coefficients using 
% Simpson's rule
%
% Done with custom function simprule.m which performs a numeric integration
% using Simpson's Rule


    %% Organize Variables
    
    % Create cell matrix to store sectional lift and drag coefficients for
    % various values of N, as well as respective theta values
    SimpMat = ones(length(N_cyl),2);
    % Orgainzed as follows
    % Column 1, TrapMat{i,1}: sectional lift coefficient
    % Column 2, TrapMat{i,2}: sectional drag coefficient


    %% Perform Integration 
    
    for i = 1:length(N_cyl)
        %% Create linear vector for theta based off N panels

        theta = linspace(0,2*pi,N_cyl(i));

        %% Calculate Coefficient of Pressure with Respective Theta Values

        c_p = CofP(theta);

        %% Calculate Integrand for Coefficient of Lift and Drag   

        % Lift
        iGrand_l = c_p .* sin(theta);    
        % Drag
        iGrand_d = c_p .* cos(theta);

        %% Perform Numerical Integration with Simpson's Rule

        % Coefficient of Lift
        c_l = -0.5 * simprule(theta,iGrand_l);

        % Coefficient of Drag
        c_d = -0.5 * simprule(theta,iGrand_d);

        %% Write Outputs to SimpMat Matrix for Plotting

        % Column 1, SimpMat(i,1): sectional lift coefficient
        SimpMat(i,1) = c_l;

        % Column 2, SimpMat(i,2): sectional drag coefficient
        SimpMat(i,2) = c_d;

    end

    %% Plotting
    
    % Figure 1: Coefficient of Lift vs. N 
    figure(1)
    
    %Simpson's Rule
    plot(N_cyl,SimpMat(:,1)); hold on
    
    %Trapezoidal Rule
    plot(N_cyl,TrapMat(:,1))
    
    xlabel("Number of Panels (N)")
    ylabel("Sectional Coefficient of Lift (c_l)")
    title("c_l vs. N")
    xlim([N_cyl(1) N_cyl(end)])
    set(gca,'fontsize',12)
    legend("Simpson's","Trapezoidal")
    hold off
    
    % Figure 2: Coefficient of Drag vs. N   
    figure(2)
    
    %Simpson's Rule
    plot(N_cyl,SimpMat(:,2)); hold on
    
    %Trapezoidal Rule
    plot(N_cyl,TrapMat(:,2))
    
    xlabel("Number of Panels (N)")
    ylabel("Sectional Coefficient of Drag (c_d)")
    title("c_d vs. N")
    xlim([N_cyl(1) N_cyl(end)])
    set(gca,'fontsize',12) 
    legend("Simpson's","Trapezoidal")
    hold off

%% Part 4: 0.2% Relative Error (Trapezoidal Rule)
%
% Considering the analytic models to be "true", find the 0.2% relative
% error when using the trapezoidal rule
%
% Iterates across values for N while calculating the error until
% error reaches the minimum acceptable value

    %% Temporary Variables
    
    Error = 0.002; %[~= 0.2%]
    c_l_error = 1;
    c_l_index = -1;
    c_d_error = 1;
    c_d_index = -1;
    error_index = -1;
    
    
    %% Iterate Across Trapezoidal Matrix
    
    for i = 1:length(N_cyl)
        c_l_error = abs((c_l_a - TrapMat(i,1)) / c_l_a);
        if c_l_error < Error
            c_l_index = i;
            break;
        end
    end
    
    for i = 1:length(N_cyl)
        c_d_error = abs((c_d_a - TrapMat(i,2)) / c_d_a);
        if c_d_error < Error
            c_d_index = i;
            break;
        end
    end
    
    
    fprintf("Smallest N for < 2/10 percent error on c_l (Trapezoidal): %.d\n",N_cyl(c_l_index));

%% Part 5: 0.2% Relative Error (Simpson's Rule)
%
% Considering the analytic models to be "true", find the 0.2% relative
% error when using Simpson's Rule
%
% Iterates across values for N while calculating the error until
% error reaches the minimum acceptable value

    %% Temporary Variables
    
    Error = 0.002; %[~= 0.2%]
    c_l_error = 1;
    c_l_index = -1;
    c_d_error = 1;
    c_d_index = -1;
    error_index = -1;
    
    
    %% Iterate Across Trapezoidal Matrix
    
    for i = 1:length(N_cyl)
        c_l_error = abs((c_l_a - SimpMat(i,1)) / c_l_a);
        if c_l_error < Error
            c_l_index = i;
            break;
        end
    end
    
    for i = 1:length(N_cyl)
        c_d_error = abs((c_d_a - SimpMat(i,2)) / c_d_a);
        if c_d_error < Error
            c_d_index = i;
            break;
        end
    end
    
    
    fprintf("Smallest N for < 2/10 percent error on c_l (Simpson's): %.d\n",N_cyl(c_l_index));


%% Problem 1 Reflection

% The trapezoidal rule actually reached the minimum required error in a
% shorter time than Simpson's Rule. This is likely due to the symmetrical
% nature of the cylinder; further evidenced by the minimum error being
% reached at a value of N = 5 for the trapezoidal rule. This would allow an
% estimate from every 90 degree angle in the circle (counting 0 and 2pi as)
% and would lead to a deceptively low error. Meanwhile, the parabolic
% approximation of Simpson's rule wouldn't make a shockingly low estimate
% quite as fast as the trapezoidal rule for counting 90 degree angles, but
% would likely have higher accuracy overall



%% Problem 2: Airflow Over NACA 0012 Airfoil

% The following parts outline the numeric integration of pressure over a
% NACA 0012 airfoil to determine coefficient of lift and total lift

%% Read in Cp Data

load("Cp.mat");


%% Part 0: Run Numeric Integration
%
% Utilize compound trapezoidal rule to numerically integrate coefficient of
% pressure on upper and lower surfaces, then convert this value into lift
%
% Conversion to lift done using provided values for freestream velocity,
% air density, and pressure


    %% Define Variables
    
    % Constants provided in the lab report:
    V_inf = 40; %[m/s]
    rho_inf = 1.225; %[kg/m^3]
    p_inf = 101.3 * 10^3; %[Pa]
    alpha = 9; %[degrees]
    
    % For use in integration:
    % 
    % Note: N = 1000 appended to vector as "truth" value
    %
    % Also: Reflection section mentions different number of port holes for
    % error. I'm iterating N by 2 so the code runs faster but the
    % reflection includes the correct value
    
    N_foil = [2:2:64, 1000]; %Big N used on both sides, little n (total points) will be 2*N 
    
    
    % For storing integrated c_p values
    cl_foil = ones(length(N_foil),1);
    cd_foil = ones(length(N_foil),1);
    
    
    %% Define Airfoil Function
    
    % Variables given as constants
    t = 12/100; %[Unitless] from NACA 0012 airfoil
    c = 1.5; %[m]
    
    syms x
    % Hard code function for airfoil
    y_t = (t/0.2) * c * (0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c)^2 + 0.2843*(x/c)^3 - 0.1036*(x/c)^4);
    % Create function for derivative of y_t function for use in
    % calculations (will be negative in lower section)
    dy_dx(x) = diff(y_t);
    % dy_dx(x) = (t/0.2) * c * (0.2969/(2*sqrt(x/c)) - 0.1260 - 0.7032*(x/c) + 3*0.2843*(x/c)^2 - 4*0.1036*(x/c)^3);
    
    
    %% Run Numerical Integration 
    
    for i = 1:length(N_foil)
        %% Create x/c values
        
        % These values are equispaced, and range from "0" to 1
        %
        % Note: Starts at 0.01 instead of 0 to avoid divide by 0 in dy/dx
        xc = linspace(0.01,1,N_foil(i));
        
        %% Create c_p values
        
        % Create empty vectors for optimization
        %
        % Has 2 Rows
        % 1) Coefficient of pressure
        % 2) Coefficient of pressure times dy/dx (negative for lower)
        C_p_up = zeros(2,N_foil(i));
        C_p_low = zeros(2,N_foil(i));
        
        % Iterate across airfoil for upper and lower pressures
        for j = 1:N_foil(i)
            C_p_up(1,j) = fnval(Cp_upper,xc(j));
            C_p_up(2,j) = fnval(Cp_upper,xc(j)) * dy_dx(xc(j) * c);

            C_p_low(1,j) = fnval(Cp_lower,xc(j));
            C_p_low(2,j) = fnval(Cp_lower,xc(j)) * -dy_dx(xc(j) * c);
            
        
        end
        
        % Run Trapezoidal Integration on c_p values
        %
        % Utilizes built in numerical integration function trapz
        %
        % The integration of coefficient of pressure over the airfoil is
        % coefficient of the normal force (Anderson eq. 1.15)
        c_n_tot = (trapz(xc,C_p_low(1,:)) - trapz(xc,C_p_up(1,:))) / c;
        
        % The integration of coefficient of pressure multiplied by the
        % derivative of y with respect to x provides the coefficient of the
        % axial force (Anderson eq. 1.16)
        c_a_tot = (trapz(xc,C_p_up(2,:)) - trapz(xc,C_p_low(2,:))) / c;
        
        % Convert normal / axial force coefficient to lift coefficient
        cl_foil(i) = c_n_tot * cosd(alpha) - c_a_tot * sind(alpha);
        
        % Convert normal / axial force into drag coefficient
        cd_foil(i) = c_n_tot * sind(alpha) + c_a_tot * cosd(alpha);
        
    end
    
    
    %% Convert C_l into Lift / Drag
    
    %Multiply by q*c to convert from coefficient
    L_foil = cl_foil * (0.5 * rho_inf * V_inf^2 * c);
    D_foil = cd_foil * (0.5 * rho_inf * V_inf^2 * c);
    
    %Print Lift / Drag per unit span
    fprintf("The Lift Per Unit Span is %.2f [N/m]\n",L_foil)
    fprintf("The Drag Per Unit Span is %.2f [N/m]\n",D_foil)
    

%% Part 1: 5.0% Relative Error

    %% Define "Truth Value"
    
    L_truth = L_foil(end);
    
    %Placeholder variable for 
    N_5 = -1;
    
    for i = 1:length(L_foil)
        % Calculate Error for Each Element
        error = abs(L_foil(i) - L_truth) / L_truth;
        % If error is within acceptable range, write it out and break loop
        if error <= 0.05
            % Multiply by 2 for top and bottom points
            N_5 = N_foil(i) * 2;
            break;
        end
    end
    
    
    fprintf("Smallest n for < 5 percent error on Lift: %.d\n",N_5);

%% Part 2: 1.0% Relative Error
%
% Repeat Process from Part 1 but with different error

    N_1 = -1;
    
    for i = 1:length(L_foil)
        error = abs(L_foil(i) - L_truth) / L_truth;
        if error <= 0.01
            N_1 = N_foil(i) * 2;
            break;
        end
    end
    
    fprintf("Smallest n for < 1 percent error on Lift: %.d\n",N_1);

%% Part 3: 0.2% Relative Error
%
% Same process as previous two parts

    N_dot2 = -1;
    
    for i = 1:length(L_foil)
        error = abs(L_foil(i) - L_truth) / L_truth;
        if error <= 0.002
            N_dot2 = N_foil(i) * 2;
            break;
        end
    end
    
    fprintf("Smallest n for < 0.2 percent error on Lift: %.d\n",N_dot2);

%% Problem 2 Reflection

% In order to recieve a 0.2% error it requires 122 seperate points (See
% note in Problem 2: Part 0: Define Variables) on the airfoil. This is a 
% lot more than could feasibly be placed on a test airfoil, not to mention 
% the difficulty of getting ports near the leading edge or the tail edge
% of the airfoil. Though the required 26 points for 5 percent error
% is still a lot, it is marginally less than that for 0.2 percent
% error. If this experimental data was to be used, it is likely
% that 5 percent error for lift would be reasonable so long as a large
% enough factor of safety was used in designing other components. To
% maximize accuracy based on port location, the leading edge and the tail
% edge could likely be left out. On top of being incredibly difficult, I
% believe the data was skewed by the two points, as the derivative for the
% y_t function divided by x and dy/dx became infinity at the leading edge.
% This is expected because the airfoil goes vertical as it switches from
% top to bottom though it still caused some difficulty in the calculations,
% causing the leading edge to ultimately be exempted anyways.

