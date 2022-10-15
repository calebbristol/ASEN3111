%% ASEN 3111 Computational Assignment 02
%  
% Problem Statement:
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 10/10/21
%      


%% Housekeeping

clc
clear
close all;


%% Part 1: Preliminary Work on Thin Airfoils


    %% Problem 1: 
    
        %% Define Constants
        c = 1.5; %[m]
        alpha_d = 6; %[degrees]
        alpha_r = deg2rad(alpha_d); %[rad]
        V_inf = 30; %[m/s]
        p_inf = 101300; %[Pa]
        rho_inf = 1.225; %[kg/m^3]
        
        %% Define N vector
        N = [50 100 250]; %This will be the number of discrete vortices being tested
        
        %% Run Plot_Airfoil_Flow
        for i = 1:length(N)
            Plot_Airfoil_Flow(c,alpha_r,V_inf,p_inf,rho_inf,N(i));
        end
        
%% Part 2: Thick Airfoils


    %% Problem 2: NACA 0012 Airfoil Lifting Flow Analysis
    
    
        %% Determine Nominal Panel Count for Accuracy
        %
        % Takes high N count as "truth data"
        %
        % Calculates error from "truth" coefficient of pressure
        %
        % Determine minimum N within 5 percent error
        
        [x_t,y_t] = NACA_Airfoils(0,0,0.12,1,1000);
        
        [cl_t,cp_t,xc_t] = Vortex_Panel(x_t,y_t,30,0,0);
        
        N_nom = 5:1000; %Tested N values for nominal accuracy
        N_5 = 0; %N_5 is the minimum N within 5 percent error
        
            %% Iterative for loop
            %
            % Checks each N_nom value for accuracy, assigns minimum
            % acceptable value to N_5, then leaves loop
            
            for i = 1:length(N_nom)
                %% Perform Vortex Panel Method with N Panels
                [x_nom,y_nom] = NACA_Airfoils(0,0,0.12,1,N_nom(i));
                [cl_nom,cp_nom,xc_nom] = Vortex_Panel(x_nom,y_nom,V_inf,0,0);

                %% Compute Error
                Err = (cl_t - cl_nom)/cl_t;

                %% Check Error for Acceptable Value
                if abs(Err) <= 0.05
                    N_5 = N_nom(i); %Assign Appropriate Value to N_5
                    break; %Leave the for loop, further Ns unnecessary
                end
            end
            
            
        %% Calculate Coefficient of Lift and Pressure
        %
        % Utilizes previously calculated acceptable N value N_5
        %
        % Note: Vortex_Panel gives really small numbers for coefficient of
        % lift and coefficient of pressure, I was unable to fix this bug
        % because NACA_Airfoils gives the correct boundaries for a NACA
        % airfoil but the coefficient of lift plots are wrong
        
        alpha_0012 = [-6 0 6 9]; %[deg]
        
        cl_0012 = zeros(1,length(alpha_0012));
        cp_0012 = cell(1,length(alpha_0012));
        xc_0012 = cell(1,length(alpha_0012));
        [x_0012,y_0012] = NACA_Airfoils(0,0,0.12,1,N_5);
        
        for i = 1:length(alpha_0012)
            [cl_0012(i),cp_0012{i},xc_0012{i}] = Vortex_Panel(x_0012,y_0012,V_inf,alpha_0012(i),1);
            hold on
            title(['C_p for NACA 0012 at \alpha = ',num2str(alpha_0012(i)),'^o'])
            hold off
        end
        
        
        %% Plotting
        %
        % Plotting coefficient of pressure already done by Vortex_Panel
        % function, plot sectional coefficient of lift vs alpha
        
        figure()
        plot(alpha_0012,cl_0012,'--b','LineWidth',3); hold on
        title('NACA 0012 Airfoil C_l vs. \alpha')
        xline(0,'k')
        yline(0,'k')
        grid on
        xlabel('\alpha [degrees]')
        ylabel('C_l')
        hold off
    
    
    %% Problem 3: Differently Shaped NACA Airfoils Comparisons

    
        %% Define Airfoil Geometry for Airfoils
        %
        % Note: NACA 0012 already defined above
        
        [x_2412,y_2412] = NACA_Airfoils(0.02,0.4,0.12,1,N_5);
        [x_4412,y_4412] = NACA_Airfoils(0.04,0.4,0.12,1,N_5);
        [x_2424,y_2424] = NACA_Airfoils(0.02,0.4,0.24,1,N_5);
        
        
        %% Find Coefficient of Lift for Different alphas
        %
        % Because the relationship should be linear, we can use the alpha
        % values defined above
        
        alpha_foil = alpha_0012;
        
        % Preallocate Vectors
        cl_2412 = zeros(1,length(alpha_foil));
        cp_2412 = cell(1,length(alpha_foil));
        xc_2412 = cell(1,length(alpha_foil));
        cl_4412 = zeros(1,length(alpha_foil));
        cp_4412 = cell(1,length(alpha_foil));
        xc_4412 = cell(1,length(alpha_foil));
        cl_2424 = zeros(1,length(alpha_foil));
        cp_2424 = cell(1,length(alpha_foil));
        xc_2424 = cell(1,length(alpha_foil));

        
        % Run Vortex Panel for all airfoils
        for i = 1:length(alpha_foil)
            [cl_2412(i),cp_2412{i},xc_2412{i}] = Vortex_Panel(x_2412,y_2412,V_inf,alpha_foil(i),0);
            [cl_4412(i),cp_4412{i},xc_4412{i}] = Vortex_Panel(x_4412,y_4412,V_inf,alpha_foil(i),0);
            [cl_2424(i),cp_2424{i},xc_2424{i}] = Vortex_Panel(x_2424,y_2424,V_inf,alpha_foil(i),0);
        end
        
        
        %% Plotting
        %
        % Note: The magnitude of the coefficient of lift is wrong as
        % outlined in the first calculation of coefficient of lift and
        % pressure, but the lift curves look correct in relation to each
        % other so the comparison should still be valid
        
        figure()
        plot(alpha_foil,cl_0012,'--c','LineWidth',2); hold on
        plot(alpha_foil,cl_2412,'--b','LineWidth',2)
        plot(alpha_foil,cl_4412,'--r','LineWidth',2)
        plot(alpha_foil,cl_2424,'--g','LineWidth',2)
        xline(0,'k')
        yline(0,'k')
        grid on
        xlabel('\alpha [degrees]')
        ylabel('C_l')
        legend('NACA 0012','NACA 2412','NACA 4412','NACA 2424')
        title('C_l vs. \alpha for NACA Airfoils')
        
        hold off;
        
        %% Estimate Lift Slope & Zero Lift Angle of Attack
        %
        % Estimate Slope with Differentials
        %
        % Estimate Zero Lift Angle of Attack by Interpolating with Slope
        
        % Calculating Slope
        dcl_0012 = mean(diff(cl_0012) ./ diff(alpha_foil));
        dcl_2412 = mean(diff(cl_2412) ./ diff(alpha_foil));
        dcl_4412 = mean(diff(cl_4412) ./ diff(alpha_foil));
        dcl_2424 = mean(diff(cl_2424) ./ diff(alpha_foil));
        
        % Knowing that second vector element is coefficient of lift at zero
        % angle of attack:
        alpha_0_0012 = -cl_0012(2) / dcl_0012;
        alpha_0_2412 = -cl_2412(2) / dcl_2412;
        alpha_0_4412 = -cl_4412(2) / dcl_4412;
        alpha_0_2424 = -cl_2424(2) / dcl_2424;
        
        % Print Results to Command Window
        fprintf("For NACA 0012:\n")
        fprintf("Lift Slope = %f\n",dcl_0012)
        fprintf("Zero Lift Angle of Attack = %f degrees\n",alpha_0_0012)
        
        fprintf("For NACA 2412:\n")
        fprintf("Lift Slope = %f\n",dcl_2412)
        fprintf("Zero Lift Angle of Attack = %f degrees\n",alpha_0_2412)
        
        fprintf("For NACA 4412:\n")
        fprintf("Lift Slope = %f\n",dcl_4412)
        fprintf("Zero Lift Angle of Attack = %f degrees\n",alpha_0_4412)
        
        fprintf("For NACA 2424:\n")
        fprintf("Lift Slope = %f\n",dcl_2424)
        fprintf("Zero Lift Angle of Attack = %f degrees\n",alpha_0_2424)
