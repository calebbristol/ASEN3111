%% ASEN 3111 Computational Assignment 01
%  
% Problem Statement:
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 10/31/21
%      


%% Housekeeping

clc
clear
close all;


%% Problem #1
%
% This problem involves numerically solving the fundamental equation for
% Prandtl's Lifting Line Theory to produce the span efficiency factor, the
% coefficient of lift, and the coefficient of induced drag
%
% Work for problem 1 is done in the file PLLT.m


%% Problem #2
%
% This problem involves an airfoil with a span of 80ft and a taper from c_r
% = 12ft to c_t = 4ft. The root airfoil is a NACA 2412 and the tip is a
% NACA 0012, with a linear spanwise lift slope and zero lift angle of
% attack. Geometric twist causes a linear variation in geometric angle of
% attack from 6 degrees at the root to 1 degree at the tip.
%
% The goal of this problem is to determine the number of odd terms required
% in the series expansion to achieve a 5, 2, and 0.2 percent error at 130
% mph, sea level.

    %% Define Terms of Problem
    b_2 = 80; %[ft]
    c_r_2 = 12; %[ft]
    c_t_2 = 4; %[ft]
    a_r_2 = 6.814225; %[/rad]
    a_t_2 = 6.836323; %[/rad]
    alpha0_r_2 = deg2rad(-2.168806); %[rad]
    alpha0_t_2 = deg2rad(-0.001472); %[rad]
    geo_r_2 = deg2rad(6); %[rad]
    geo_t_2 = deg2rad(1); %[rad]

    %% Run PLLT.m for "truth" value at high N
    [e_2,c_L_2,c_Di_2] = PLLT(b_2,a_t_2,a_r_2,c_t_2,c_r_2,alpha0_t_2,alpha0_r_2,geo_t_2,geo_r_2,1000);

    %% Convert Coefficients into Lift/Drag Values
    rho_2 = 2.3769e-3; %[slug/ft^3]
    V_inf_2 = 130 * (5280/3600); %[ft/s]
    q_inf_2 = 0.5 * rho_2 * V_inf_2^2;
    S_2 = (b_2/2) * (c_t_2 + c_r_2);

    L_2 = c_L_2 * q_inf_2 * S_2;
    Di_2 = c_Di_2 * q_inf_2 * S_2;

    %% Iterate across values of N to find errors
    %
    % Error values calculated using coefficients of lift and induced drag
    % becuase conversion to full lift and induced drag would be multiplication
    % by the same scalar of q_inf*S and wouldn't change the error value

    N_L_5 = -1;
    N_Di_5 = -1;
    N_L_2 = -1;
    N_Di_2 = -1;
    N_L_dot2 = -1;
    N_Di_dot2 = -1;

    
    for i = 2:999
        [e_i,c_L_i,c_Di_i] = PLLT(b_2,a_t_2,a_r_2,c_t_2,c_r_2,alpha0_t_2,alpha0_r_2,geo_t_2,geo_r_2,i);
        error_L = 100 * abs((c_L_2 - c_L_i) / c_L_2);
        error_Di = 100 * abs((c_Di_2 - c_Di_i) / c_Di_2);
        if N_L_5 == -1 && error_L <= 5
            N_L_5 = i;
        end
        if N_Di_5 == -1 && error_Di <= 5
            N_Di_5 = i;
        end
        if N_L_2 == -1 && error_L <= 2
            N_L_2 = i;
        end
        if N_Di_2 == -1 && error_Di <= 2
            N_Di_2 = i;
        end        
        if N_L_dot2 == -1 && error_L <= 0.2
            N_L_dot2 = i;
        end 
        if N_Di_dot2 == -1 && error_Di <= 0.2
            N_Di_dot2 = i;
        end     
        if N_L_dot2 ~= -1 && N_Di_dot2 ~= -1
            break;
        end

    end
    
    %% Print Values to Command Window
    fprintf("The N values required for c_L with 5 percent error: ")
    fprintf(append(num2str(N_L_5),'\n'));
    
    fprintf("The N values required for c_L with 2 percent error: ")
    fprintf(append(num2str(N_L_2),'\n'));
    
    fprintf("The N values required for c_L with 0.2 percent error: ")
    fprintf(append(num2str(N_L_dot2),'\n'));
    
    fprintf("The N values required for c_Di with 5 percent error: ")
    fprintf(append(num2str(N_Di_5),'\n'));

    fprintf("The N values required for c_Di with 2 percent error: ")
    fprintf(append(num2str(N_Di_2),'\n'));
    
    fprintf("The N values required for c_Di with 0.2 percent error: ")
    fprintf(append(num2str(N_Di_dot2),'\n'));

%% Problem #3
%
% Problem 3 involves plotting span efficiency factor against taper ratio
% for a thin wing with no geometric or aerodynamic twist. It will be done
% at AR's of 4, 6, 8, and 10, with at least twenty odd terms in the series
% expansion.


    %% Define terms of problem
    AR = [4 6 8 10];
    c_r_3 = 1;
    a_r_3 = 2*pi; %[/rad]
    a_t_3 = 2*pi; %[/rad]
    alpha0_r_3 = 0; %[rad]
    alpha0_t_3 = 0; %[rad]
    geo_r_3 = 3; %[rad]
    geo_t_3 = 3; %[rad]

    %% Define functions to help with wing geometry
    tap_rat = linspace(0,1,100);
    c_t_func = @(rat) c_r_3 * rat;
    S_func = @(rat) (b_func) * c_r_3 * (1 + rat) / 2;
    b_func = @(AR,rat) AR * c_r_3 * (1 + rat) / 2;
    % Define Wing tip vector
    c_t_3 = c_t_func(tap_rat);


    %% Solve for Span Efficiency Factors
    
    % AR = 4
    b_3_4 = b_func(AR(1),tap_rat);
    e_3_4 = zeros(100,1);
    c_L_3_4 = zeros(100,1);
    c_Di_3_4 = zeros(100,1);
    for i = 1:100
        [e_3_4(i),~,~] = PLLT(b_3_4(i),a_t_3,a_r_3,c_t_3(i),c_r_3,alpha0_t_3,alpha0_r_3,geo_t_3,geo_r_3,20);
    end

    % AR = 6
    b_3_6 = b_func(AR(2),tap_rat);
    e_3_6 = zeros(100,1);
    c_L_3_6 = zeros(100,1);
    c_Di_3_6 = zeros(100,1);
    for i = 1:100
        [e_3_6(i),~,~] = PLLT(b_3_6(i),a_t_3,a_r_3,c_t_3(i),c_r_3,alpha0_t_3,alpha0_r_3,geo_t_3,geo_r_3,20);
    end

    % AR = 8
    b_3_8 = b_func(AR(3),tap_rat);
    e_3_8 = zeros(100,1);
    c_L_3_8 = zeros(100,1);
    c_Di_3_8 = zeros(100,1);
    for i = 1:100
        [e_3_8(i),~,~] = PLLT(b_3_8(i),a_t_3,a_r_3,c_t_3(i),c_r_3,alpha0_t_3,alpha0_r_3,geo_t_3,geo_r_3,20);
    end

    % AR = 10
    b_3_10 = b_func(AR(4),tap_rat);
    e_3_10 = zeros(100,1);
    c_L_3_10 = zeros(100,1);
    c_Di_3_10 = zeros(100,1);
    for i = 1:100
        [e_3_10(i),~,~] = PLLT(b_3_10(i),a_t_3,a_r_3,c_t_3(i),c_r_3,alpha0_t_3,alpha0_r_3,geo_t_3,geo_r_3,20);
    end


    %% Plotting

    figure()
    plot(tap_rat,e_3_4); hold on
    plot(tap_rat,e_3_6)
    plot(tap_rat,e_3_8)
    plot(tap_rat,e_3_10)
    xlabel('Taper Ratio (c_t/c_r)')
    ylabel('Span Efficiency Factor e')
    title('Span Efficiency Factor as Function of Taper Ratio')
    legend('AR = 4','AR = 6','AR = 8','AR = 10')
    grid on
    hold off
    
    
    %% Cross Reference Anderson Figure 5.20
    %
    % Anderson's figure uses induced drag factor instead of span efficiency
    % factor so it requires a conversion
    delta_3_4 = 1./e_3_4 - 1;
    delta_3_6 = 1./e_3_6 - 1;
    delta_3_8 = 1./e_3_8 - 1;
    delta_3_10 = 1./e_3_10 - 1;
    
    %Plot
    figure()
    plot(tap_rat,delta_3_4); hold on
    plot(tap_rat,delta_3_6)
    plot(tap_rat,delta_3_8)
    plot(tap_rat,delta_3_10)
    xlabel('Taper Ratio (c_t/c_r)')
    ylabel('Induced Drag Factor \delta')
    title('\delta as Function of Taper Ratio')
    legend('AR = 4','AR = 6','AR = 8','AR = 10')
    grid on
    hold off
