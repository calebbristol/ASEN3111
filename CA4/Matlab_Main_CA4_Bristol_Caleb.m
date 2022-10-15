%% ASEN 3111 Computational Assignment 04
%  
% Problem Statement:
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 10/31/21
%      


%% Workspace Cleaning

clc
clear
close all;

%% Problem #1
%
% This problem involves using built in MATLAB functions to reproduce
% similar numbers to Anderson examples from chapter 8-10 that use the
% appendices A, B, and C, or figure 9.9.
%
% Exempted Examples: 8.10,  8.16,  8.17,  8.19, 8.20, 8.22, 8.24

    %% Example 8.8
    fprintf('Example 8.8: \n')
    fprintf('  Anderson Provided: App A, M = 3.5 \n')
    fprintf('    p_0 / p = 76.27 \n')
    fprintf('    T_0 / T = 3.45 \n')
    fprintf('  MATLAB Calculated: M = 3.5 \n')
    [p0op,t0ot,~] = isentropic(3.5);
    fprintf('    p_0 / p = %.2f \n',p0op)
    fprintf('    T_0 / T = %.2f \n',t0ot)
      
    %% Example 8.9
    fprintf('Example 8.9: \n')
    fprintf('  Anderson Provided: App A, M_inf = 0.6 \n')
    fprintf('    p_0 / p = 1.276 \n')
    fprintf('    M_1 = 0.9 \n')
    fprintf('  MATLAB Calculated: M_inf = 0.6 \n')
    [p0op,~,~] = isentropic(0.6);
    M = isentropicFindM(1.4,1.69);
    fprintf('    p_0 / p = %.2f \n',p0op)
    fprintf('    M_1 = %.2f \n',M)
    
    %% Example 8.11
    fprintf('Example 8.11: \n')
    fprintf('  Anderson Provided: App B, M_1 = 2 \n')
    fprintf('    p_2 / p_1 = 4.5 \n')
    fprintf('    T_2 / T_1 = 1.687 \n')
    fprintf('    M_2 = 0.5774 \n')
    fprintf('  MATLAB Calculated: M_1 = 2 \n')
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2);
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    T_2 / T_1 = %.2f \n',t2ot1)
    fprintf('    M_2 = %.2f \n',M2n)
    
    %% Example 8.12
    fprintf('Example 8.12: \n')
    fprintf('  Anderson Provided: App A,B, M_1 = 2 \n')
    fprintf('    p0_1 / p_1 = 7.824 \n')
    fprintf('    p0_2 / p0_1 = 0.7209 \n')
    fprintf('  Anderson Provided: App A,B, M_1 = 4 \n')
    fprintf('    p0_1 / p_1 = 151.8 \n')
    fprintf('    p0_2 / p0_1 = 0.1388 \n')
    fprintf('  MATLAB Calculated: M_1 = 2 \n')
    [p0op,t0ot,rho0orho] = isentropic(2);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)
    fprintf('  MATLAB Calculated: M_1 = 4 \n')
    [p0op,t0ot,rho0orho] = isentropic(4);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(4);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)
    
    %% Example 8.13
    fprintf('Example 8.13: \n')
    fprintf('  Anderson Provided: App A,B, M_inf = 2 \n')
    fprintf('    p0_inf / p_inf = 7.824 \n')
    fprintf('    t0_inf / t_inf = 1.8 \n')
    fprintf('    p0_1 / p0_inf = 0.7209 \n')
    fprintf('  Anderson Provided: App A, M_2 = 0.2 \n')
    fprintf('    p_2 / p0_2 = 1 / 1.028 = 0.973 \n')
    fprintf('    t_2 / t0_2 = 1 / 1.008 = 0.992 \n')
    fprintf('  MATLAB Calculated: M_inf = 2 \n')
    [p0op,t0ot,rho0orho] = isentropic(2);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2);
    fprintf('    p0_inf / p_inf = %.2f \n',p0op)
    fprintf('    t0_inf / t_inf = %.2f \n',t0ot)
    fprintf('    p0_1 / p0_inf = %.2f \n',p02op01)
    fprintf('  MATLAB Calculated: M_2 = 0.2 \n')
    [p0op,t0ot,rho0orho] = isentropic(0.2);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(0.2);
    fprintf('    p_2 / p0_2 = %.2f \n',1/p0op)
    fprintf('    t_2 / t0_2 = %.2f \n',1/t0ot)
    
    %% Example 8.14
    fprintf('Example 8.14: \n')
    fprintf('  Anderson Provided: App A,B, M_inf = 10 \n')
    fprintf('    p0_inf / p_inf = 42440 \n')
    fprintf('    t0_inf / t_inf = 21 \n')
    fprintf('    p0_1 / p0_inf = 0.003045 \n')
    fprintf('  Anderson Provided: App A, M_2 = 0.2 \n')
    fprintf('    p_2 / p0_2 = 1 / 1.028 = 0.973 \n')
    fprintf('    t_2 / t0_2 = 1 / 1.008 = 0.992 \n')
    fprintf('  MATLAB Calculated: M_inf = 10 \n')
    [p0op,t0ot,rho0orho] = isentropic(10);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(10);
    fprintf('    p0_inf / p_inf = %.2f \n',p0op)
    fprintf('    t0_inf / t_inf = %.2f \n',t0ot)
    fprintf('    p0_1 / p0_inf = %.6f \n',p02op01)
    fprintf('  MATLAB Calculated: M_2 = 0.2 \n')
    [p0op,t0ot,rho0orho] = isentropic(0.2);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(0.2);
    fprintf('    p_2 / p0_2 = %.2f \n',1/p0op)
    fprintf('    t_2 / t0_2 = %.2f \n',1/t0ot)
    
    %% Example 8.15
    fprintf('Example 8.15: \n')
    fprintf('  Anderson Provided: App B, M_inf = 2 \n')
    fprintf('    M_1 = 2 \n')
    fprintf('    M_2 = 0.5774 \n')
    fprintf('    rho_2 / rho_1 = 2.667 \n')
    fprintf('    t_2 / t_1 = 1.687 \n')
    fprintf('  MATLAB Calculated: M_inf = 2 \n')
    M1 = shock_calcGetM1(4.5);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(M1);
    fprintf('    M_1 = %.2f \n',M1)
    fprintf('    M_2 = %.2f \n',M2n)
    fprintf('    rho_2 / rho_1 = %.2f \n',rho2orho1)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    
    %% Example 8.18
    fprintf('Example 8.18: \n')
    fprintf('  Anderson Provided: App B, M_1 = 3.5 \n')
    fprintf('    M_2 = 0.4512 \n')
    fprintf('    t_2 / t_1 = 3.315 \n')
    fprintf('  MATLAB Calculated: M_1 = 3.5 \n')
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(3.5);
    fprintf('    M_2 = %.2f \n',M2n)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    
    %% Example 8.21
    fprintf('Example 8.21: \n')
    fprintf('  Anderson Provided: App B, M_1 = 3.53 \n')
    fprintf('    M_2 = 0.45 \n')
    fprintf('  MATLAB Calculated: M_1 = 3.53 \n')
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(3.53);
    fprintf('    M_2 = %.2f \n',M2n)
    
    %% Example 8.23
    fprintf('Example 8.23: \n')
    fprintf('  Anderson Provided: App A,B, M_1 = 8 \n')
    fprintf('    p0_1 / p_1 = 9763 \n')
    fprintf('    p0_2 / p0_1 = 0.008488 \n')
    fprintf('  MATLAB Calculated: M_1 = 8 \n')
    [p0op,t0ot,rho0orho] = isentropic(8);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(8);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    fprintf('    p0_2 / p0_1 = %.6f \n',p02op01)
    
    %% Example 9.2
    fprintf('Example 9.2: \n')
    fprintf('  Anderson Provided: App B, M_1 = 1.6, 2, Fig. 9.9, M_1 = 2, theta = 20 \n')
    fprintf('    Beta = 53.4 \n')
    fprintf('    Mn_2 = 0.6684 \n')
    fprintf('    p_2 / p_1 = 2.82 \n')
    fprintf('    t_2 / t_1 = 1.388 \n')
    fprintf('    p0_2 / p0_1 = 0.8952 \n')
    fprintf('    p0_1 / p_1 = 7.824 \n')
    fprintf('    t0_1 / t_1 = 1.8 \n')
    fprintf('  MATLAB Calculated: M_1 = 1.6, 2, theta = 20 \n')
    [p0op,t0ot,rho0orho] = isentropic(2);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.6);
    [ beta ] = rad2deg(BTMeq(deg2rad(20),2));
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    fprintf('    t0_1 / t_1 = %.2f \n',t0ot)
    
    %% Example 9.3
    fprintf('Example 9.2: \n')
    thetaeq=@(b,M1) atand((2*cot(b).*(M1.^2.*(sin(b)).^2-1)./...
    (M1.^2.*(1.4+cos(2.*b))+2)));
    theta = thetaeq(deg2rad(30),2.4);
    fprintf('    theta = %.2f \n',theta)
    
    %% Example 9.4
    fprintf('Example 9.4: \n')
    [ M1 ] = shock_calcGetM1( 3 );
    fprintf('    Mn_1 = %.2f \n',M1)
    
    %% Example 9.5
    fprintf('Example 9.5: \n')
    theta = thetaeq(deg2rad(40),3);
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.93);
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)
    fprintf('    theta = %.2f \n',theta)
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.9);
    fprintf('    p0_3 / p0_2 = %.2f \n',p02op01)
    
    
    %% Example 9.6
    fprintf('Example 9.6: \n')
    [ beta ] = rad2deg(BTMeq(deg2rad(15),5));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2.05);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    
    %% Example 9.7
    fprintf('Example 9.7: \n')
    [ beta ] = rad2deg(BTMeq(deg2rad(10),3.6));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.464);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.358);
    fprintf('    Mn_3 = %.2f \n',M2n)
    fprintf('    p_3 / p_2 = %.2f \n',p2op1)
    fprintf('    t_3 / t_2 = %.2f \n',t2ot1)
    
    %% Example 9.8
    fprintf('Example 9.8: \n')
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(8);
    fprintf('  a) \n')
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(6.9);
    fprintf('  b) \n')
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    t_2 / t_1 = %.2f \n',t2ot1)
    
    %% Example 9.9
    fprintf('Example 9.9: \n')
    [p0op,t0ot,rho0orho] = isentropic(1.5);
    [ v ] = rad2deg(PMeq_findv( 1.5 ));
    fprintf('    v_1 = %.2f \n',v)
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    fprintf('    t0_1 / t_1 = %.2f \n',t0ot)
    [p0op,t0ot,rho0orho] = isentropic(2);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    fprintf('    t0_2 / t_2 = %.2f \n',t0ot)
    
    %% Example 9.10
    fprintf('Example 9.10: \n')
    [p0op,t0ot,rho0orho] = isentropic(10);
    [ v ] = rad2deg(PMeq_findv( 10 ));
    fprintf('    v_1 = %.2f \n',v)
    [ M ] = PMeq_findM(deg2rad(87.3),5);
    fprintf('    M_2 = %.2f \n',M)
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    [p0op,t0ot,rho0orho] = isentropic(6.4);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    
    %% Example 9.11
    fprintf('Example 9.11: \n')
    [ beta ] = rad2deg(BTMeq(deg2rad(15),10));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(3.42);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    Mn_2 = %.2f \n',M2n)
    fprintf('    p_2 / p_1 = %.2f \n',p2op1)
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)
    [p0op,t0ot,rho0orho] = isentropic(10);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    [p0op,t0ot,rho0orho] = isentropic(5.22);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    
    %% Example 9.12
    fprintf('Example 9.12: \n')
    [p0op,t0ot,rho0orho] = isentropic(3);
    [ v ] = rad2deg(PMeq_findv( 3 ));
    fprintf('    v_1 = %.2f \n',v)
    [ M ] = PMeq_findM(deg2rad(54.76),3);
    fprintf('    M_2 = %.2f \n',M)
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    [p0op,t0ot,rho0orho] = isentropic(3.27);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    [ beta ] = rad2deg(BTMeq(deg2rad(5),3));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.177);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    p_3 / p_1 = %.2f \n',p2op1)
    
    %% Example 9.13
    fprintf('Example 9.13: \n')
    fprintf('  a) \n')
    [ v ] = rad2deg(PMeq_findv( 7 ));
    fprintf('    v_1 = %.2f \n',v)
    [ M ] = PMeq_findM(deg2rad(100.97),7);
    fprintf('    M_2 = %.2f \n',M)
    [p0op,t0ot,rho0orho] = isentropic(7);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    [p0op,t0ot,rho0orho] = isentropic(9.56);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    [ beta ] = rad2deg(BTMeq(deg2rad(10),7));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(1.99);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    p_3 / p_1 = %.2f \n',p2op1)
    fprintf('  b) \n')
    [ v ] = rad2deg(PMeq_findv( 7 ));
    fprintf('    v_1 = %.2f \n',v)
    [ M ] = PMeq_findM(deg2rad(95.97),7);
    fprintf('    M_2 = %.2f \n',M)
    [p0op,t0ot,rho0orho] = isentropic(8.1);
    fprintf('    p0_2 / p_2 = %.2f \n',p0op)
    [p0op,t0ot,rho0orho] = isentropic(7);
    fprintf('    p0_1 / p_1 = %.2f \n',p0op)
    [ beta ] = rad2deg(BTMeq(deg2rad(15),7));
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2.79);
    fprintf('    Beta = %.2f \n',beta)
    fprintf('    p_3 / p_1 = %.2f \n',p2op1)
    
    %% Example 10.1
    fprintf('Example 10.1: \n')
    [Msup, Msub] = AoverAstar(1.4,10.25);
    fprintf('    M_e = %.2f \n',Msup)
    [p0op,t0ot,rho0orho] = isentropic(3.95);
    fprintf('    p_e / p_0 = %.2f \n',1/p0op)
    fprintf('    t_e / t_0 = %.2f \n',1/t0ot)
    
    %% Example 10.2
    fprintf('Example 10.2: \n')
    fprintf('  a) \n')
    [Msup, Msub] = AoverAstar(1.4,2);
    fprintf('    M_e = %.2f \n',Msup)
    [p0op,t0ot,rho0orho] = isentropic(2.2);
    fprintf('    p_e / p_0 = %.2f \n',1/p0op)
    fprintf('    t_e / t_0 = %.2f \n',1/t0ot)
    fprintf('  b) \n')
    [Msup, Msub] = AoverAstar(1.4,2);
    fprintf('    M_e = %.2f \n',Msub)
    [p0op,t0ot,rho0orho] = isentropic(0.3);
    fprintf('    p_e / p_0 = %.2f \n',1/p0op)
    fprintf('    t_e / t_0 = %.2f \n',1/t0ot)
    
    %% Example 10.3
    fprintf('Example 10.3: \n')
    [Msub] = isentropicFindM(1.4,1.028);
    fprintf('    M_e = %.2f \n',Msub)
    [Msup, Msub] = AoverAstar(1.4,1.482);
    fprintf('    M_e = %.2f \n',Msub)
    
    %% Example 10.6
    fprintf('Example 10.6: \n')
    [ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(2);
    fprintf('    p0_2 / p0_1 = %.2f \n',p02op01)


%% Problem #2
%
% This problem involves solving Anderson problem 9.14 through the use of
% built in MATLAB functions

    %% Constants
    M_1 = 3;
    epsilon = 10;
    alpha = 15;
    [p01op1,t01ot1,rho01orho1] = isentropic(M_1);

    %% Region 2 (NW)
    v_1 = rad2deg(PMeq_findv(3));
    v_2 = v_1 + alpha - epsilon;
    M_2 = PMeq_findM(deg2rad(v_2),3.1);
    [p02op2,t02ot2,rho02orho2] = isentropic(M_2);
    
    %% Region 3 (NE)
    v_3 = v_2 + 2*epsilon;
    M_3 = PMeq_findM(deg2rad(v_3),3.1);
    [p03op3,t03ot3,rho03orho3] = isentropic(M_3);
    
    %% Region 4 (SW)
    [beta_4] = rad2deg(BTMeq(deg2rad(alpha + epsilon),3));
    Mn_1 = M_1 * sind(beta_4);
    [ M4n,p4op1,rho4orho1,t4ot1,deltasoR,p04op01] = shock_calc(Mn_1);
    M_4 = M4n / sind(beta_4 - (epsilon + alpha));
    [p04op4,t04ot4,rho04orho4] = isentropic(M_4);
    
    %% Region 5 (SE)
    v_4 = rad2deg(PMeq_findv(M_4));
    v_5 = v_4 + 2*epsilon;
    M_5 = PMeq_findM(deg2rad(v_5),M_4);
    [p05op5,t05ot5,rho05orho5] = isentropic(M_5);
    
    %% Pressure Ratios
    %
    % All processes are isentropic except for the oblique shock into region
    % 4. That pressure ratio was calculated accordingly but the rest is a
    % simple algebra problem where p01 = p02 = p03 and p04 = p05
    p2op1 = p01op1 / p02op2;
    p3op1 = p2op1 * p02op2 / p03op3;
    p5op1 = p4op1 * p04op4 / p05op5;
    
    %% Sectional Lift Coeff
    %
    % This formula is just some simple algebra where the lift and drag
    % forces are pressure times area times cosine or sine of the surface's
    % angle. The length cancels out when the mach version of q_inf is used
    % and the resulting equation is below
    c_l = (1 / (1.4*M_1^2*cosd(epsilon))) * ...
        ((p4op1 - p3op1)*cosd(epsilon+alpha) + (p5op1-p2op1)*cosd(alpha-epsilon));
    
    %% Sectional Drag Coeff
    c_d = (1 / (1.4*M_1^2*cosd(epsilon))) * ...
        ((p4op1 - p3op1)*sind(epsilon+alpha) + (p5op1-p2op1)*sind(alpha-epsilon));
    
    %% Print Results
    fprintf('\nAnderson Problem 9.14: \n')
    fprintf('  Sectional Lift Coefficient: %.3f \n',c_l)
    fprintf('  Sectional Drag Coefficient: %.3f \n',c_d)

%% Problem #3
%
% This problem involves reproducing Anderson plot 9.9 through the use of
% built in MATLAB functions

    %% Produce Beta and Mach Vectors
    
    % As defined in Anderson:
    M = [0.05:0.05:0.5 0.6:0.1:2 2.2:0.2:4 4.5 5 6 8 10 20];
    
    % Redefine with atan2 to ensure quadrant correctness
    thetaeq=@(b,M1) atan2((2*cot(b).*(M1^2.*(sin(b)).^2-1)),((M1.^2.*(1.4+cos(2.*b))+2)));
    
    % Create Approximate minimum beta values
    %
    % They're all about 0 but this just stops errors later
    beta_0 = zeros(length(M),1);
    for i = 1:length(M)
        beta_0(i) = BTMeq(0.01,M(i));
    end
    
    % Define beta vectors from 0 to 90
    beta = zeros(length(M),1000);
    for i = 1:length(M)
        beta(i,:) = linspace(beta_0(i),90,1000);
    end
    
    %% Solve for theta with Mach and Beta
    theta = zeros(length(M),1000);
    for i = 1:length(M)
        for j = 1:1000
            theta(i,j) = rad2deg(thetaeq(deg2rad(beta(i,j)),M(i)));
        end
    end

    %% Plotting    
    figure()
    hold on
    for i = 1:length(M)
        plot(theta(i,:),beta(i,:),'k')
    end
    xlim([0 58])
    ylim([0 90])
    xlabel('\theta [degrees]')
    ylabel('\beta [degrees]')
    title('M - \beta - \theta Plot')
    set(gca,'FontSize',14)
    grid on
    hold off