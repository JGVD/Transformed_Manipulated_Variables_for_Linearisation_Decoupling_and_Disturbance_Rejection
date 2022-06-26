%% *J.G. van Dyk (vandyk@sun.ac.za): Final Year Project 2022: Case Study 1*
% Project Title: _Transformed manipulated variables for linearisation, 
% decoupling, and disturbance rejection._ 
% Project Supervisor: Prof S. M. Bradshaw (smb@sun.ac.za).
% Case Study 1: Mixing tank.

%% *Acknowledgement of work of authors and previous work of J. G. van Dyk*

% The author, J. G. van Dyk, acknowledges the use of his previous work, and
% also acknowledges and appreciates the MATLAB and Simulink model inputs of
% Prof S. M. Bradshaw in said work:
% * Bradshaw, S. M. & van Dyk, J. G., 2020. Chemical Engineering 426 
% Assignment 2: Base, PI controller, and PID controller simulation of a 
% three-tank mixing process, Stellenbosch: s.n.
%% *Clear MATLAB workspace, command window and memory, and reset RNG* 

clc;
clf;
close;
clear;
rng('shuffle');
%% *Parameters and inputs*

% System properties
p.V = 35;     % m^3, volume of mixing tank 
p.valve = 50; % % (percentage), valve position
time = 200;   % min, simulation runtime

% System properties at steady state conditions
p.F1 = 0.14;    % m^3/min, flowrate of stream 1
p.F2 = 6.9;     % m^3/min, flowrate of stream 2
p.T1 = 350;     % K, temperature of stream 1
p.T2 = 340;     % K, temperature of stream 2
p.m = 100;      % kg/s, mass flow rate

p.x_A_stream_1 = 100; % %A, concentration of A in stream A 
p.x_A_stream_2 = 1;   % %A, concentration of A in stream B
%% *Calculations* 

% Material balance over mixing tank
xA_0 = ((p.F2*p.x_A_stream_2) + (p.F1*p.x_A_stream_1))/(p.F1 + p.F2); % %A, 
                                         % initial concentration of A
F = p.F1 + p.F2;                         % m^3/h, outlet flowrate
T0 = ((p.T1+p.T2)/(p.F1+p.F2))*p.m;      % K, initial temperature

xA = ((p.F1 + p.F2)*xA_0)/(p.F1 + p.F2); % %A, concentration of A in mixing
% tank

% Solve non-linear differential equations
[t_nl,y_nl] = ode45(@(t,y) func_dydt_nl(y,p), [0 time],T0); % Solve ODE
T_nl = y_nl(:,1);               % Steady state variable data, CV
cv.X0 = 0.1452;                 % m^3/min, F1 initial
cv.X_final = 0.1559;            % m^3/min, F1 final
cv.deltaX = cv.X_final-cv.X0;
cv.Y0 = 2.974;                  % % A, xA initial
cv.Y_final = 3.163;             % % A, xA final
cv.deltaY = cv.Y_final-cv.Y0;
cv.Y_28 = cv.Y0 + (cv.deltaY*0.28);
cv.Y_63 = cv.Y0 + (cv.deltaY*0.63)

% Steady state variable data, DV
dv.X0 = 1;      % % A, xA stream 2 initial
dv.X_final = 2; % A, xA stream 2 final
dv.deltaX = dv.X_final-dv.X0;
dv.Y0 = 2.974;      % A; xA initial
dv.Y_final = 3.951; % A, xA final
dv.deltaY = dv.Y_final-dv.Y0;
dv.Y_28 = dv.Y0 + (dv.deltaY*0.28);
dv.Y_63 = dv.Y0 + (dv.deltaY*0.63)

% Time data
cv.t_step = 20;                 % mins, time of stepchange in MV
cv.t_28 = 21.969- cv.t_step;
cv.t_63 = 24.850 - cv.t_step; 
dv.t_step = 20;                 % mins, time of stepchange in DV
dv.t_28 = 21.690 - dv.t_step;
dv.t_63 = 24.875 - dv.t_step

% Gains
kP = cv.deltaY/cv.deltaX        % Process gain; MV
kD =  dv.deltaY/dv.deltaX       % Disturbance gain; DV
kV = 0.5;

% Time constants
tauP = 1.5*(cv.t_63-cv.t_28)    % mins, time constant for process; MV
tauD = 1.5*(dv.t_63-dv.t_28)    % mins, time constant for disturbance; DV

% FOPDT values for controller fine-tuning
thetaP = cv.t_63-tauP       % mins, intercept of maximum slope with 
%                           initial value; MV
thetaD= dv.t_63-tauD        % mins, intercept of maximum slope with 
%                           initial value; DV
    
% Signal noise generation
varPI = 5E-5;   % Variance for random number generator for signal noise
seedPI = 3;     % Seed for random number generator for signal noise

%% *Feedback controller tuning parameters: PID controller*

deadtimeSP = thetaP/(thetaP+tauP)   % -, fraction of deadtime for step 
%                                   change in SP
kC_PID = 1/(kP*kV)                  % Controller gain for step change in SP
TI_PID = 0.82*(thetaP+tauP)         % mins, integral time for step 
%                                   change in SP
TD_PID = 0.05*(thetaP+tauP)         % mins, derivative time for step 
%                                   change in SP
alphaD = 0.2;                       % Specified derivative filter value

%% *ODE functions for non-linear and linear models of non-isothermal CSTR*

% Non-linear ODE: CA and T
function dydt = func_dydt_nl(y,p)
    T  = y(1);                  
    
    dTdt = (p.F1*(p.T1-T)+p.F2*(p.T2-T))/p.m; % ODE for time-dependent T 
  
    dydt = dTdt;            % Output func_dydt_nl
end % Function func_dydt_nl ends.