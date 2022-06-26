%% *J.G. van Dyk (vandyk@sun.ac.za): Final Year Project 2022: Case Study 2*
% Project Title: _Transformed manipulated variables for linearisation, 
% decoupling, and disturbance rejection._ 
% Project Supervisor: Prof S. M. Bradshaw (smb@sun.ac.za).
% Case Study 2: Non-isothermal CSTR.

%% *Acknowledgement of work of authors and previous work of J. G. van Dyk*

% The author, J. G. van Dyk, acknowledges and is grateful for the adapted 
% use of the work of the following authors:
% * Louw, T. M., 2019. Chemical Engineering 344 Tutorial 6 - Linearization. 
% [Online] Available at: 
% https://1drv.ms/b/s!Ar_fp7YSv427o3xDrs7qnfCgJF81?e=4TRFdM
% * Louw, T. M., 2019. Chemical Engineering 344 Tutorial 9 - Modelling 
% examples w/ solutions. [Online] Available at: 
% https://1drv.ms/b/s!Ar_fp7YSv427pjGRK7EBGmXq9hP1?e=sDipO5
% * Marlin, T. E., 2000. Process Control: Designing Processes and 
% Control Systems for Dynamic Performance. 2nd ed. s.l.:McGraw-Hill.

% The author, J. G. van Dyk, acknowledges the use of his previous work:
% * van Dyk, J. G., 2021. Comparison of fault detection techniques in a 
% simulated reaction process - Progress Report , Stellenbosch: s.n.
%% *Clear MATLAB workspace, command window and memory, and reset RNG*

clc;
clf;
close;
clear;
rng('shuffle')
%% *Parameters and inputs*

% Physical properties
p.k0  = 1e10;               % 1/min, pre-exponential coefficient
p.EA_div_R  = 8330.1;       % K, activation energy
p.R   = 1.987;              % cal/mol.K, Universal gas constant
p.Hrx = -130e6;             % cal/kmole, heat of reaction (exothermic)
p.rho = 1e6;                % g/m3, liquid density
p.cpA  = 1;                 % cal/g.K, liquid heat capacity

% System properties
p.CA0 = 2;                  % kmole/m3, inlet concentration
p.T0  = 330;                % K, inlet temperature
p.Tcin = 365;               % K, coolant inlet temperature
p.Fcin = 15;                % m3/min, coolant flowrate
p.F   = 1;                  % m3/min, inlet flowrate
p.V   = 1;                  % m3, reactor volume
p.UA0 = 6e5;                % cal/m3.K, flowrate-dependent heat transfer 
                            % coefficient
p.tau = (p.CA0*p.V)/p.F;    % 1/min, reactor space time
global time;                % min, simulation observation time
time = 10;                  % min, simulation observation time 

% Sensor noise
p.NoiseNLVariance = 0.02;   % Define variance for sensor noise
%% *Calculations*

% Solve non-linear differential equations
[t_nl,y_nl] = ode45(@(t,y) func_dydt_nl(y,p), [0 time], [p.CA0 p.T0]); % ..
                                                % Solve ODE
CA_nl = y_nl(:,1); 
T_nl  = y_nl(:,2);
X_MB_plot = (p.tau*p.k0*exp(-p.EA_div_R./(T_nl)))./...
    (1+p.tau*p.k0*exp(-p.EA_div_R./(T_nl)));    % Conversion i.t.o mole 
                                                % balance on reactant A

% Introduce noise to ODE solutions
p.NoiseNL = p.NoiseNLVariance*randn(size(y_nl)); % Pseudorandom noise
CA_nl = CA_nl + p.NoiseNL(:,1);                  % Add noise to ODE soln.
T_nL = T_nl + (p.NoiseNL(:,2)*200);              % Add noise to ODE soln.

% Determine sensor sampling frequencies, non-linear model
samplingFrequencyT_nL = time/size(T_nl,1)*60;    % Hz, temperature sampling 
                                                 % frequency
samplingFrequencyCA_nL = time/size(CA_nl,1)*60;  % Hz, concentration 
%                                                sampling frequency

% Steady state values
p.Ts = T_nl(size(T_nl,1),1);    % K, set dynamic steady state temperature
p.CAs = CA_nl(size(CA_nl,1),1); % kmole/m3, set dynamic concentration

% Linearisation
p.L_C1 = -(p.F/p.V+p.k0*exp(-p.EA_div_R/p.Ts)*p.V);
p.L_T1 = (p.k0*(p.EA_div_R/p.Ts^2)*exp(-p.EA_div_R/p.Ts)*p.CAs);
p.L_C2 = -(p.UA0*(p.Ts-p.Tcin)*p.Fcin/(p.rho*p.cpA*p.V) - ...
    (p.Hrx*p.k0*exp(-p.EA_div_R/p.Ts)/(p.rho*p.V)));
p.L_T2 = -(p.F/p.V + p.UA0*p.Fcin/(p.rho*p.cpA*p.V) - ...
    (p.Hrx*p.k0*p.EA_div_R/p.Ts^2*exp(-p.EA_div_R/p.Ts)*p.CAs)/...
    (p.rho*p.cpA));

% Solve linear differential equations
[t_dv,y_dv] = ode45(@(t,y) func_dydt_dv(y,p), [0 time], ...
    [p.CA0-p.CAs p.T0-p.Ts]);
CA_l = y_dv(:,1) + p.CAs; 
T_l  = y_dv(:,2) + p.Ts;

% Determine sensor sampling frequencies, linearised model
samplingFrequencyT_L = time/size(T_l,1)*60;   % Hz, temperature sampling 
                                              % frequency (linear)
samplingFrequencyCA_L = time/size(CA_l,1)*60; % Hz, concentration sampling 
                                              % frequency (linear)

% Plot results
subplot(2,2,1); % Plot outlet concentration, CA, as a function of time.
plot(t_nl,CA_nl, '.', t_dv, CA_l,'LineWidth',1.5);
xlabel('Time (min)'); ylabel('C_A (mol/m^3)');
legend('C_A (non-linear)', 'C_A (linearised)', ...
       'Location', 'best');

subplot(2,2,2); % Plot reactor temperature, T, as a function of time.
plot(t_nl,T_nl,'.', t_dv, T_l,'LineWidth',2);
xlabel('Time (min)'); ylabel('Temperature (K)');
legend('T (non-linear)', 'T (linearised)', ...
       'Location', 'best');

subplot(2,2,3); % Plot conversion, X, as a function of time.
plot(t_nl,X_MB_plot, 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('Conversion (-)');

subplot(2,2,4); % Plot conversion, X, as a function of temperature.
plot(T_nl,X_MB_plot, 'LineWidth',1.5);
xlabel('Temperature (K)'); ylabel('Conversion (-)');

%% *System parameters: PID controller*

% Steady state variable data, CV
cv.X0 = 15;         % m^3/min, Fcin at initial steady state
cv.X_final = 16.5;  % m^3/min, Fcin at final steady state after SP change
cv.deltaX = cv.X_final-cv.X0;
cv.Y0 = 381.4;      % K, T at initial steady state
cv.Y_final = 379.6; % K, T at final steady state after SP change
cv.deltaY = cv.Y_final-cv.Y0;
cv.Y_28 = cv.Y0 + (cv.deltaY*0.28);
cv.Y_63 = cv.Y0 + (cv.deltaY*0.63)

% Steady state variable data, DV
dv.X0 = 330;        % K, T0 at initial steady state
dv.X_final = 297;   % K, T0 at final steady state
dv.deltaX = dv.X_final-dv.X0;
dv.Y0 = 381.4;      % K, T at initial steady state
dv.Y_final = 376.8; % K, T at final steady state after DV at t = 20
dv.deltaY = dv.Y_final-dv.Y0;
dv.Y_28 = dv.Y0 + (dv.deltaY*0.28);
dv.Y_63 = dv.Y0 + (dv.deltaY*0.63)

% Time data
cv.t_step = 20;                 % mins, time of stepchange in MV
cv.t_28 = 20.036 - cv.t_step; 
cv.t_63 = 20.081 - cv.t_step; 
dv.t_step = 20                  % mins, time of stepchange in DV
dv.t_28 = 20.041 - dv.t_step;
dv.t_63 = 20.089 - dv.t_step

% Gains
kP = cv.deltaY/cv.deltaX;       % Process gain; MV
kD =  dv.deltaY/dv.deltaX;      % Disturbance gain; DV
kV = 0.5;                       % Valve gain

% Time constants
tauP = 1.5*(cv.t_63-cv.t_28);   % mins, time constant for process; MV
tauD = 1.5*(dv.t_63-dv.t_28);   % mins, time constant for disturbance; DV

% FOPDT values for controller fine-tuning
thetaP = cv.t_63-tauP           % mins, intercept of maximum slope with 
%                               initial value; MV
thetaD= dv.t_63-tauD;           % mins, intercept of maximum slope with 
%                               initial value; DV
    
% Signal noise generation
varPI = 5E-5;   % Variance for random number generator for signal noise
seedPI = 3;     % Seed for random number generator for signal noise

%% *Feedback controller tuning parameters: PID controller*

deadtimeSP = thetaP/(thetaP+tauP)   % -, fraction of deadtime for step 
%                                   change in SP
kC_PID = 1/(kP*kV*100)              % Controller gain for step change in SP
TI_PID = 0.82*(thetaP+tauP)         % mins, integral time for step change 
%                                   in SP
TD_PID = 0.05*(thetaP+tauP)         % mins, derivative time for step change
%                                   in SP
alphaD = 0.2;                       % Derivative filter value
%% *ODE functions for non-linear and linear models of non-isothermal CSTR*

% Non-linear ODE: CA and T
function dydt = func_dydt_nl(y,p)
    global X_MB
    CA = y(1);
    T  = y(2);
    k = p.k0*exp(-p.EA_div_R/(T));       % Reaction rate at temperature T
    UA = p.UA0*p.Fcin;                   % Heat transfer term                      
    
    dCAdt = p.F*(p.CA0 - CA)/p.V - k*CA; % ODE for time-dependent 
                                         % concentration of A
    dTdt  = p.F*(p.T0 - T)/p.V - (p.Hrx*k*CA)/(p.rho*p.cpA*p.V) - ...
        UA*(T-p.Tcin)/(p.rho*p.cpA);    % ODE for time-dependent reactor T
    X_MB = (p.tau*k)/(1 + p.tau*k);     % Conversion at time t
    dydt = [dCAdt; dTdt];               % Output func_dydt_nl
end % Function func_dydt_nl ends.

% Linear ODE: CA and T
function dydtX = func_dydt_dv(y,p)
    CAd = y(1);
    Td  = y(2);
    
    dCAddt = p.L_C1*CAd + p.L_T1*Td; % Linearised ODE for time-dependent 
                                     % concentration of A
    dTddt  = p.L_C2*CAd + p.L_T2*Td; % Linearised ODE for time-dependent 
                                     % reactor T
    
    dydtX = [dCAddt; dTddt];         % Output func_dydt_dv
end % Function func_dydt_dv ends.