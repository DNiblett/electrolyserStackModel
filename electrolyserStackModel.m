%% Electrolyser Stack Model %%
% Developed by Dr Daniel Niblett, 13/11/2025 for EPSRC Ocean Refuel Project
% A 0D model for an electrolyser stack to predict the transient response due to temperature
% Consists of:
% 1 - Stack temperature T_s
% 2 - Combined electrolyte streams T_e
% 4 - Heat transfer to surroundings T_ambient
% 5 - Drift velocity for two-phase flow
% 6 - Mixture approach to properties
% 7 - Temperature dependant properties
%
%       STACK CONFIGURATION
%               P_in
%                |
%             ___v____
%  T_e_in -> |       | T_e ->
%            |  T_s  | 
%  T_c_in -> |_______| T_c -> (not implemented at the moment)
%                |
%                v
%            T_ambient

% Tasks still to do:
% 1 - implement PEM stack option (Completed 10/12/2025)
% 2 - implement minimum load clipping (Completed 07/12/2025)
% 3 - implement cooling water balance (in progress - for now artificially increase electrolyte flow rate)
% 4 - document and automatic plotting (in progress)

clc
clear
close all hidden

%% ---------------------------------------------------- %%
% Simulation Settings - Only Change This Part For Input 
% ----------------------------------------------------- %%

T_set = 50 +273;                        % Desired operating temperature of stack [K]
Power = 2e+6;                           % Constant Power Applied [W]
N_stacks = 3;                           % Number of stacks
electrolyte.flowRatePerCell = 40;       % [L/min] - lumped electrolyte/cooling water.
electrolyte.inletTemperature = 20 +273; % [K] intial temperature of inlet electrolyte/water.

sim.electrolyserType = 'Alkaline';      % Electrolyser Type (Alkaline or PEM)

sim.feedbackLoop = 'on';         % (on) to control the inlet electrolyte/water temperature with feedback loop
sim.minLoadClipping = 'on';      % (on) to clip the power to the electrolyser by the minimum load
sim.readExternalPower = 'on';    % (on) to read external power data from .csv file
sim.interpolatePower = 'on';     % (on) to interpolate power between timesteps in .csv file
sim.temporalScheme = 'implicit'; % explicit (forward Euler), implicit (backwards Euler)
sim.alphaScheme = 'implicit';    % explicit (forward Euler), implicit (backwards Euler)
sim.timeStep = 1;                % timestep for the simulation (default = 1) - slower = more stable
sim.endTime = 60*60;             % Limit for the simulation time.
sim.iterationLimit = 1e+6;       % Limit on the number of iterations for the simulation
twoPhase.active = 'on';          % To make cell include two-phase effects

% PowerData - Replace name here for the external file.
externalData = csvread('Data_7MW_Turbine.csv');

% Extract the time and power from the relevant columns.
powerData_Time = externalData(:,1); % Convert to s
powerData_Power = externalData(:,3).*1e6; % Convert to W

%% ---------------------------------------------------- %%
% END of Simulation Settings 
% ----------------------------------------------------- %%

% Initialise recording of time taken.
tic

%% ---------------------------------------------------- %%
% Simulation Properties
% ----------------------------------------------------- %%

% Two-Phase Drift Flux Properties
twoPhase.C0 = 1;
twoPhase.driftVelocity = 0.2;

% Physical Properties of Materials and Fluids
water.density = 998;            % [kg/m3]
membrane.density = 1200;        % [kg/m3]
steel.density = 7850;           % [kg/m3]
steel.heatCapacity = 460;       % [J/(kg*K)]
h2.heatCapacity = 143;          % [J/(kg*K)]
o2.heatCapacity = 92;           % [J/(kg*K)]
constants.Faraday = 96485;      % [C/mol]
constants.gasConstant = 8.314;  % [J/K*mol]
constants.pressure = 101325;    % [Pa]

if string(sim.readExternalPower) == 'on'
sim.endTime = 111111;
sim.endTime = sim.endTime;
end


%% ---------------------------------------------------- %%
% Stack Properties
% ----------------------------------------------------- %%

if string(sim.electrolyserType) == 'Alkaline'
% Stack input properties
stack.ratedCapacity = 2.2e+6;       % [W]
stack.ratedPotential = 2.0;         % [V]
stack.ratedCurrentDensity = 2500;   % [A/m2]
stack.minimumLoad = 0.2;            % [% of rated capacity]
stack.cellArea = 1.3831*1.3831;     % [m2]
stack.numberOfCells = round(stack.ratedCapacity./(stack.ratedPotential.*stack.ratedCurrentDensity.*stack.cellArea));

% Cell properties Alkaline
cell.PTL.thickness = 500e-6;        % [m]
cell.PTL.porosity = 0.5;            % [-]
cell.BPP.thickness = 0.002;         % [m]
cell.BPP.porosity = 0.2;            % [-] i.e. 20% of volume is channels
cell.BPP.channelHeight = 0.001;     % [m] - height for fluid flow
cell.membrane.thickness = 500e-6;   % [m]
cell.membrane.porosity = 0.5;       % [m] % 0 for PEM, 0.5 for Alkaline

electrolyte.KOHfraction = 0.3; % [% wt] 

% Electrochemical properties
cathode.j0 = 0.1;      % [A/m2]
cathode.b = 120/2303;  % [V]
anode.j0 = 0.05;       % [A/m2]
anode.b = 45/2303;     % [V]

end

if string(sim.electrolyserType) == 'PEM'
% Approximate Stack Mass = 1020 kg per stack (https://nelhydrogen.com/wp-content/uploads/2024/03/PSM-Series_PD-0600-0141-Rev-D.pdf)
% Stack input properties
stack.ratedCapacity = 1.25e+6;       % [W]
stack.ratedPotential = 1.8;         % [V]
stack.ratedCurrentDensity = 20000;  % [A/m2]
stack.minimumLoad = 0.1;            % [% of rated capacity]
stack.cellArea = 0.5*0.5;     % [m2]
stack.numberOfCells = round(stack.ratedCapacity./(stack.ratedPotential.*stack.ratedCurrentDensity.*stack.cellArea));

% Cell properties Alkaline
cell.PTL.thickness = 500e-6;        % [m]
cell.PTL.porosity = 0.5;            % [-]
cell.BPP.thickness = 0.0012;         % [m]
cell.BPP.porosity = 0.2;            % [-] i.e. 20% of volume is channels
cell.BPP.channelHeight = 0.001;     % [m] - height for fluid flow
cell.membrane.thickness = 50e-6;   % [m]
cell.membrane.porosity = 0;       % [m] % 0 for PEM, 0.5 for Alkaline

electrolyte.KOHfraction = 0; % [% wt] 

cathode.j0 = 0.1;            % [A/m2]
cathode.b = 30/2303;        % [V]
anode.j0 = 0.001;           % [A/m2]
anode.b = 40/2303;          % [V]

end

%% ---------------------------------------------------- %%
% Simulation Functions and Conversions
% ----------------------------------------------------- %%

powerData = [powerData_Time powerData_Power]; % convert to W
powerData(:,1) = powerData(:,1) - min(powerData(:,1));

% Inlet Electrolyte/Water Conditions
electrolyte.flowRate = ((electrolyte.flowRatePerCell.*stack.numberOfCells)/60)./1e3; % [m3/s] per stack

% Physical Properties Functions
rho_electrolyte = @(T,x) (997 - (T-298.15).*0.3) + (1030.*x);
mu_electrolyte = @(T,x) (2.414e-5 .* 10.^(247.8./(T-140))).*(1+3.*x);
k_electrolyte = @(T,C_KOH) (-2.041.*(C_KOH) + -0.0028.*(C_KOH.^2) + 0.005332.*(C_KOH.*T) ...
    + 207.2.*(C_KOH/T) + 0.001043.*(C_KOH.^3)+ -0.0000003.*(C_KOH.^2 * T.^2)).*100; % (doi:10.1016/j.ijhydene.2006.10.062)
Cp_electrolyte = @(T,x) (4180 + 0.5.*(T-298.15))*(1-0.5.*x);
k_PEM = @(T) 9.*exp((-12000/8.314) .*((1/T) - (1/353)));            % [S/m] conductivity
mu_h2 = @(T) 8.411e-5 .* (T./273).^(3/2) .* ((273+97)./(T+97));     % [kg m-2 s-1] % check these
mu_o2 = @(T) 1.919e-5 .* (T./273).^(3/2) .* ((273+139)./(T+139));   % [kg m-2 s-1]
rho_gas = @(M,T,p) p.*M./(8.314.*T); % [kg/m3] 
j0_HER_Ni = @(j0,T) j0.*exp((-27000./8.314).*(1/T - 1/333)); % [A/m2]
j0_OER_Ni = @(j0,T) j0.*exp((-65000./8.314).*(1/T - 1/333)); % [A/m2]
j0_HER_PEM = @(j0,T) j0.*exp((-30000./8.314).*(1/T - 1/353)); % [A/m2]
j0_OER_PEM = @(j0,T) j0.*exp((-45000./8.314).*(1/T - 1/353)); % [A/m2]
E0_rev = @(T) 1.229 - 8.5e-4 .*(T - 298); % from https://www.sciencedirect.com/science/article/pii/S0360319924016471


%% ---------------------------------------------------- %%
%       Stack sizing and composition calculations
% ----------------------------------------------------- %%

% Pre-calculations
electrolyte.density = rho_electrolyte(electrolyte.inletTemperature,electrolyte.KOHfraction);
electrolyte.molarity = (electrolyte.KOHfraction.*electrolyte.density)./56.11;
cell.inletArea = cell.BPP.channelHeight.*sqrt(stack.cellArea); % [m2] - based on channel height and cell area
cell.electrolyteFlow = electrolyte.flowRate./stack.numberOfCells; % divide flow rate by number of cells
cell.BPP.length = sqrt(stack.cellArea);
cell.interfaceArea = 2.*(2.*(cell.BPP.channelHeight + cell.BPP.length).*cell.BPP.length); % for cathode and anode
cell.interfaceArea = cell.interfaceArea.*1.5; % increase for the presence of PTL

% Stack size
cell.solidVolume = (cell.BPP.thickness.*(1-cell.BPP.porosity) + cell.PTL.thickness.*(1-cell.PTL.porosity) ...
          + (cell.membrane.thickness.*(1-cell.membrane.porosity))./2 ).*2.*stack.cellArea;
cell.thickness = (cell.BPP.thickness + cell.PTL.thickness + cell.membrane.thickness./2).*2;
cell.liquidVolume = cell.thickness.*stack.cellArea - cell.solidVolume;
stack.totalVolume = cell.thickness.*stack.numberOfCells.*stack.cellArea;
stack.solidVolume = cell.solidVolume.*stack.numberOfCells;
stack.mass = stack.solidVolume.*steel.density;
stack.length = cell.thickness.*stack.numberOfCells;
stack.height = cell.BPP.length;
stack.externalArea = stack.cellArea.*2 + (stack.length.*stack.height).*4;
stack.internalArea = cell.interfaceArea.*stack.numberOfCells;
stack.BOPutilisation = 0.038.*stack.ratedCapacity; % fraction of rated power used by stack

% Pre-allocate variables
T_e = electrolyte.inletTemperature;
T_c = 273+25;
T_s = 273+25;
T_ext = 273+25;

% update physical properties
x_KOH = electrolyte.KOHfraction;
rho_e = rho_electrolyte(T_e,x_KOH);
mu_e = mu_electrolyte(T_e,x_KOH);
cp_e = Cp_electrolyte(T_e,x_KOH);

% electrochem properties
k_e = k_electrolyte(T_e,electrolyte.molarity);
k_mem = k_e.*cell.membrane.porosity.^(2); % for alkaline porous separator
delta_mem = cell.membrane.thickness;

if string(sim.electrolyserType) == 'Alkaline'
j0a = j0_OER_Ni(anode.j0,T_e);   
j0c = j0_HER_Ni(cathode.j0,T_e);
end

if string(sim.electrolyserType) == 'PEM'
j0a = j0_OER_PEM(anode.j0,T_e);   
j0c = j0_HER_PEM(cathode.j0,T_e);
end

ba = anode.b;
bc = cathode.b;
E0a = 1.23;      % [V] equalibrium potential
E_thermo = 1.48; % [V] thermoneutral potential

%% ---------------------------------------------------- %%
% Simulation Calculating Heat Transfer Properties
% ----------------------------------------------------- %%

Q_e_in = cell.electrolyteFlow;
A_e_in = cell.inletArea;
cp_s = steel.heatCapacity;
dh_e = 4.*cell.inletArea./(2.*(cell.BPP.channelHeight + cell.BPP.length)); % hydraulic diameter
u_e_in = Q_e_in./A_e_in;
Pr_e = cp_e.*mu_e./0.6;
Re_e = rho_e.*u_e_in.*dh_e./mu_e;
Nu = 0.023.* Re_e.^0.8 .* Pr_e.^0.4;
h_e = Nu.*0.6./dh_e;
%h_e = (3.* 0.6)./dh_e; % Nusselt number of 4.5 estimate
A_e = stack.internalArea;
h_ext = 10; %W/m.K approximate for now.
A_ext = stack.externalArea;
m_s = stack.mass;
m_e = (stack.totalVolume-stack.solidVolume).*rho_e;
T_e_in = electrolyte.inletTemperature;
N_cells = stack.numberOfCells;
A_cell = stack.cellArea;

% electrochemistry functions
%f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a).*j - Power;
%df = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + E0a + (ba + bc) + 2.*j.*delta_mem./k_mem);

j = 1000;
t = 0;
potential = ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a;
dt = sim.timeStep;
h2ProductionCumlative = 0;
Q_g_out = 0;
I = 0;
alpha = 0;
endTime = sim.endTime;
P_orig = Power;

% Pre-allocate saving lists
result = zeros(1e5,10);
metrics = zeros(1e5,4);
result(1,:) = [0 T_s T_e potential j T_e_in 0 0 0 0];

%% ---------------------------------------------------- %%
%       Start of simulation iteration time loop
% ----------------------------------------------------- %%

% Start of simulation
for iteration = 1:sim.iterationLimit


% ====== Power Conditions ====== %%
Power = P_orig;

% Linear Interpolation For Current Time Power
if string(sim.readExternalPower) == 'on'

if string(sim.interpolatePower) == 'on'
% Get points between current time
p1 = powerData((powerData(:,1) <= t),:);
p1 = p1(end,:);
p2 = powerData((powerData(:,1) > t),:);
p2 = p2(end,:);
powerGradient = (p2(2) - p1(2)) ./ (p2(1)-p1(1));
p3 = powerGradient.*(t - p1(1)) + p1(2);
Power = p3;
end

if string(sim.interpolatePower) == 'off'
sampledPower = powerData(powerData(:,1) > (t - dt.*10) & powerData(:,1) < (t + dt.*10),:);
Power = abs(sampledPower(1,2));
end

end

% Subtract energy for BOP (i.e. heat exchangers, pumps etc.) 
Power = Power - stack.BOPutilisation;
Power(Power<0) = 0; % power supplied cannot be negative

% If power is less than 0, then turn off electrolyte flow rate and temperature control.
if Power == 0
Q_e_in = 0;
else
Q_e_in = cell.electrolyteFlow;
end

% Check if the power is above the minimum load, if so clip power.
if string(sim.minLoadClipping) == 'on'
powerPerStack = Power./N_stacks;
stackMinimumLoad = stack.minimumLoad.*stack.ratedCapacity;
if (powerPerStack) < (stackMinimumLoad)
    Power = 0;
end
end

% if electrolyser temperature gets too hot then turn off power
if T_s >= (90+278)
Power = 0;
end


% ====== End of Power Conditions ====== %

j_old = j;

% update physical properties
x_KOH = electrolyte.KOHfraction;
rho_e = rho_electrolyte(T_e,x_KOH);
mu_e = mu_electrolyte(T_e,x_KOH);
cp_e = Cp_electrolyte(T_e,x_KOH);

% electrochem properties function of temperature
if string(sim.electrolyserType) == 'Alkaline'
k_e = k_electrolyte(T_e,electrolyte.molarity);
k_mem = k_e.*cell.membrane.porosity.^(2); % for alkaline porous separator
j0a = j0_OER_Ni(anode.j0,T_e);   
j0c = j0_HER_Ni(cathode.j0,T_e);
end

if string(sim.electrolyserType) == 'PEM'
k_mem = k_PEM(T_s); % for PEM
j0a = j0_OER_PEM(anode.j0,T_e);   
j0c = j0_HER_PEM(cathode.j0,T_e);
end

delta_mem = cell.membrane.thickness;
E0 = E0_rev(T_e);

% redefine electrochemistry functions
f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a).*j - Power;
df = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + E0a + (ba + bc) + 2.*j.*delta_mem./k_mem);

% Amend for two-phase flow
if string(twoPhase.active) == 'on'
% bubble coverage
bubbleCoverage = 0.023.*j.^0.3; % Vogt https://www.sciencedirect.com/science/article/pii/S001346860400948X
%bubbleCoverage = alpha;
j0a = j0_OER_Ni(anode.j0,T_e).*(1-bubbleCoverage);   
j0c = j0_HER_Ni(cathode.j0,T_e).*(1-bubbleCoverage);

% update properties based on two-phase flow (mixtures)
mu_e = mu_e.*(1-alpha) + mu_h2(T_e).*alpha;
rho_e = rho_e.*(1-alpha) + rho_gas((0.002+0.032)./2,T_e,constants.pressure).*alpha;
cp_e = cp_e.*(1-alpha) + ((h2.heatCapacity + o2.heatCapacity)./2).*alpha;

% added resistance from bubbles (from membraneless paper) - https://www.sciencedirect.com/science/article/pii/S1385894725042780
V_bubble = 2.*(cell.PTL.thickness./k_e) .*(1./(1-alpha.^3) - 1) .*j; % ohm cm2
if string(sim.electrolyserType) == 'PEM'
bubbleCoverage = 0;
V_bubble = 0;
j0a = j0_OER_PEM(anode.j0,T_e);   
j0c = j0_HER_PEM(cathode.j0,T_e);
end

f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a + V_bubble).*j - Power;
df = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + E0a + (ba + bc) + 2.*j.*delta_mem./k_mem);

end


% update time step
t = t + dt;

%V_cell = ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a;

% solve electrochemical model
j = j_old;
for i = 1:10
f1 = f(j);
df1 = df(j);
jn = j - 1.*(f1./df1);
residual = abs(jn-j);
j = jn;
% figure(2)
% plot(i,residual,'o')
% set(gca,'Yscale','log')
% hold on
end

% Damping due to other balance of plant things
% if iteration >1
% tau_lim = 1;
% j_actual = j_old + (dt./tau_lim).*(j - j_old);
% j = j_actual;
% end

tau_lim = dt;
beta = (dt./tau_lim);
jn = (1-beta).*j_old + beta.*j;
j = jn;

potential = Power./(j*A_cell.*N_cells.*N_stacks);
stackCurrent = j.*A_cell; % current (A)
stackVoltage = potential.*N_cells; % voltage (V)

W_elec = stackCurrent.*(potential-E_thermo).*N_cells; % because the current is transferred across each cell

% Solve Temperature Transport Equations
if string(twoPhase.active)=='on'
if iteration > 2
Q_mix_out = Q_e_out + Q_g_out; % add both electrolyte and gas phases
end
end

% Explicit (Forward Euler)
if string(sim.temporalScheme) == 'explicit'
T_s_n = T_s + (-h_e.*A_e*(T_s - T_e) -h_ext.*A_ext.*(T_s - T_ext) + W_elec).*(dt./(m_s.*cp_s));
T_e_n = T_e + (Q_e_in.*N_cells.*rho_e.*cp_e.*T_e_in - Q_e_in.*N_cells.*rho_e.*cp_e.*T_e + h_e.*A_e*(T_s - T_e) ).*(dt./(m_e.*cp_e));
end

% Implicit (Backwards Euler)
if string(sim.temporalScheme) == 'implicit'
F = @(T2,T1) -T2 + T1 + (-h_e.*A_e*(T2 - T_e) -h_ext.*A_ext.*(T2 - T_ext) + W_elec).*(dt./(m_s.*cp_s));
dF = @(T1) -1 + (-h_e.*A_e - h_ext.*A_ext).*(dt./(m_s.*cp_s));

F_e = @(T2,T1) -T2 + T1 + (Q_e_in.*N_cells.*rho_e.*cp_e.*T_e_in - Q_e_in.*N_cells.*rho_e.*cp_e.*T2 + h_e.*A_e*(T_s - T2) ).*(dt./(m_e.*cp_e));
dF_e = @(T1) -1 + (0 - Q_e_in.*N_cells.*rho_e.*cp_e -h_e.*A_e).*(dt./(m_e.*cp_e));

T_s_old = T_s;
T_e_old = T_e;

for k = 1:10
T_s_n = T_s_old - F(T_s_old,T_s)./dF(T_s_old);
T_e_n = T_e_old - F_e(T_e_old,T_e)./dF_e(T_e_old);
res_T_s = abs(T_s_n - T_s_old);
res_T_e = abs(T_e_n - T_e_old);
T_s_old = T_s_n;
T_e_old = T_e_n;

% figure(5)
% plot(k,res_T_s,'o')
% hold on
% plot(k,res_T_e,'o')
% hold on
% pause(0.001)
% set(gca,'Yscale','log')

end
end

T_s = T_s_n;
T_e = T_e_n;

if isnan(T_s)
    break
end

if string(twoPhase.active) == 'on'
% ____ Two-Phase Flows ____ %

h2_volumetric_flow = ((j.*A_cell)./(2.*constants.Faraday)).*((constants.gasConstant.*T_e)./constants.pressure); % [m3/s]
o2_volumetric_flow = ((j.*A_cell)./(4.*constants.Faraday)).*((constants.gasConstant.*T_e)./constants.pressure); % [m3/s/cell]

% (1/V_cell) * d(alpha_g)/dt = Q_in_g - Q_out_g + Q_generate
% d(alpha_g)/dt = (Q_in_g - Q_out_g)./V_cell + Q_generate


Q_generate = h2_volumetric_flow + o2_volumetric_flow;
Q_e_out = Q_e_in;       % assume electrolyte outlet equal to inlet
Q_mix_out = Q_e_out + Q_g_out; % add both electrolyte and gas phases
U_mix_out = Q_mix_out./A_e_in;
Q_g_out = alpha.*(twoPhase.C0.*Q_mix_out + twoPhase.driftVelocity.*A_e_in); % update the Q_gas at outlet

%alpha_n = alpha + (Q_generate - Q_g_out).*(dt./cell.liquidVolume);

if string(sim.alphaScheme) == 'explicit'
alpha_n = alpha + (Q_generate - Q_g_out).*(dt./cell.liquidVolume);
end

%Implicit
if string(sim.alphaScheme) == 'implicit'
F_alpha = @(alpha_2,alpha_1) - alpha_2 + alpha_1 + (Q_generate - alpha_2.*(twoPhase.C0.*Q_mix_out + twoPhase.driftVelocity.*A_e_in)).*(dt./cell.liquidVolume);
dF_alpha = @(alpha_1) - 1 + (0 - (twoPhase.C0.*Q_mix_out + twoPhase.driftVelocity.*A_e_in)).*(dt./cell.liquidVolume);
alpha_n_old = alpha;
for k = 1:5
alpha_n = alpha_n_old - 1.* (F_alpha(alpha_n_old,alpha)./dF_alpha(alpha_n_old));
res_alpha = abs(alpha_n_old - alpha_n);
alpha_n_old = alpha_n;
 % figure(5)
 % plot(k,res_alpha,'o')
 % hold on
% plot(k,res_T_e,'o')
% hold on
% pause(0.001)
% set(gca,'Yscale','log')
end
end
alpha = alpha_n;

% 
% figure(3)
% plot(iteration,alpha,'o')
% hold on

end


% % Control Feedback Loop to adjust inlet electrolyte temperature to meet set-point stack temperature.
% if Power > 0
% if string(sim.feedbackLoop) == 'on'
% % Control Loop for T_e_in so that the electrolyser temperature is met
% T_e_in = T_e_in + (T_set - T_s).*(0.0001);
% end
% end

if Power > 0 && strcmpi(sim.feedbackLoop,'on')
    Kp   = 0.005;         % [1/s] effective gain
    Tmin = 20+273; Tmax = 80+273;

    e = T_set - T_s;
    T_e_in = T_e_in + Kp*e*dt;              % note the *dt
    T_e_in = min(max(T_e_in,Tmin),Tmax);    % clamp
end




% 
% figure(1)
% subplot(2,2,1)
% plot(result(:,1),result(:,2),'k-','LineWidth',1)
% plot(result(:,1),result(:,2),'k-','LineWidth',1)
% hold on
% plot(result(:,1),result(:,3),'b-','LineWidth',1)
% xlabel('Time (s)')
% ylabel('Temperature (K)')
% set(gca,'fontsize',16,'LineWidth',1)
% yyaxis right

% plot(result(:,1),result(:,2),'k-','LineWidth',1)
% hold on
% plot(result(:,1),result(:,3),'b-','LineWidth',1)
% xlabel('Time (s)')
% ylabel('Temperature (K)')
% set(gca,'fontsize',16,'LineWidth',1)
% yyaxis right
% plot(result(:,1),result(:,4),'r-','LineWidth',1)
% plot(result(:,1),result(:,5)./10000,'r-','LineWidth',1)
% ylabel('Potential (V) / currentDensity (A/cm2)')
% %grid on


% ---- Calculate Production Parameters ---- %

h2ProductionRateMoles = (stackCurrent./(2.*constants.Faraday)) .* N_cells .*N_stacks; %I/nF [mol s-1] % because the reaction happens in series
h2ProductionRateKg = h2ProductionRateMoles.*0.002; % kg s-1
h2ProductionRateKg_hour = h2ProductionRateKg.*3600; % kg h-1
h2ProductionCumlative = h2ProductionCumlative + h2ProductionRateKg.*dt;
h2P = h2ProductionCumlative;


result(iteration+1,:) = [t T_s T_e potential j T_e_in h2P Power stackCurrent stackVoltage];

metrics(iteration,:) = [t h2P h2ProductionRateKg_hour alpha ];


if t >= endTime
    break
end


end

result = result(any(result,2),:);
metrics = metrics(any(metrics,2),:);

disp(['Total H2 Produced: ' num2str(h2ProductionCumlative(end)) ' kg' ])
timeTaken = toc;
disp(['Simulation Time: ' num2str(timeTaken) ' s' ])
%legend('T_s','T_e')

% 2 - Power Supplied to Electrolysers (after minimum load cut off)
% 3 - Current Per Stack
% 4 - Voltage Per Stack
% 5 - Capacity used Per Stack
% 6 - Efficiency Per Stack

capacityUsedPerStack = ((result(:,9).*result(:,10))./stack.ratedCapacity).*100;

%% ---------------------------------------------------- %%
%       Plotting of Data in Subplot Arrangement
% ----------------------------------------------------- %%

%figure('Units','normalized','Position',[0 0 1 1]); % full screen

sgtitle(['Hydrogen Produced = ' num2str(h2ProductionCumlative(end)) ' kg' ','...
         ' N_{stack} =  ' num2str(N_stacks) ' ,'...
         ' P_{stack} = ' num2str(stack.ratedCapacity.*1e-6) ' MW' ',' ...
         ' Electrolyser = ' sim.electrolyserType])
% --- Subplot 1 (Power Supplied) ---
subplot(3,4,1)
plot(result(:,1)./60,result(:,8).*1e-6,'k-','LineWidth',1); hold on
plot(result(:,1)./60,(result(:,8).*1e-6)./N_stacks,'b-','LineWidth',1)
xlabel('Time (min)')
ylabel('Power Supplied (MW)')
set(gca,'fontsize',14,'LineWidth',1)
legend('Total','Per Stack')
grid on

% --- Subplot 2 (Current Per Stack) ---
subplot(3,4,2)
plot(result(:,1)./60,result(:,9),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Current Per Stack (A)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 3 (Voltage Per Stack) ---
subplot(3,4,3)
plot(result(:,1)./60,result(:,10),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Voltage Per Stack (V)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 4 (Efficiency Per Stack) ---
efficiencyPerStack = (1.48./(result(:,10)./N_cells)).*100;
subplot(3,4,4)
plot(result(:,1)./60,efficiencyPerStack,'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Efficiency Per Stack (%)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 5 (Capacity Used Per Stack) ---
subplot(3,4,5)
plot(result(:,1)./60,capacityUsedPerStack,'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Capacity Used (%)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 6 (Temperatures Ts and Te) ---
subplot(3,4,6)
plot(result(:,1)./60,result(:,2)-273,'k-','LineWidth',1); hold on
plot(result(:,1)./60,result(:,3)-273,'b-','LineWidth',1)
xlabel('Time (min)')
ylabel('Temperature (°C)')
set(gca,'fontsize',14,'LineWidth',1)
legend('T_s','T_e')
grid on

% --- Subplot 7 (Potential) ---
subplot(3,4,7)
plot(result(:,1)./60,result(:,4),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Potential (V)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 8 (Current Density) ---
subplot(3,4,8)
plot(result(:,1)./60,result(:,5)./10000,'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Current Density (A cm^{-2})')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 9 (Inlet Temperature) ---
subplot(3,4,9)
plot(result(:,1)./60,result(:,6)-273,'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Inlet Electrolyte Temp (°C)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 10 (H2 Production) ---
subplot(3,4,10)
plot(metrics(:,1)./60,metrics(:,2),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('H_2 Production (kg)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 11 (H2 Rate) ---
subplot(3,4,11)
plot(metrics(:,1)./60,metrics(:,3),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('H_2 Rate (kg h^{-1})')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 12 (Gas Fraction) ---
subplot(3,4,12)
plot(metrics(:,1)./60,metrics(:,4),'k-','LineWidth',1)
xlabel('Time (min)')
ylabel('Gas Fraction')
ylim([0 1])
set(gca,'fontsize',14,'LineWidth',1)
grid on

