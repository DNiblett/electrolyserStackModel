%% Dynamic Electrolyser Stack Model
% Developed by Dr Daniel Niblett, 13/11/2025 for EPSRC Ocean Refuel Project
% Updated by Dr Daniel Niblett, 14/07/2026 for EPSRC Ocean Refuel and ION-H2
%
% A 0D lumped thermal-electrochemical model for an electrolyser stack.
% The model resolves:
% 1 - Stack solid temperature, T_s
% 2 - Stack electrolyte outlet temperature, T_e
% 3 - Reservoir electrolyte temperature, T_r
% 4 - Stack inlet electrolyte temperature after pipe losses, T_e_in
% 5 - Heat transfer to ambient air
% 6 - Drift-flux two-phase gas hold-up
% 7 - Temperature-dependent electrolyte and electrochemical properties
%
%     STACK CONFIGURATION

%    Q_W    P_HE             P_S
%     |      |                   |
%   __v__    v                ___v____
%  |     |-> HE -> T_e_in -> |       | T_e ->
%  | T_r |                   |  T_s  | 
%  |_____|<- T_e             |_______| 
%                                |
%                                v
%                            T_ambient

% Notes:
% - PEM, alkaline, and AEM parameter sets are included.
% - The reservoir heater is treated as an immersed heat source.
% - The stack inlet temperature includes an effective pipe heat-loss model.
% - External power and experimental data are optional local inputs.
% - Isothermal condition fixes temperatures to set values.

clc;
clear;
close all hidden;

%% ---------------------------------------------------- %%
% Simulation Settings - Only Change This Part For Input 
% ----------------------------------------------------- %%
PowerLoop = (500:200:20000).*1.87.*(63.62./10000);

for timeLoop = 1

clearvars -except timeLoop PowerLoop timeLoopSave;

% ======== START - Standard Settings to Change for Renewable Energy Integration ========= %

sim.electrolyserType = 'PEM';      % Electrolyser type: 'Alkaline', 'PEM', or 'AEM'
N_stacks = 55*14;                      % Number of stacks
Power = 960*1e6;                       % Operating power if constant current
T_set = 55 + 273;                  % Desired operating temperature of stack [K]
sim.feedbackLoop = 'on';           % (on) to control the inlet electrolyte/water temperature with feedback loop
sim.minLoadClipping = 'on';        % (on) to clip the power to the electrolyser by the minimum load
sim.readExternalPower = 'on';     % (on) to read external power data from .csv file
sim.timeStep = 2;                  % timestep for the simulation (default = 1) - slower = more stable
sim.endTime = 60*60*24*1;           % Limit for the simulation time.
sim.writeTime = 10;                % Resolution of the simulation timescale
sim.isothermal = 'off';            % If isothermal, simulation ignores temperature dynamics
degradation.active = 'on';         % Turn on if degradation model active.
externalPowerFile = 'Wind_power_2020_AR1_1s.mat'; % External file with power data (time|power [W])
plotting.figureScale = 'days';    % What scale to plot figures. 'min' 's' 'hours' 'days'

if strcmpi(sim.readExternalPower,'on')
                % to min   to hour to day
    %sim.endTime = 60 * 60 * 24 * (365/12);
    sim.endTime = 60 * 60 * 24 * 30;
     sim.endTime = 60 * 60 * 24*365;
end

% ======== END - Standard Settings to Change for Renewable Energy Integration ========= %


%Power = PowerLoop(timeLoop);                           % Constant Power Applied [W]
%electrolyte.flowRatePerCell = 4.2/22;       % [L/min] - lumped electrolyte/cooling water.

electrolyte.flowRatePerCell = 20;       % [L/min] - lumped electrolyte/cooling water.
electrolyte.inletTemperature = 55 + 273; % [K] initial temperature of inlet electrolyte/water.

sim.interpolatePower = 'on';     % (on) to interpolate power between timesteps in .csv file
sim.temporalScheme = 'implicit'; % explicit (forward Euler), implicit (backwards Euler)
sim.alphaScheme = 'implicit';    % explicit (forward Euler), implicit (backwards Euler)
sim.iterationLimit = 100e+6;       % Limit on the number of iterations for the simulation
sim.type = 'large';               % large adds directly the amount of water that is consumed
sim.imposeFluctuating = 'off'     % During constant current operation, how to impose fluctuating conditions
twoPhase.active = 'on';          % To make cell include two-phase effects
usePorousElectrodeCorrelation = 'no'; % To account for electrolyte conduction in catalyst layer (AEM only)


% External power data. This is only loaded if sim.readExternalPower = 'on'.
powerData_Time  = [];
powerData_Power = [];

% externalPowerFile = 'Data_7MW_Turbine.csv';
% 
% if strcmpi(sim.readExternalPower,'on')
%     if exist(externalPowerFile,'file') ~= 2
%         error('External power file not found: %s', externalPowerFile);
%     end
% 
%     externalData = readmatrix(externalPowerFile);
% 
%     % Extract the time and power from the relevant columns.
%     powerData_Time  = externalData(:,1);       % [s]
%     %powerData_Power = externalData(:,3).*1e6;  % [W]
%     powerData_Power = (externalData(:,3)./max(externalData(:,3))).*2300;  % [W]
% 
% end


if strcmpi(sim.readExternalPower,'on')
    if exist(externalPowerFile,'file') ~= 2
        error('External power file not found: %s', externalPowerFile);
    end

    externalData = load(externalPowerFile);

    % Extract the time and power from the relevant columns.
    powerData_Time  = double(externalData.time_s);       % [s]
    %powerData_Power = externalData(:,3).*1e6;  % [W]
    %powerData_Power = double((externalData.power_MW./max(externalData.power_MW)).*2300);  % [W]
    powerData_Power = double((externalData.power_MW)).*1e6.*(960/15);  % [W]

end



if strcmpi(sim.readExternalPower,'off')
  powerData_Time = linspace(0,sim.endTime,sim.endTime)';
  powerData_Power = ones(size(powerData_Time,1),1).*Power;

    powerData(:,1) = powerData_Time;
    powerData(:,2) = powerData_Power;  
end

  




% Balance of Plant Settings

sim.waterToReservoir = 'off';       % Toggle on/off
tank.waterFillRate = 0.42;          % [L/min]
tank.waterFillStartTime = 120.*60;  % [s] time when cooling-water injection is first allowed to start. Change this to your measured start time.
tank.waterFillFrequency = 38.*60;   % [s] how often the water is turned on after tank.waterFillStartTime
tank.waterFillDuration = 0.63.*60;  % [s] for how long the water feeds the reservoir each cycle

% Cooling Fan Power

sim.inletCoolingFan = 'on';         % Change the cooling 

% Pipe heat loss between electrolyte reservoir and stack inlet
UA_pipe = 0.6;                      % [W/K] effective pipe heat loss, tune from experiment

% Heater Settings

startUpHeaterTime = 60*60*2;        % startup heater duration [s]
clipPowerDuringStartup = 'off';
alwaysClipPower = 'off';            % To clip power always.
sim.startUpHeater = 'off';
startUpHeaterRatedPower = 135;      % [W] measured electrical heater power, V*I
startUpHeaterMaxPower = startUpHeaterRatedPower; % [W] retained for compatibility
startUpHeaterGain = Inf;           % [W/K] Inf gives on/off rated-power heater control



%% ---------------------------------------------------- %%
%                 Simulation Properties
% ----------------------------------------------------- %%

% Initialise recording of time taken.
tic

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
tank.density = 1800;            % [kg/m3]



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

% Degradation Properties

        % Mature commercial alkaline stack:
        % DOE 2022 status = 3.2 mV per 1000 h = 3.2 microV/h
        degradation.steadyRate       = 3.2e-6/3600; % [V/s]
        degradation.currentExponent  = 1.5;         % [-]

        % Dynamic alkaline degradation data are sparse.
        % Use conservative sensitivity assumptions.
        degradation.lowLoadRate      = 2.0e-6/3600; % [V/s], extra at near-zero load
        degradation.rampLoss         = 0.25e-6;     % [V] per full-rated-current change
        degradation.startLoss        = 1.0e-6;      % [V/start]
        degradation.stopLoss         = 1.0e-6;      % [V/stop]

        degradation.lowLoadLimit     = 0.20;         % fraction rated current
        degradation.onLimit          = 0.01;
        degradation.rampDeadband     = 1e-3;         % ignore <0.1% rated changes

        % High catalyst-loading baseline from review
        degradation.steadyRate      = 5.0e-6/3600;
        degradation.currentExponent = 2.0;

        degradation.lowLoadRate     = 5.0e-6/3600;
        degradation.rampLoss        = 0.5e-6;

        % Calibrated for one complete ON/OFF cycle per hour
        % and 30 microV/h cycling degradation
        degradation.startLoss       = 15e-6;
        degradation.stopLoss        = 15e-6;

        degradation.highCurrentThreshold = 2.0e4;
        degradation.highCurrentReference = 3.0e4;
        degradation.highCurrentExtraRate = ...
            45e-6/3600;
        stack.BOPutilisation = (1/0.985)*0.07*stack.ratedCapacity; % fraction of rated power used by stack 7% for Alkaline and 5% for PEM

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

%- Lumped degradation parameters -

        % DOE-style commercial system degradation assumption:
        % approximately 2.3 microV/h
        degradation.steadyRate       = 2.3e-6/3600; % [V/s]
        degradation.currentExponent  = 1.5;         % [-]

        % Dynamic operation can raise measured PEM voltage degradation.
        degradation.lowLoadRate      = 2.0e-6/3600; % [V/s], additional low-load loss
        degradation.rampLoss         = 1.0e-6;      % [V] per full-rated-current change
        degradation.startLoss        = 2.0e-6;      % [V/start]
        degradation.stopLoss         = 2.0e-6;      % [V/stop]

        degradation.lowLoadLimit     = 0.10;         % PEM tolerates lower loads
        degradation.onLimit          = 0.01;
        degradation.rampDeadband     = 1e-3;


        % High catalyst-loading baseline from review
        degradation.steadyRate      = 5.0e-6/3600;
        degradation.currentExponent = 2.0;

        degradation.lowLoadRate     = 5.0e-6/3600;
        degradation.rampLoss        = 0.5e-6;

        % Calibrated for one complete ON/OFF cycle per hour
        % and 30 microV/h cycling degradation
        degradation.startLoss       = 15e-6;
        degradation.stopLoss        = 15e-6;

        degradation.highCurrentThreshold = 2.0e4;
        degradation.highCurrentReference = 3.0e4;
        degradation.highCurrentExtraRate = ...
            45e-6/3600;
   stack.BOPutilisation = (1/0.985)*0.05*stack.ratedCapacity; % fraction of rated power used by stack 7% for Alkaline and 5% for PEM

end

if string(sim.electrolyserType) == 'AEM'
% Stack input properties
stack.ratedCapacity = 2300;       % [W]
stack.ratedPotential = 1.87;         % [V]
stack.ratedCurrentDensity = 8400;   % [A/m2]
stack.minimumLoad = 0.4;            % [fraction of rated capacity]
stack.cellArea = 63.62./10000;     % [m2]
stack.numberOfCells = round(stack.ratedCapacity./(stack.ratedPotential.*stack.ratedCurrentDensity.*stack.cellArea));

% Cell properties Alkaline
cell.PTL.thickness = 500e-6;        % [m]
cell.PTL.porosity = 0.5;            % [-]
cell.BPP.thickness = 0.002;         % [m]
cell.BPP.porosity = 0.4;            % [-] i.e. 20% of volume is channels
cell.BPP.channelHeight = 0.001;     % [m] - height for fluid flow
cell.membrane.thickness = 50e-6;   % [m]
cell.membrane.porosity = 0;       % [m] % 0 for PEM, 0.5 for Alkaline

electrolyte.KOHfraction = 0.01; % [% wt] 

% Electrochemical properties
cathode.j0 = 5.5;      % [A/m2]
cathode.b = 70/2303;  % [V]
anode.j0 = 5.5;       % [A/m2]
anode.b = 70/2303;     % [V]

%- Lumped degradation parameters -

        % DOE-style commercial system degradation assumption:
        % approximately 2.3 microV/h
        degradation.steadyRate       = 2.3e-6/3600; % [V/s]
        degradation.currentExponent  = 1.5;         % [-]

        % Dynamic operation can raise measured PEM voltage degradation.
        degradation.lowLoadRate      = 2.0e-6/3600; % [V/s], additional low-load loss
        degradation.rampLoss         = 1.0e-6;      % [V] per full-rated-current change
        degradation.startLoss        = 2.0e-6;      % [V/start]
        degradation.stopLoss         = 2.0e-6;      % [V/stop]

        degradation.lowLoadLimit     = 0.10;         % PEM tolerates lower loads
        degradation.onLimit          = 0.01;
        degradation.rampDeadband     = 1e-3;
         stack.BOPutilisation = (1/0.985)*0.05*stack.ratedCapacity; % fraction of rated power used by stack 7% for Alkaline and 5% for PEM

end

%% ---------------------------------------------------- %%
% Simulation Functions and Conversions
% ----------------------------------------------------- %%

if strcmpi(sim.readExternalPower,'on')
    powerData = [powerData_Time powerData_Power];
    powerData(:,1) = powerData(:,1) - min(powerData(:,1));
else
    %powerData = [];
end

% Inlet Electrolyte/Water Conditions
electrolyte.flowRate = ((electrolyte.flowRatePerCell.*stack.numberOfCells)/60)./1e3; % [m3/s] per stack

% Physical Properties Functions
rho_electrolyte = @(T,x) (997 - (T-298.15).*0.3) + (1030.*x);
mu_electrolyte = @(T,x) (2.414e-5 .* 10.^(247.8./(T-140))).*(1+3.*x);
k_electrolyte = @(T,C_KOH) (-2.041.*(C_KOH) + -0.0028.*(C_KOH.^2) + 0.005332.*(C_KOH.*T) ...
    + 207.2.*(C_KOH/T) + 0.001043.*(C_KOH.^3)+ -0.0000003.*(C_KOH.^2 * T.^2)).*100; % (doi:10.1016/j.ijhydene.2006.10.062)
%Cp_electrolyte = @(T,x) (4180 + 0.5.*(T-298.15))*(1-0.5.*x);
Cp_electrolyte = @(T,x) 4180 - 30.*(x.*100) - 0.10.*(x.*100).^2 + 2.*((T-298.15)-25);
k_PEM = @(T) 9.*exp((-12000/8.314) .*((1/T) - (1/353)));            % [S/m] conductivity
k_AEM = @(T) (0.524.*18 - 0.318).*exp((-1270) .*((1/T) - (1/353)));            % [S/m] conductivity
mu_h2 = @(T) 8.411e-5 .* (T./273).^(3/2) .* ((273+97)./(T+97));     % [kg m-2 s-1] % check these
mu_o2 = @(T) 1.919e-5 .* (T./273).^(3/2) .* ((273+139)./(T+139));   % [kg m-2 s-1]
rho_gas = @(M,T,p) p.*M./(8.314.*T); % [kg/m3] 
j0_HER_Ni = @(j0,T) j0.*exp((-27000./8.314).*(1/T - 1/333)); % [A/m2]
j0_OER_Ni = @(j0,T) j0.*exp((-65000./8.314).*(1/T - 1/333)); % [A/m2]
j0_HER_PEM = @(j0,T) j0.*exp((-30000./8.314).*(1/T - 1/353)); % [A/m2]
j0_OER_PEM = @(j0,T) j0.*exp((-45000./8.314).*(1/T - 1/353)); % [A/m2]
j0_HER_AEM = @(j0,T) j0.*exp((-24890./8.314).*(1/T - 1/303)); % [A/m2]
j0_OER_AEM = @(j0,T) j0.*exp((-31000./8.314).*(1/T - 1/303)); % [A/m2]
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
%stack.BOPutilisation = 0.038.*stack.ratedCapacity; % fraction of rated power used by stack
% add mass of bipolar plates and the manifolds? 20/04/26
stack.BPPthickness = 0.05 .* 2;
stack.mass = stack.mass + stack.cellArea.*stack.BPPthickness.*steel.density;

% Electrolyte reservoir details
tank.lx = 0.25;
tank.ly = 0.105;
tank.lz = 0.16;

if string(sim.type) == 'large'
tank.lx = 0.25 .*5;
tank.ly = 0.105 .*5;
tank.lz = 0.16 .*5;
end

tank.volume = tank.lx.*tank.ly.*tank.lz;   % total electrolyte reservoir volume
tank.fillFraction = 0.65; %0.71;     % how much of tank is filled.
tank.shellThickness = 0.005; % thickness of tank shell
tank.shellVolume = tank.lx.*tank.ly.*tank.shellThickness.*2 + tank.lx.*tank.lz.*tank.shellThickness.*2 + ...
tank.ly.*tank.lz.*tank.shellThickness.*2;
tank.shellMass = tank.shellVolume.*tank.density;
tank.totalMass = tank.volume.*tank.fillFraction.*water.density + tank.shellMass;
tank.mass = tank.volume.*tank.fillFraction.*water.density;
tank.liquidVolume = tank.volume.*tank.fillFraction;
tank.shellArea = 2.*tank.lx.*tank.ly + 2.*tank.lx.*tank.lz + 2.*tank.ly.*tank.lz;
A_tank = tank.shellArea;

T_all = 19;
% Pre-allocate variables
T_e = electrolyte.inletTemperature; % electrolyte inlet temperature
T_c = 273+55;   % cell temperature
T_s = 273+55;   % stack temperature
T_ext = 273+17.5; % external temperature
T_r = 273+54;   % reservoir temperature
T_w_in = 273+17.5; % External Water Conditions


% update physical properties
x_KOH = electrolyte.KOHfraction;
rho_e = rho_electrolyte(T_e,x_KOH);
mu_e = mu_electrolyte(T_e,x_KOH);
cp_e = Cp_electrolyte(T_e,x_KOH);

% electrochem properties
k_e = k_electrolyte(T_e,electrolyte.molarity);
k_mem = k_e.*cell.membrane.porosity.^(2); % for alkaline porous separator
delta_mem = cell.membrane.thickness;

% if string(sim.electrolyserType) == 'Alkaline' 
% j0a = j0_OER_Ni(anode.j0,T_e);   
% j0c = j0_HER_Ni(cathode.j0,T_e);
% end
% 
% if string(sim.electrolyserType) == 'PEM'
% j0a = j0_OER_PEM(anode.j0,T_e);   
% j0c = j0_HER_PEM(cathode.j0,T_e);
% end
% 
% if string(sim.electrolyserType) == 'AEM'
% j0a = j0_OER_AEM(anode.j0,T_e);   
% j0c = j0_HER_AEM(cathode.j0,T_e);
% end

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

if string(sim.electrolyserType) == 'AEM'
k_mem = k_AEM(T_s); % for AEM
j0a = j0_OER_AEM(anode.j0,T_e);   
j0c = j0_HER_AEM(cathode.j0,T_e);
end


ba = anode.b;
bc = cathode.b;
E0a = 1.23;      % [V] equilibrium potential
E_thermo = 1.48; % [V] thermoneutral potential

% --- Lumped degradation model (cell-voltage loss) --- %
% Rates are initial literature-informed sensitivity values and should be
% calibrated for the selected electrolyser and operating profile.
% https://www.sciencedirect.com/science/article/pii/S1364032125008433

% degradation.steadyRate = 10e-6/3600;     % [V/s] at rated current (10 uV/h)
% degradation.currentExponent = 1.5;       % current dependence
% degradation.rampLoss = 2e-6;             % [V] per rated-current change
% degradation.startLoss = 2e-6;            % [V/start]
% degradation.stopLoss = 1e-6;             % [V/stop]
% degradation.lowLoadRate = 10e-6/3600;    % [V/s] maximum extra low-load loss
% degradation.lowLoadLimit = 0.2;          % fraction of rated current
% degradation.onLimit = 0.01;              % fraction of rated current
% degradation.rampDeadband = 1e-5;         % minimum |Delta j|/j_rated

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
h_ext = 8; %W/m.K approximate for now.
h_ext_tank = 5;
A_ext = stack.externalArea;
m_s = stack.mass;
m_e = (stack.totalVolume-stack.solidVolume).*rho_e;
T_e_in = electrolyte.inletTemperature;
N_cells = stack.numberOfCells;
A_cell = stack.cellArea;

% electrochemistry functions
%f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a).*j - Power;
%df = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + E0a + (ba + bc) + 2.*j.*delta_mem./k_mem);

j = 10;
j_old = j;
t = 0;
potential = ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a;
dt = sim.timeStep;
h2ProductionCumulative = 0;
Q_g_out = 0;
I = 0;
alpha = 0;
endTime = sim.endTime;
P_orig = Power;
FanPower = 0;
heaterPower = 0;
Q_consumed = 0;
t_write = 0;
writeIter = 0;
power_old = Power;

% Degradation state variables
degradePotential = 0;                    % [V/cell]
degradeSteady = 0;                       % [V/cell]
degradeLowLoad = 0;                      % [V/cell]
degradeRamp = 0;                         % [V/cell]
degradeCycles = 0;                       % [V/cell]
j_degrade_old = 0;
wasRunning = false;
numberOfStarts = 0;
numberOfStops = 0;
totalPower = 0;

V_r = tank.liquidVolume;
mass_r = V_r.*water.density;

% Pre-allocate saving lists
Q_e_out = Q_e_in;
% Results columns:
% 1 time [s], 2 T_s [K], 3 T_e,out [K], 4 cell potential [V],
% 5 current density [A/m2], 6 T_e,in [K], 7 cumulative H2 [kg],
% 8 power [W], 9 stack current [A], 10 stack voltage [V], 11 T_r [K],
% 12 total degradation [V/cell], 13 steady, 14 ramp, 15 low-load,
% 16 start/stop degradation [V/cell]
result = zeros(1e7,17);

% Metrics columns:
% 1 time [s], 2 cumulative H2 [kg], 3 H2 rate [kg/h],
% 4 gas fraction [-], 5 inlet fan/cooling power [W]
metrics = zeros(1e7,6);

result(1,:) = [0 T_s T_e potential j T_e_in 0 0 0 0 T_r 0 0 0 0 0 0];

% Water filling parameters
% waterFilltimer/waterFillCounter are retained for compatibility, but the
% active water-fill logic now uses tank.waterFillStartTime and mod().
waterFilltimer = 0;
waterFillCounter = 0;
waterFill = 'off';

% Start up heater properties
startUpHeaterActive = 'off';


%% ---------------------------------------------------- %%
%       Start of simulation iteration time loop
% ----------------------------------------------------- %%

powerTime = powerData(:,1);
powerValue = powerData(:,2);

powerIndex = 1;
nPower = numel(powerTime);

% Start of simulation
for iteration = 1:sim.iterationLimit


% ====== Power Conditions ====== %%
Power = P_orig;

% Linear Interpolation For Current Time Power
if string(sim.readExternalPower) == 'on'

% if string(sim.interpolatePower) == 'on'
% % Get points between current time
% p1 = powerData((powerData(:,1) <= t),:);
% p1 = p1(end,:);
% p2 = powerData((powerData(:,1) > t),:);
% p2 = p2(end,:);
% powerGradient = (p2(2) - p1(2)) ./ (p2(1)-p1(1));
% p3 = powerGradient.*(t - p1(1)) + p1(2);
% Power = p3;
% end


if strcmpi(sim.interpolatePower,'on')

    % Move index forward only when t passes the next data point
    while powerIndex < nPower - 1 && t >= powerTime(powerIndex + 1)
        powerIndex = powerIndex + 1;
    end

    % Handle times before or after the available dataset
    if t <= powerTime(1)

        Power = powerValue(1);

    elseif t >= powerTime(end)

        Power = powerValue(end);

    else

        t1 = powerTime(powerIndex);
        t2 = powerTime(powerIndex + 1);

        P1 = powerValue(powerIndex);
        P2 = powerValue(powerIndex + 1);

        fraction = (t - t1) / (t2 - t1);

        Power = P1 + fraction * (P2 - P1);
    end
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

% if startUpHeater is on, then turn off power to electrolyser
if string(clipPowerDuringStartup) == 'on'

    if t < startUpHeaterTime
    Power = 0;
    end

    if t >= startUpHeaterTime 
    Power = P_orig;
    end

    if string(alwaysClipPower) == 'on'
    Power = 0;
    end

end


if Power > 0
    if j_old <= 0 
    j_old = 1000;
    j = j_old;
    end
end


% Fluctuating Power Profiles

%Power = %2300 + 2000*sign(sin(2*pi*t/1800)); % switches every 15 min
%Power = P_on*(mod(t,t_on+t_off) < t_on);
%Power = 2300 *(mod(t,(60*10) + (60*30)) < 60*10);

% P_on  = 500;      % W


if string(sim.imposeFluctuating) == 'on'

% ----- V1 - On/Off ------ %

 t_on  = 60*1;      % For how many seconds is it on
 t_off = 60*1;     % For how many seconds is it off
 isOn = (mod(t,t_on+t_off) < t_on);
 Power = P_orig*isOn;

 % ----- V2 - sine fluctuating ----- %

 
%  P_rated = P_orig;
% P_mean  = P_orig;
% 
% Power = P_mean ...
%       + 250*sin(2*pi*t/(45*60)) ...
%       + 120*sin(2*pi*t/(12*60)) ...
%       + 60*sin(2*pi*t/(3*60));
% 
% Power = max(0,min(Power,P_rated));
end
% ----- V3 -------
% 
% P_rated = 2000; % W
% 
% Power = 0.35*P_rated ...
%       + 0.25*P_rated*sin(2*pi*t/(4*3600)) ...
%       + 0.15*P_rated*sin(2*pi*t/(45*60)) ...
%       + 0.08*P_rated*sin(2*pi*t/(8*60));
% 
% Power = max(0,min(Power,P_rated));
% 
% isOn = Power > 0;

% 
% 
 if (Power>0 && power_old == 0)
  j = 1000; % A/m^2
    potential = ba.*log(j./j0a) ...
              + bc.*log(j./j0c) ...
              + j.*delta_mem./k_mem ...
              + E0a;
else
    j = j;
    potential = potential;   % open circuit voltage
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

if string(sim.electrolyserType) == 'AEM'
k_mem = k_AEM(T_s); % for AEM
j0a = j0_OER_AEM(anode.j0,T_e);   
j0c = j0_HER_AEM(cathode.j0,T_e);
end

delta_mem = cell.membrane.thickness;
E0 = E0_rev(T_e);

% redefine electrochemistry functions
f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a + degradePotential).*j - Power;
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

if string(sim.electrolyserType) == 'AEM'
V_bubble = 2.*((cell.PTL.thickness)./k_e) .*(1./(1-alpha.^3) - 1) .*j; % ohm cm2
j0a = j0_OER_AEM(anode.j0,T_e).*(1-bubbleCoverage);   
j0c = j0_HER_AEM(cathode.j0,T_e).*(1-bubbleCoverage);
end


if string(usePorousElectrodeCorrelation) == 'yes'
catalystThickness = 10e-6;
Gamma    = (1e+7.*j0c.* (catalystThickness.^2)) ./ (k_mem .* cathode.b);                % = 1/W
Lambda  = (j .* catalystThickness) ./ (k_mem .* cathode.b);

% Calculate dimensionless parameters
cathodeOverpotential = (2.2333.*(Gamma.^-0.0006).*(Lambda.^0.2347)-1.7933 -0.9988.*log(Gamma) + 0.7876.*log(Lambda)).*cathode.b;
%cathodeOverpotential = (0.0381.*(Gamma.^-0.0013).*(Lambda.^0.709)+0.8966-0.9971.*log(Gamma)+1.2523.*log(Lambda)).*bc;

catalystThickness = 10e-6;
Gamma    = (1e+7.*j0a.* (catalystThickness.^2)) ./ (2.*k_e .* anode.b);                % = 1/W
Lambda  = (j .* catalystThickness) ./ (2.*k_e .* anode.b);
anodeOverpotential = (2.2333.*(Gamma.^-0.0006).*(Lambda.^0.2347)-1.7933 -0.9988.*log(Gamma) + 0.7876.*log(Lambda)).*anode.b;

else

anodeOverpotential = 0;
cathodeOverpotential = 0;



end

V_bubble = V_bubble + anodeOverpotential + cathodeOverpotential;



f = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + j.*delta_mem./k_mem + E0a + V_bubble + degradePotential).*j - Power;
df = @(j) N_stacks.*N_cells.*A_cell.*(ba.*log(j./j0a) + bc.*log(j./j0c) + E0a + (ba + bc) + j.*delta_mem./k_mem);

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

if Power <=0
    potential = 0;
    j = 0;
end

% ---- Lumped degradation update ----
if strcmpi(degradation.active,'on')
    isRunning = j >= degradation.onLimit.*stack.ratedCurrentDensity;
    startEvent = isRunning && ~wasRunning;
    stopEvent = ~isRunning && wasRunning;
    loadFraction = max(j,0)./stack.ratedCurrentDensity;
    rampFraction = abs(j-j_degrade_old)./stack.ratedCurrentDensity;

    rSteady = isRunning.*degradation.steadyRate.*loadFraction.^degradation.currentExponent;
    rLow = isRunning.*(loadFraction < degradation.lowLoadLimit).* ...
        degradation.lowLoadRate.*max(1-loadFraction./degradation.lowLoadLimit,0);
    dVramp = degradation.rampLoss.*max(rampFraction-degradation.rampDeadband,0);
    dVcycle = startEvent.*degradation.startLoss + stopEvent.*degradation.stopLoss;

    degradeSteady = degradeSteady + rSteady.*dt;
    degradeLowLoad = degradeLowLoad + rLow.*dt;
    degradeRamp = degradeRamp + dVramp;
    degradeCycles = degradeCycles + dVcycle;
    degradePotential = degradeSteady + degradeLowLoad + degradeRamp + degradeCycles;

    numberOfStarts = numberOfStarts + startEvent;
    numberOfStops = numberOfStops + stopEvent;
    j_degrade_old = j;
    wasRunning = isRunning;
end

stackCurrent = j.*A_cell; % current (A)
stackVoltage = potential.*N_cells; % voltage (V)

W_elec = stackCurrent.*(potential-E_thermo).*N_cells; % because the current is transferred across each cell

% Solve Temperature Transport Equations
if string(twoPhase.active)=='on'
if iteration > 2
Q_mix_out = Q_e_out + Q_g_out; % add both electrolyte and gas phases
end
end


% Cooling-water injection into the reservoir. Injection is disabled until
% tank.waterFillStartTime. After this, it follows a repeating
% on/off cycle defined by tank.waterFillFrequency and tank.waterFillDuration.
if string(sim.waterToReservoir) == 'on' && t >= tank.waterFillStartTime
    timeSinceWaterStart = t - tank.waterFillStartTime;
    timeInWaterCycle = mod(timeSinceWaterStart,tank.waterFillFrequency);

    if timeInWaterCycle <= tank.waterFillDuration
        waterFill = 'on';
    else
        waterFill = 'off';
    end
else
    waterFill = 'off';
end

% Check if water flow rate from external is on
if string(waterFill) == 'on'
    Q_w_in = ((tank.waterFillRate) /60)/1000;
else
    Q_w_in = 0;
end

% Startup heater control: immersed heater inside the electrolyte reservoir.
% This heat source is added directly to the reservoir energy balance below,
% instead of artificially increasing the inlet-line temperature T_e_in.
heaterPower = 0;
if t < startUpHeaterTime && strcmpi(sim.startUpHeater,'on')
    startUpHeaterActive = 'on';
    if T_r < T_set.*0.95
        heaterPower = startUpHeaterRatedPower;                % [W] measured immersed heater power
    else
        heaterPower = 0;                                      % [W] thermostat off at/above setpoint
    end
else
    startUpHeaterActive = 'off';
end

% % ONLY for validation purposes
% if heaterPower ==0
%     Q_e_in = Q_e_in.*0.01;
% end

% Explicit (Forward Euler)
if string(sim.temporalScheme) == 'explicit'
T_s_n = T_s + (-h_e.*A_e*(T_s - T_e) -h_ext.*A_ext.*(T_s - T_ext) + W_elec).*(dt./(m_s.*cp_s));
T_e_n = T_e + (Q_e_in.*N_cells.*rho_e.*cp_e.*T_e_in - Q_e_in.*N_cells.*rho_e.*cp_e.*T_e + h_e.*A_e*(T_s - T_e) ).*(dt./(m_e.*cp_e));

% Reservoir
T_r_n = T_r + (Q_e_out.*N_cells.*rho_electrolyte(T_e,x_KOH).*Cp_electrolyte(T_e,x_KOH).*T_e - Q_e_in.*N_cells.*rho_electrolyte(T_r,x_KOH).*Cp_electrolyte(T_r,x_KOH).*T_r - h_ext_tank.*A_tank*(T_r - T_ext) + ...
        - Q_w_in.*rho_electrolyte(T_w_in,0).*Cp_electrolyte(T_w_in,0).*T_w_in + heaterPower).*(dt./(mass_r.*Cp_electrolyte(T_r,x_KOH)));

% reservoir volume balance
V_r_n = V_r + (Q_e_out - Q_e_in + Q_w_in).*dt;

V_r = V_r_n;
mass_r = V_r.*electrolyte.density;

end

% Implicit (Backwards Euler)
if string(sim.temporalScheme) == 'implicit'
F = @(T2,T1) -T2 + T1 + (-h_e.*A_e*(T2 - T_e) -h_ext.*A_ext.*(T2 - T_ext) + W_elec).*(dt./(m_s.*cp_s));
dF = @(T1) -1 + (-h_e.*A_e - h_ext.*A_ext).*(dt./(m_s.*cp_s));

F_e = @(T2,T1) -T2 + T1 + (Q_e_in.*N_cells.*rho_e.*cp_e.*T_e_in - Q_e_in.*N_cells.*rho_e.*cp_e.*T2 + h_e.*A_e*(T_s - T2) ).*(dt./(m_e.*cp_e));
dF_e = @(T1) -1 + (0 - Q_e_in.*N_cells.*rho_e.*cp_e -h_e.*A_e).*(dt./(m_e.*cp_e));

% Reservoir implicit balance. Fluid properties are evaluated explicitly at
% the old temperatures, while T_r is solved implicitly in the dominant tank
% outflow and ambient-loss terms. heaterPower acts directly in the tank.
rho_r_old = rho_electrolyte(T_r,x_KOH);
cp_r_old  = Cp_electrolyte(T_r,x_KOH);
rho_e_old = rho_electrolyte(T_e,x_KOH);
cp_e_old  = Cp_electrolyte(T_e,x_KOH);
rho_w_old = rho_electrolyte(T_w_in,0);
cp_w_old  = Cp_electrolyte(T_w_in,0);
F_r = @(T2,T1) -T2 + T1 + (Q_e_out.*N_cells.*rho_e_old.*cp_e_old.*T_e - Q_e_in.*N_cells.*rho_r_old.*cp_r_old.*T2 - h_ext_tank.*A_tank.*(T2 - T_ext) + ...
        - Q_w_in.*rho_w_old.*cp_w_old.*T_w_in + heaterPower).*(dt./(mass_r.*cp_r_old));
dF_r = @(T1) -1 + (- Q_e_in.*N_cells.*rho_r_old.*cp_r_old - h_ext_tank.*A_tank).*(dt./(mass_r.*cp_r_old));

T_s_old = T_s;
T_e_old = T_e;
T_r_old = T_r;

T_s_n = T_s_old;
T_e_n = T_e_old; 
T_r_n = T_r_old; 

if string(sim.isothermal) == 'off'
for k = 1:10
T_s_n = T_s_old - F(T_s_old,T_s)./dF(T_s_old);
T_e_n = T_e_old - F_e(T_e_old,T_e)./dF_e(T_e_old);
T_r_n = T_r_old - F_r(T_r_old,T_r)./dF_r(T_r_old);
res_T_s = abs(T_s_n - T_s_old);
res_T_e = abs(T_e_n - T_e_old);
res_T_r = abs(T_r_n - T_r_old);
T_s_old = T_s_n;
T_e_old = T_e_n;
T_r_old = T_r_n;

% figure(5)
% plot(k,res_T_s,'o')
% hold on
% plot(k,res_T_e,'o')
% hold on
% pause(0.001)
% set(gca,'Yscale','log')

end
end

% reservoir volume balance (for AEM specific stack)
V_r_n = V_r + (Q_e_out - Q_e_in + Q_w_in).*dt;

% for larger scale stacks
if string(sim.type)== 'large'
V_r_n = V_r + (Q_e_out - Q_e_in + Q_w_in + Q_consumed).*dt;
end

V_r = V_r_n;
mass_r = V_r.*electrolyte.density;


end


T_s = T_s_n;
T_e = T_e_n;
T_r = T_r_n;

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

% Stack inlet temperature from reservoir after heat loss through feed pipe.
% Q_e_in is per cell in this model, so the total electrolyte flow to the
% stack is Q_e_in*N_cells. The exponential form prevents unphysical
% over-cooling at low but non-zero flow rates.
Q_e_pipe = max(Q_e_in.*N_cells,eps);
T_e_in_noCooling = T_ext + (T_r - T_ext).*exp(-UA_pipe./(rho_e.*Q_e_pipe.*cp_e));

if Power > 0 && strcmpi(sim.feedbackLoop,'on')
    Kp   = 0.005;         % [1/s] effective gain
    Tmin = 20+273; Tmax = 80+273;

    e = T_set - T_s;
    T_e_in = T_e_in + Kp*e*dt;              % note the *dt
    T_e_in = min(max(T_e_in,Tmin),Tmax);    % clamp
else
    % Without an active inlet-temperature controller, the stack inlet is the
    % reservoir outlet after pipe heat loss to ambient.
    T_e_in = T_e_in_noCooling;
end

if string(sim.inletCoolingFan) == 'on'
  
    if string(startUpHeaterActive) == 'off'
        if T_s > (20 + 273) % limit when it can be active

    Kp   = 5;          % proportional gain
    Tmin = 10 + 273;      % minimum possible cooled inlet [K]
    Tmax = T_e_in_noCooling; % fan can only cool below the pipe-loss inlet temperature

    % positive error means stack is too hot
    e = T_s - T_set;

    % only cool if stack is above setpoint
    coolingSignal = max(e,0);

    % update inlet temperature
    T_e_in = T_e_in_noCooling - Kp.*coolingSignal;

    % clamp: cannot cool below Tmin, cannot heat above feed temperature
    T_e_in = min(max(T_e_in,Tmin),Tmax);

    % calculate required cooling power [W]
    FanPower = Q_e_pipe.*electrolyte.density .* cp_e .* (T_e_in_noCooling - T_e_in);
        end
    end

end

% Startup heater is handled before the temperature equations as an
% immersed heat source in the reservoir energy balance.



%% ---- Consumption of Electrolyte ---- %%

% Consumed Water
Q_consumed = ((stackCurrent.*N_stacks.*0.018)./(2.*constants.Faraday.*water.density)).*dt;

% mass balance of electrolyte
Q_e_out = Q_e_in - Q_consumed;

% Water that is evaporated in the oxygen separator

% ---- Calculate Production Parameters ---- %

h2ProductionRateMoles = (stackCurrent./(2.*constants.Faraday)) .* N_cells .*N_stacks; %I/nF [mol s-1] % because the reaction happens in series
h2ProductionRateKg = h2ProductionRateMoles.*0.002; % kg s-1
h2ProductionRateKg_hour = h2ProductionRateKg.*3600; % kg h-1
h2ProductionCumulative = h2ProductionCumulative + h2ProductionRateKg.*dt;
h2P = h2ProductionCumulative;

totalPower = totalPower + Power.*dt; %[Ws];



%result(iteration+1,:) = [t T_s T_e potential j T_e_in h2P Power stackCurrent stackVoltage T_r];
%metrics(iteration,:) = [t h2P h2ProductionRateKg_hour alpha FanPower heaterPower];

% ---- Save every X time step ----- %
t_write = t_write + dt;

if t_write > sim.writeTime
writeIter = writeIter + 1;
result(writeIter,:) = [t T_s T_e potential j T_e_in h2P Power stackCurrent stackVoltage T_r ...
    degradePotential degradeSteady degradeRamp degradeLowLoad degradeCycles totalPower];
metrics(writeIter,:) = [t h2P h2ProductionRateKg_hour alpha FanPower heaterPower];
t_write = 0;

end





%V_r_save(iteration) = V_r;

if t >= endTime
    break
end

if string(sim.isothermal) == 'on'

    T_s = T_set;
    T_e = T_set;
    T_e_in = T_set;
    T_r = T_set;
    

end

power_old = Power;
end

result = result(any(result,2),:);
metrics = metrics(any(metrics,2),:);
%%

totalTime = (result(end,1)./3600) %h % 
h2perS = h2ProductionCumulative(end)./(result(end,1)); %kg/s
h2perY = h2perS.*3600.*24.*365; %kg/year
CF= (totalPower./1e6 ./3600)/(Power*(totalTime)*1e-6)

disp(['Electrolyser Type: ' sim.electrolyserType ', N_stacks: ' num2str(N_stacks)] )
disp(['Total Power Supplied: ' num2str(totalPower./1e6 ./3600) ' MWh'])
disp(['Capacity Factor: ' num2str(CF) ' (-)'])
disp(['Total H2 Produced: ' num2str(h2ProductionCumulative(end)) ' kg' ])
disp(['Total H2 Produced (kg/hr): ' num2str(h2ProductionCumulative(end)/(totalTime)) ' kg/hr' ])
disp(['Total H2 Production Rate: ' num2str(h2perY./1000) ' t/year' ])
disp(['Total Energy Consumption: ' num2str((totalPower./1e3 ./3600)./h2ProductionCumulative(end) ) ' kWh/kg  '])
timeTaken = toc;
disp(['Simulation Time: ' num2str(timeTaken) ' s' ])
disp(['Degradation: ' num2str(degradePotential.*1e3,'%.4f') ' mV'])
disp(['Starts / stops: ' num2str(numberOfStarts) ' / ' num2str(numberOfStops)])
disp(['Averaged Degradation Rate : ' num2str(degradePotential.*1e6 ./ totalTime ,'%.4f') ' µV/h'])
disp(['Averaged Degradation Rate : ' num2str(degradePotential ./ (totalTime./(24.*365)) ,'%.4f') ' V/year'])


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

plotScale = 1;
plotXlabel = 'Time (s)';
if strcmpi(plotting.figureScale,'days')
plotScale = 3600*24;
plotXlabel = 'Time (days)';
end
if strcmpi(plotting.figureScale,'hours')
plotScale = 3600;
plotXlabel = 'Time (h)';
end
if strcmpi(plotting.figureScale,'min')
plotScale = 60;
plotXlabel = 'Time (min)';
end


%figure('Units','normalized','Position',[0 0 1 1]); % full screen

sgtitle(['Hydrogen Produced = ' num2str(h2ProductionCumulative(end)) ' kg' ','...
         ' N_{stack} =  ' num2str(N_stacks) ' ,'...
         ' P_{stack} = ' num2str(stack.ratedCapacity.*1e-6) ' MW' ',' ...
         ' Electrolyser = ' sim.electrolyserType]);
% --- Subplot 1 (Power Supplied) ---
subplot(3,4,1)
plot(result(:,1)./plotScale,result(:,8).*1e-6,'k-','LineWidth',1); hold on
plot(result(:,1)./plotScale,(result(:,8).*1e-6)./N_stacks,'b-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Power Supplied (MW)')
set(gca,'fontsize',14,'LineWidth',1)
legend('Total','Per Stack')
grid on

% --- Subplot 2 (Current Per Stack) ---
subplot(3,4,2)
plot(result(:,1)./plotScale,result(:,9),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Current Per Stack (A)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 3 (Voltage Per Stack) ---
subplot(3,4,3)
plot(result(:,1)./plotScale,result(:,10),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Voltage Per Stack (V)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 4 (Efficiency Per Stack) ---
efficiencyPerStack = (1.48./(result(:,10)./N_cells)).*100;
subplot(3,4,4)
plot(result(:,1)./plotScale,efficiencyPerStack,'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Efficiency Per Stack (%)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 5 (Capacity Used Per Stack) ---
subplot(3,4,5)
plot(result(:,1)./plotScale,capacityUsedPerStack,'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Capacity Used (%)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 6 (Temperatures) ---
subplot(3,4,6)
plot(result(:,1)./plotScale,result(:,2)-273,'k-','LineWidth',1); hold on
plot(result(:,1)./plotScale,result(:,3)-273,'b-','LineWidth',1)
plot(result(:,1)./plotScale,result(:,6)-273,'m-','LineWidth',1)
plot(result(:,1)./plotScale,result(:,11)-273,'r-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Temperature (°C)')
set(gca,'fontsize',14,'LineWidth',1)
legend('T_s','T_{e,out}','T_{e,in}','T_r','Location','best')
grid on

% --- Subplot 7 (Potential) ---
subplot(3,4,7)
plot(result(:,1)./plotScale,result(:,4),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Potential (V)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 8 (Current Density) ---
subplot(3,4,8)
plot(result(:,1)./plotScale,result(:,5)./10000,'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Current Density (A cm^{-2})')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 9 (Inlet Temperature) ---
subplot(3,4,9)
plot(result(:,1)./plotScale,result(:,6)-273,'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Inlet Electrolyte Temp (°C)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 10 (H2 Production) ---
subplot(3,4,10)
plot(metrics(:,1)./plotScale,metrics(:,2),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('H_2 Production (kg)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 11 (H2 Rate) ---
subplot(3,4,11)
plot(metrics(:,1)./plotScale,metrics(:,3),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('H_2 Rate (kg h^{-1})')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% --- Subplot 12 (Gas Fraction) ---
subplot(3,4,12)
plot(metrics(:,1)./plotScale,metrics(:,4),'k-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Gas Fraction')
ylim([0 1])
set(gca,'fontsize',14,'LineWidth',1)
grid on


timeLoopSave(timeLoop,:) = [j stackVoltage];

end

%%
figure(2)
% Stack Solid
plot(result(:,1)./plotScale,result(:,2)-273,'m-','LineWidth',1); hold on
% Stack electrolyte outlet T
plot(result(:,1)./plotScale,result(:,3)-273,'k-','LineWidth',1)
% Stack electrolyte inlet T
plot(result(:,1)./plotScale,result(:,6)-273,'b-','LineWidth',1)
% Stack reservoir internal T
plot(result(:,1)./plotScale,result(:,11)-273,'r-','LineWidth',1)
xlabel(plotXlabel)
ylabel('Temperature (°C)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

% yyaxis right
% %V_r_save(1:end+1) = [0 V_r_save];
% fancoolpower = metrics(:,5);
% fancoolpower(1:end+1) = [0;fancoolpower];
% plot(result(:,1)./60,fancoolpower,'-','lineWidth',1)
% ylabel ('Fan Power [W]')

legend('T_s','T_{e,out}','T_{e,in}','T_r','Location','best')

yyaxis right
plot(result(1:end,1)./plotScale,result(:,8),'k--','LineWidth',1); hold on
legend('T_s','T_{e,out}','T_{e,in}','T_r','Power [W]','Location','best')
%set(gca,'Xscale','log')
set(gca,'Ycolor','k')
ylabel('Power [W]')

%xlim([40 50])




% 
% figure(3)
% plot(metrics(:,1),metrics(:,6)./result(1:end-1,8))

% % Optional experimental comparison data.
% if exist('exp_data.m','file') == 2
%     run('exp_data.m')
% 
%     plot(exp_T_stack_outlet(:,1)./60,exp_T_stack_outlet(:,2),'ko','LineWidth',1)
%     plot(exp_T_res(:,1)./60,exp_T_res(:,2),'ro','LineWidth',1)
%     plot(exp_T_e_in(:,1)./60,exp_T_e_in(:,2),'bo','LineWidth',1)
% 
%     legend('T_s','T_{e,out}','T_{e,in}','T_r','exp - T_s','exp - T_r','exp - T_{e,in}','Location','best')
% else
%     legend('T_s','T_{e,out}','T_{e,in}','T_r','Location','best')
% end

%% ============== Current Voltage ================= %% 
figure(3)
% Stack Solid
plot(result(1:end,1)./plotScale,result(:,9),'k-','LineWidth',1); hold on
% Stack electrolyte outlet T
ylabel('Current (A)')
yyaxis right
plot(result(1:end,1)./plotScale,result(:,10),'b-','LineWidth',1)

xlabel(plotXlabel)

ylabel('Potential (V)')
set(gca,'fontsize',14,'LineWidth',1)
grid on

%result(iteration+1,:) = [t T_s T_e potential j T_e_in h2P Power stackCurrent stackVoltage T_r];
%metrics(iteration,:) = [t h2P h2ProductionRateKg_hour alpha FanPower];


%% ============ efficiency per stack ============= %%

% % % --- Subplot 4 (Efficiency Per Stack) ---
% efficiencyPerStack = (1.48./(result(:,10)./N_cells)).*100;
% plot(result(:,1)./plotScale,efficiencyPerStack,'k-','LineWidth',1)
% xlabel(plotXlabel)
% ylabel('Efficiency Per Stack (%)')
% set(gca,'fontsize',14,'LineWidth',1)
% grid on

% 
% % --- Subplot 5 (Capacity Used Per Stack) ---
% plot(result(:,1)./60,capacityUsedPerStack,'k-','LineWidth',1)
% xlabel('Time (min)')
% ylabel('Capacity Used (%)')
% set(gca,'fontsize',14,'LineWidth',1)
% grid on
% 
