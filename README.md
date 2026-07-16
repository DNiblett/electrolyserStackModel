# Dynamic Electrolyser Stack Model

A zero-dimensional, lumped thermal–electrochemical model for simulating the dynamic operation of alkaline, proton-exchange-membrane (PEM), and anion-exchange-membrane (AEM) electrolyser stacks in MATLAB.

The model resolves stack electrochemistry, heat transfer, electrolyte circulation, reservoir temperature, gas hold-up, hydrogen production, balance-of-plant loads, variable electrical power input, and an optional empirical voltage-degradation model.

## Features

- Alkaline, PEM, and AEM parameter sets
- Dynamic power input from constant, synthetic, or external profiles
- Galvanostatic operating point calculated from imposed electrical power
- Temperature-dependent exchange current density, conductivity, and electrolyte properties
- Lumped stack, electrolyte, inlet-line, and reservoir temperatures
- Explicit or implicit thermal time integration
- Two-phase gas hold-up using a drift-flux approximation
- Bubble-coverage and bubble-resistance corrections
- Hydrogen production and water-consumption calculations
- Minimum-load clipping
- Stack start-up, shutdown, heater, cooling, and reservoir refill logic
- Optional current-history-dependent degradation model


  ## Model overview

<p align="center">
  <img src="figures/stackModelOutput_1.png" width="700">
</p>

## Model Structure

The model represents the electrolyser system using four principal temperatures:

- `T_s`: stack solid temperature
- `T_e`: stack electrolyte outlet temperature
- `T_r`: reservoir electrolyte temperature
- `T_e_in`: stack inlet electrolyte temperature after pipe losses and cooling

The electrochemical cell voltage is calculated from:

$$
V_{\mathrm{cell}} = E_{\mathrm{rev}} + \eta_{\mathrm{anode}} + \eta_{\mathrm{cathode}} + \eta_{\mathrm{ohmic}} + \eta_{\mathrm{bubble}} + \Delta V_{\mathrm{deg}}
$$

The stack current density is obtained iteratively so that:

$$
P_{\mathrm{input}} = N_{\mathrm{stacks}}N_{\mathrm{cells}}A_{\mathrm{cell}}jV_{\mathrm{cell}}
$$

## Requirements

- MATLAB
- No additional MATLAB toolboxes are required for the core model
- Optional external power data must be supplied as a `.mat` or compatible data file

The supplied example expects external wind-power data in variables similar to:

```matlab
time_s
power_MW
```

## Running the Model

1. Download or clone the repository.
2. Open MATLAB in the repository folder.
3. Open `dynamicElectrolyserStackModel_V3_degradation.m`.
4. Edit the simulation settings near the top of the file.
5. Run the script.

Example:

```matlab
sim.electrolyserType = 'PEM';
sim.readExternalPower = 'off';
sim.isothermal = 'on';

Power = 1e6;
sim.timeStep = 2;
sim.endTime = 60*60*100;
```

## Main Simulation Settings

### Electrolyser type

```matlab
sim.electrolyserType = 'Alkaline';
```

Available options:

```text
Alkaline
PEM
AEM
```

### Power input

For constant power:

```matlab
sim.readExternalPower = 'off';
Power = 1e6;
```

For external time-dependent power:

```matlab
sim.readExternalPower = 'on';
externalPowerFile = 'Wind_power_2020_AR1_1s.mat';
```

Interpolation can be enabled using:

```matlab
sim.interpolatePower = 'on';
```

### Numerical settings

```matlab
sim.temporalScheme = 'implicit';
sim.alphaScheme = 'implicit';
sim.timeStep = 2;
sim.endTime = 60*60*100;
sim.writeTime = 10;
```

### Minimum-load clipping

```matlab
sim.minLoadClipping = 'on';
```

Power below the technology-specific minimum load is set to zero.

### Isothermal operation

```matlab
sim.isothermal = 'on';
```

When enabled, the stack, electrolyte, inlet, and reservoir temperatures are fixed at the selected set temperature.

### Two-phase model

```matlab
twoPhase.active = 'on';
```

The model estimates gas hold-up using a lumped drift-flux relation and modifies electrochemical performance through bubble coverage and additional resistance.

## Degradation Model

The degradation model adds a cumulative cell-voltage penalty:

$$
V_{\mathrm{cell,degraded}} = V_{\mathrm{cell,fresh}} + \Delta V_{\mathrm{deg}}
$$

The total degradation can include:

$$
\Delta V_{\mathrm{deg}} = \Delta V_{\mathrm{steady}} + \Delta V_{\mathrm{low-load}} + \Delta V_{\mathrm{ramp}} + \Delta V_{\mathrm{cycle}}
$$

### Steady degradation

$$
\dot V_{\mathrm{steady}} = r_{\mathrm{ref}}\left(\frac{j}{j_{\mathrm{rated}}}\right)^m
$$

### Ramp degradation

$$
\Delta V_{\mathrm{ramp}} = k_{\mathrm{ramp}}\max\left(\frac{|j-j_{\mathrm{previous}}|}{j_{\mathrm{rated}}}-\epsilon_{\mathrm{ramp}},0\right)
$$

### Start-stop degradation

$$
\Delta V_{\mathrm{cycle}} = N_{\mathrm{start}}\Delta V_{\mathrm{start}} + N_{\mathrm{stop}}\Delta V_{\mathrm{stop}}
$$

The degradation model is empirical and intended for comparing operating strategies, renewable-following profiles, relative lifetime impacts, and system-level techno-economic studies. It should not currently be interpreted as a mechanistic prediction of catalyst dissolution, membrane thinning, passivation, corrosion, or PTL oxidation.

## Example Degradation Parameters

```matlab
switch lower(sim.electrolyserType)

    case 'alkaline'
        degradation.steadyRate      = 3.2e-6/3600;
        degradation.currentExponent = 1.5;
        degradation.lowLoadRate     = 2.0e-6/3600;
        degradation.rampLoss        = 0.25e-6;
        degradation.startLoss       = 1.0e-6;
        degradation.stopLoss        = 1.0e-6;

    case 'pem'
        degradation.steadyRate      = 5.0e-6/3600;
        degradation.currentExponent = 2.0;
        degradation.lowLoadRate     = 5.0e-6/3600;
        degradation.rampLoss        = 0.5e-6;
        degradation.startLoss       = 15e-6;
        degradation.stopLoss        = 15e-6;
end
```

These values are illustrative. Event losses depend on the assumed cycling frequency and should be documented alongside any reported simulation results.

## Outputs

Typical columns in the main result matrix are:

| Column | Quantity | Units |
|---:|---|---|
| 1 | Time | s |
| 2 | Stack solid temperature | K |
| 3 | Electrolyte outlet temperature | K |
| 4 | Cell potential | V |
| 5 | Current density | A m\(^{-2}\) |
| 6 | Electrolyte inlet temperature | K |
| 7 | Cumulative hydrogen production | kg |
| 8 | Supplied electrical power | W |
| 9 | Stack current | A |
| 10 | Stack voltage | V |
| 11 | Reservoir temperature | K |
| 12 | Total degradation voltage | V cell\(^{-1}\) |
| 13 | Steady degradation contribution | V cell\(^{-1}\) |
| 14 | Ramp degradation contribution | V cell\(^{-1}\) |
| 15 | Low-load degradation contribution | V cell\(^{-1}\) |
| 16 | Start-stop degradation contribution | V cell\(^{-1}\) |

The script also calculates hydrogen production rate, cumulative hydrogen production, gas fraction, cooling power, heater power, stack utilisation, and stack efficiency.

## External Power Profiles

External power data are loaded using:

```matlab
externalData = load(externalPowerFile);
powerData_Time = double(externalData.time_s);
powerData_Power = double(externalData.power_MW).*1e6;
```

Power is interpolated using a moving index rather than repeatedly searching the full dataset, which is important for long, high-resolution renewable-power profiles.

## Model Assumptions

- Spatially uniform, lumped stack
- Uniform current density across each cell
- Identical cells within a stack
- Common current for cells connected electrically in series
- Lumped gas hold-up
- Negligible pressure dynamics
- Empirical bubble effects
- Empirical degradation represented as an added voltage loss
- No detailed species crossover model
- No explicit mechanical compression model
- No spatial catalyst, membrane, separator, or PTL degradation

## Known Limitations

- Some material and electrochemical parameters are approximate.
- The AEM parameter set is less mature than the alkaline and PEM sets.
- Dynamic degradation coefficients are not universal and require calibration.
- External balance-of-plant models are simplified.
- The model does not yet distinguish reversible conditioning from irreversible degradation.
- At very low current density, logarithmic activation expressions require careful numerical treatment.
- High-current operation may require an additional degradation term rather than a single current exponent.
- The two-phase model is lumped and does not resolve bubble-size distributions or spatial gas transport.

## Recommended Validation

The model should be validated against:

1. steady-state polarisation curves
2. stack current and voltage under imposed power
3. measured stack and electrolyte temperatures
4. reservoir heating and cooling transients
5. hydrogen production rate
6. gas hold-up or outlet gas-flow measurements
7. constant-current degradation tests
8. repeated load-cycle and start-stop experiments

## Suggested Repository Structure

```text
.
├── dynamicElectrolyserStackModel_V3_degradation.m
├── README.md
├── data/
│   └── example_power_profile.mat
├── examples/
│   ├── constant_power_case.m
│   └── fluctuating_power_case.m
├── figures/
└── LICENSE
```

## Citation

When using this model in academic work, please cite the repository and associated publication once available.

Suggested temporary citation:

```text
Niblett, D. Dynamic Electrolyser Stack Model.
Newcastle University, 2026. GitHub repository.
```

## Author

**Dr Daniel Niblett**, **Dr Majid Rahgoshay**, **Professor Mohamed Mamlouk**

Newcastle University

Developed initially for the UKRI EPSRC Ocean Refuel project and subsequently updated for Ocean Refuel and EPSRC ION-H2.

## Contributing

Issues, validation datasets, bug reports, and model improvements are welcome.

Useful contributions include additional electrolyser parameter sets, experimentally calibrated degradation coefficients, alternative gas-hold-up models, improved property correlations, pressure and crossover models, numerical acceleration, and example renewable-power profiles.

## Disclaimer

This model is a research tool. It is not certified for electrolyser design, plant safety assessment, warranty prediction, operational control, or commercial performance guarantees.
