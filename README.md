# electrolyserStackModel
A dynamic lumped model of an electrolyser stack (PEM and Alkaline) allowing prediction of stack and electrolyte/water temperatures as well as current-voltage output from power input. Implicit and explicit time integration to allow coupling with renewable energy time series.
1. To run with constant power input, set sim.readExternalPower = 'off'
2. Select power input 'Power', number of stacks 'N_stacks', desired stack temperature 'T_set'.
3. To change stack starting temperature, external air temperature etc go to lines 210 - 213.
1. To run only change simulation settings between 37 - 61.
2. Plotting of figures lines 620 +

![Alt text](https://github.com/DNiblett/electrolyserStackModel/blob/main/electrolyserStackModelExample.pdf "optional image title")
