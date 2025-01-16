# 2D Ising Model

Simulation of the 2D Ising Model using the Metropolis Monte Carlo algorithm.

## Compilation Instructions
```bash
gfortran lanier_ising.f90 -o ising
./ising

Outputs:
IdealEnergyVsTemperature.txt: Average energy as a function of temperature.
IdealMagnetizationVsTemperature.txt: Average magnetization as a function of temperature.
IdealEnergyVsMagneticField.txt: Average energy as a function of the external magnetic field.
IdealMagnetizationVsMagneticField.txt: Average magnetization as a function of the external magnetic field.
IdealSpecificHeatVsTemperature.txt: Specific heat as a function of temperature.
