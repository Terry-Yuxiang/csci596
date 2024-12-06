# CSCI 596 Final Project

## Dynamic Visualization of Molecular Temperature and Heat Transfer in Systems

This project focuses on implementing a molecular dynamics (MD) simulation to dynamically visualize temperature variations and heat transfer across molecular layers in a fluid system under various conditions. The primary goal is to observe and analyze the propagation of temperature changes and their effect on molecular behavior, providing an intuitive understanding of thermal dynamics.

## Project Overview

This project aims to simulate and visualize temperature transfer in a molecular system using Molecular Dynamics (MD). We create a closed system where:

1. **Temperature Gradient Setup**:
   - First layer: High temperature
   - Last layer: High temperature
   - Other layers: Low temperature
   - Observing heat transfer between these extremes
   
2. **System Size Analysis**:
   - Adjustable system size through unit cell parameters (X, Y, Z directions)
   - Enables simulation of different molecular quantities
   - Allows analysis of how system size affects heat transfer dynamics

3. **Key Observations**:
   - Temperature propagation patterns across molecular layers
   - Molecular motion and behavior during heat transfer
   - System equilibrium achievement process
   - Impact of system size on heat transfer rates

### Simulation Features
- Closed system simulation (constant volume, particle number)
- Customizable initial conditions through `md.in` file
- Real-time visualization of:
  - Molecular movements
  - Temperature distribution (Deep blue [coolest] → Light blue → White → Light red → Deep red [hottest])
  - Heat transfer pathways
---

## How to Run

### Prerequisites
Ensure you have the following installed:
- GCC Compiler (`gcc`)
- A Unix-like environment (Linux/macOS/Windows with WSL or MinGW)
- A terminal or command-line interface

### Steps to Compile and Run
1. **Navigate to the project directory**:
   ```bash
   cd /path/to/your/project
   Replace /path/to/your/project with the actual directory containing md_temp.c and md.in.
2.	**Compile the program**:
    ```bash
    gcc -o md_temp md_temp.c -lm
3.	**Prepare the input file (md.in)**:
Ensure a file named md.in exists in the same directory and follows the required format. Refer to the Input Parameters section below for details.
4.  **Run the program**:
    ```bash
    ./md_temp < md.in
5.	**Check the output**:
•	Results are displayed in the terminal.
•	Simulation data is also saved in results.txt in the same directory.

### Examples of `md.in` file
```bash
4 4 4         # Number of unit cells in the X, Y, and Z directions
0.8           # Density of particles in the system (Density)
1.5           # Initial temperature of the system (InitTemp)
0.005         # Time step size for the simulation (DeltaT)
10000         # Total number of simulation steps to execute (StepLimit)
100           # Number of steps between each data output (StepAvg)
```

### Visualization with PyMOL
1. **Generate PDB file**:
   ```bash
   python generate_water_matrix.py
   ```

2. **Open PyMOL and run simulation**:
   ```bash
   cd /path/to/project
   pymol
   run simulate_molecular_motion.py
   ```

3. **Generate video from frames**:
   ```bash
   cd frames
   ffmpeg -framerate 24 -i frame_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p simulation.mp4
   ```

The final video will be saved as `simulation.mp4` in the frames directory.


## Initial effect demonstration
Temperature scale: Deep blue (coolest) → Light blue → White → Light red → Deep red (hottest)
### Four layers of molecular temperature transfer 
![temp](./asset/show1.gif)


### Ten layers of molecular temperature transfer 
![temp](./asset/show2.gif)
