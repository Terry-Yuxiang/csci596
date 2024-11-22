# CSCI 596 Final Project

## Dynamic Visualization of Temperature and Pressure Fields in Fluid Systems

This project implements a molecular dynamics (MD) simulation to dynamically visualize the temperature and pressure distribution in a fluid system under various conditions.

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