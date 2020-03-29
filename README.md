# CircuitModelLi-ionBattery
<pre>
The matlab scripts implement the Newman's P2D model to simulate a lithium ion battery cell.

2020-03-24   
Contact: zeyang.geng@chalmers.se  
         siyang.wang@mdh.se  
Folder  
 |  
 |- Main.m  
 |   |  
 |   |- Initiate the state of charge for each electrode with a number between 0 and 1.   
 |   |- Set up the simulation time and voltage limitations. The simulation will stop when it either reaches the maximum time or the voltage limitations.  
 |   |- Chose the time step. In the current version a constant time step is used and it is possible to implement variable time step if needed.  
 |   |- Chose the space discretization order.   
 |   |   |   
 |   |   |- The 4th order discretization uses SBP_variabl_4.m and SBP4.m  
 |   |   |- The 2nd order discretization uses SBP_variabl_2.m and SBP2.m   
 |   |   |- The 4th order discretization gives better result especially when the system is unstable (voltage profile changes rapidly).   
 |   |   +- The 4th order discretization requires minimum 8 meshing elements (9 grid points).  
 |   |     
 |   |- Set the current. Right now only constant current is used. Discharge current is positive and charge current is negative.  
 |   |- Setup_parameters (See explanations below)  
 |   |- SpatialDiscretization (See explanations below)  
 |   |- Initialization (See explanations below)     
 |   |- Simulation_loop (See explanations below)     
 |   +- Plot the discharge voltage profile.  
 |   
 |- Setup_parameters.m  
 |   |  
 |   |- The simulation is performed at 20 degree C.  
 |   |- The geometry of the cell (except the electrode area) and the electrochemical parameters used in this script are from "Schmalstieg, Johannes, et al. 2018".  
 |   |- Most of the electrochemical parameters are concentration dependent in reality however only constant numbers are used in this script.  
 |   |- The voltage profile for the electrodes are look up tables from Graphite_OCV.txt and NMC_OCV.txt  
 |   |- The number of the mesh element in the cell is 10 in each domain (11 grid points in each domain) by default.   
 |   |- The number of the mesh element in the particle is 30 by default.  
 |   +- A finer mesh increases both the accuracy and the calculation time.  
 |   
 |- SpatialDiscretization.m  
 |   |  
 |   +- The spatial derivatives are approximated by summation by parts (SBP) finite difference operators. Boundary and material interface conditions are imposed by the simultaneous approximation term (SAT). 
 |     
 |- Initialization.m      
 |   |  
 |   |- The initial electrolyte concentration is 1000 mol/m^3 (1 mol/dm^3).  
 |   |- The electrode concentration is initialized based on the initial condition in Main.m and thus the electrode potential.  
 |   |- The electronic resistance is initialized with the electronic conductivity and it is not updated in the later simulation.  
 |   |- The electrolyte resistance is initialized with the electrolyte conductivity before a concentration gradient is built up.     
 |   |- The charge transfer resistance is initialized with a linearization with the Butler-Volmer equation first.  
 |   |- Since the current density is different in each branch and it will affect the charge transfer resistance, a few initial loops are needed to stabilize the initial condition.   
 |   +- When converting the current density (per electrode area) to the local current density (per active surface area), the two branches at each boundary only take half of the space.  
 |  
 |- Simulation_loop.m  
 |   |  
 |   |- Set the start time to be zero and create matrixes/vectors to store the data.  
 |   |- Solve the electrolyte and electrode concentrations with the finite difference method.  
 |   |- The electrolyte concentration is a 1*n row vector where n is the total number of the grid points in the cell (1*33 by default).   
 |   |- The electrode concentration is a m*n matrix where m is the number of grid points in the particle (31*33 by default). There is one matrix for each electrode.  
 |   |- The first row in the electrode concentration matrix is the center concentration. The last row in the electrode concentration matrix is the surface concentration.  
 |   |- The electrode potential at the equilibrium state is calculated based on the average concentration but it is not used in the simulation and only for analysis.    
 |   |- The resistance matrix is updated based on the new concentration distribution and then used to solve the new current distribution.  
 |   +- Check the stop criteria in each step and give an indication about the stop reason.    
 |     
 |- SBP*.m  
 |   |  
 |   +- Used to generate the Laplacian operators.  
 |     
 +- *_OCV.txt  
     |  
     +- Loop up tables for the electrode potentials  
 <pre>


