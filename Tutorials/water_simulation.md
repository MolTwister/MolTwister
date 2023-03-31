# Water simulation, thermalization, RDF, MSD, VACF, and VDOS

## Introduction

In this tutorial we will create a single water molecule, then make copies of this to create a water phase. Subsequently, we assign force field parameters, perform an energy optimization, and then a simulation.

Once the simulation is done, we will plot the temperature as a function of time to determine the point of thermalization. This is followed by calculating

* the radial distribution function,
* the mean square displacement and self-diffusion,
* the velocity auto correlation function,
* and the vibraional density of states.

## Create a water molecule

Invoke

    add atom Hw1 at 0 0 0

to add a hydrogen atom called Hw1. Then add the oxygen, Ow, by

    add atom Ow at 0.9572 0 0

where the O-H bond lenght is 0.9572 Å. Finally, we want to add an atom at a 120 degree angle. To learn how this is done, invoke

    help
    help atom

Looking at this, we see that

    add atom Hw2 from bond 0 1 angle 104.52 0 dist 0.9572

does what we want. We should also set appropriate boundaries to properly view the atom. Lets create a boundary that extends 10 Å in each direction from the coordinate origin by

    mod userdefpbc -10 10 -10 10 -10 10
    set userdefpbc on

and then lets make sure we view the entire boundary by invoking

    genbonds

This can be stored for safe keeping to a pdb and xyz file by

    print pdb > single_water.pdb
    print xyz > single_water.xyz

## Assigning force field parameters

We can now assign force field parameters to the newly created water molecule. We could enter the commands directly, but we may need to assign the force fields several times during a study. Hence, we instead create a file called 'water_ff.script' using our favourite text editor (e.g., from within MolTwister you can call vim or nano). The file should contain

    mod mass name Hw1 to 1.00784;
    mod mass name Hw2 to 1.00784;
    mod mass name Ow to 15.999;

    add mdbond Ow Hw1 Harm all_r_less_than 1.1 1882.8 0.9572;
    add mdbond Ow Hw2 Harm all_r_less_than 1.1 1882.8 0.9572;
    add mdangle Hw1 Ow Hw2 Harm all_r_less_than 1.1 230.12 104.52;

which defines TIP3P water molecules.

## Clear all contents

Many times we want to perform subtasks in MolTwister, where we end up with a script, pdb files, or similar that can be reloaded later. After creating these, we often want to delete everything we did from the MolTwister memory. This is done in two steps. First we delete what we see in the 3D view by

    sel all
    del sel

Subsequently, we delete all force field parameters by

    del ff

On some systems, it is also possible to open several instances of MolTwister, which can be useful in some situations, but here we will proceede with clearing all contents.

## Investigate the force field

Note that for the below Python scripts to work, MatplotLib must be installed on the system, which can be done by invoking

    pip install matplotlib

Python must be installed for MolTwister to work, so this should already be present on the system.

We can output data to show plots of the assigned force field potentials, but first we need to know the indices of the assigned force fields, which can be done by invoking

    list ff

Invoke

    help moldyn

to investigate the possible plots that can be made (i.e., bondforceprofile, angleforceprofile, dihedralforceprofile, and nonbondforceprofile). We can for example look at the Lennard-Jones potential between the water oxygens, which is the non-bonded potential at index 0. This can be done by

    moldyn ff nonbondforceprofile 0 2.9 10.0 100 > OwOw_LJ.dat

where the result is stored in OwOw_LJ.dat file. This can be plotted using the following Python script, with name 'plot_non_bonded.py':

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("OwOw_LJ.dat", "r")
startOfList = False
for line in f:
    if startOfList:
        s = line.split()
        if len(s) >= 3:
            listX.append(float(s[0])) # Extract r in Å
            listY.append(float(s[2])) # Extract U in kJ/mol

    if "Potential" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='r [Å]', ylabel='U [kJ/mol]', title='Ow-Ow Lennard-Jones potential')
plt.show()
```

To plot the Lennard-Jones potential, invoke

    python plot_non_bonded.py

or

    python3 plot_non_bonded.py

depending on your system.

## Create liquid water

We now want to use the single_water.pdb file that we created and construct a water phase containing 512 water molecules, randomly placed, without collision. Lets first calculate the approximate box size if the density is 997 kg/m^3, with two hydrogen masses of 1.0078 and one oxygen mass of 15.999:

    calculate volumefromdensity 997 1.0078,1.0078,15.999 512

This yields a side length of approximately 25 Å for a cubic simulation box.

We should now make sure we start with a clean MolTwister:

    sel all
    del sel
    del ff

First we load the water molecule:

    load single_water.pdb

and say 'del' to delete anything previously there. Subsequently, select all atoms and make 512 copies that are placed randomly within a cubic box with 25 Å side lengths and a minimum distance between molecular atoms of 1.5 Å:

    copy sel random 25 25 25 512 1.5

Check if we have managed to create 512 molecules by selecting all water oxygens (Ow) and counting them:

    sel none
    sel name Ow
    measure count sel

The MD simulator only accepts systems with geometric centers placed at the origin. Hence, we need to displace the system to this position by invoking

    mod atompos all geomcent 0 0 0
    autoscale

We now save the system to a PDB file by

    print pdb > water.pdb

## Set up an appropriate periodic boundary condition (PBC)

We now want to set up appropriate periodic boundary conditions (PBC) for the simulation. We know that all atoms are within a box of 25 Å from the above construction, but we can also measure the PBC by

    measure pbc

Then, we here add a sufficient amount on each side of the PBC to make sure that molecules do not come to close to the molecules on the other side of the PBC. Thus, we invoke

    mod userdefpbc -12.5 12.5 -12.5 12.5 -12.5 12.5

to set a 25 by 25 by 25 Å boundary box and we invoke

    set userdefpbc on

to activate it. It can be smart to create a file called 'pbc.script' and insert the two last commands, just so that we do not have to manually adjust the PBC every time, and instead just write

    load pbc.script

when needed.

## Perform energy optimization

The molecules were placed randomly in the liquid phase that we created, so it is a good idea to minimize forces between the atoms by performing an energy optimization. We will do this by using the default maximum number of steps for the optimization and the default learning rate, but we will set that the algorithm will stop if we reach an energy change of less than 10 kJ/mol.

To see the default values, invoke

    moldyn cfg get

To set the accuracy of 10 kJ/mol, invoke

    moldyn cfg set optaccuracy 10

We now start from a fresh MolTwister:

    sel all
    del sel
    del ff

Next, invoke

    load water.pdb
    load pbc.script
    load water_ff.script

Then, start the energy optimization by

    moldyn cfg set outstride 10
    moldyn optimizeenergy

and wait until it terminates.

This will produce a file called 'traj.dcd' and 'out.txt'. 

By creating a Python script called 'plot_energy.py', containing
```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("out.txt", "r")
startOfList = False
for line in f:
    if startOfList and "Searching" not in line:
        s = line.split()
        if len(s) >= 6:
            listX.append(float(s[0])) # Extract time step column
            listY.append(float(s[5])) # Extract energy column

    if "Step," not in line and "Step" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='Step', ylabel='U [kJ/mol]', title='Energy optimization')
plt.show()
```

we can view the energy plot of the optimization process by

    python plot_energy.py

If we are happy with the results, we can export the last frame of the optimization process to a PDB file. First load the trajectory by

    load traj.dcd

where you can view the optimization process by clicking the 3D view and using the left and right arrow keys to step through. You can also enable fog to better the depth visualization by

    set fog on

We now export the optimizated frame (frame 499 in our case) by

    print pdb frame 499 > initial.pdb

We will use this as our starting frame for the simulations.

If you wish, rename the traj.dcd and out.txt files from the energy optimization for safe keeping. They will be deleted in a later step.

## Perform molecular dynamics simulation

We start again by clearing MolTwister memory:

    sel all
    del sel
    del ff

Now remove any left-over traj.dcd and out.txt files.

We reload the initial system, the PBC, and the force field parameters:

    load initial.pdb
    load pbc.script
    load water_ff.script

We configure to run 100 000 time steps at 0.5 fs, yielding a total of 50 ps, and we set the stride to 10 frames:

    moldyn cfg set timestep 0.5
    moldyn cfg set timesteps 100000
    moldyn cfg set outstride 10

Now we can run the simulation as

    moldyn run

On a relatively new CPU this simulation will take a few hours and produce a traj.dcd trajectory file and an out.txt output file.

Once the simulation is done, we can load the trajectory file by

    load traj.dcd

and analyze the simulation by using the left and right keyboard arrows. Invoke help to see how to adjust the number of steps per keypress, as well as how to show/hide bonds across boundaries, and switch between perspective and orthographic view.

## Plot temperature as function of simulation time

The temperature as function of time step is available in 'out.txt' and we will use Python to plot the results. Create the file 'plot_temp.py' and enter the following script:

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("out.txt", "r")
startOfList = False
for line in f:
    if startOfList and "Searching" not in line:
        s = line.split()
        if len(s) >= 3:
            listX.append(float(s[0])) # Extract time step column
            listY.append(float(s[2])) # Extract temperature column

    if "Timestep," not in line and "Timestep" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='Step', ylabel='T [K]', title='Temperature')
plt.show()
```

Now invoke

    python plot_temp.py

which should show clearly the thermalization period of the system. This will be important later to determine at which points in the simulation we find thermalized data to perform calculations on. In this tutorial we will assume that the system was thermalized after 8000 time steps, which corresponds to frame index 800 in the DCD file, since we output every 10 timestep.

## Calculate and plot a pair correlation function

We now want to calculate the pair correlation function, or radial distribution function (RDF), between oxygen atoms of water (Ow). We will use the data between frame 800 and 9999, which is the last frame of the DCD file, and create a plot of 100 points between 0.4 and 12 Å. This can be done by first making sure that we have loaded the system (i.e., the PDB file in our case, not the entire trajectory) and then invoking

    calculate paircorrelation traj.dcd 800 9999 Ow Ow 100 0.4 12.0 > pair_corr.dat

which can be viewed using the following python script

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("pair_corr.dat", "r")
startOfList = False
for line in f:
    if startOfList:
        s = line.split()
        if len(s) >= 4:
            listX.append(float(s[0])) # Extract r in Å
            listY.append(float(s[3])) # Extract RDF

    if "-----------" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='r [Å]', ylabel='RDF', title='Radial distribution function')
plt.show()
```

## Calculate and plot mean square displacement

We want to calculate the mean square displacement (MSD). However, to use this command we first need to add residue names to our PDB file. We can now start with a clean MolTwister:

    sel all
    del sel
    del ff

Then, load the system

    load initial.pdb

We assign the residue names to each atom name by

    mod resname name Ow to water
    mod resname name Hw1 to water
    mod resname name Hw2 to water

where the water residue name was set to 'water'. We now generate a new PDB file containing residue names by

    print pdb > initial.pdb

Since the mean square displacement attempts to trace the molecules movements over time we need to unwrap the molecules accross the periodic boundaries before we can start the calculation. This is done by invoking

    dcd unwrap traj.dcd

which will produce a file called 'traj_mtunwrap.dcd'.

To test the unwraped system, clear the MolTwister memory by

    sel all
    del sel
    del ff

and then invoke

    load initial.pdb
    load pbc.script
    load traj_unwrap.dcd

to load the system into the 3D viewer. In the 3D viewer, hit Alt+O to switch into orthographic view and then use the arrow keys to scroll though the frames, checking that the molecules move appropriately outside the periodic boundaries.

Once this is done we can calulate the mean square displacement by

    calculate msd traj_mtunwrap.dcd 800 9999 resname water > msd.dat

Note that the 'traj_mtunwrap.dcd' file does *not* need to be loaded into MolTwister memory in order to run the above command.

which can be plotted using the following Python script:

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("msd.dat", "r")
startOfList = False
for line in f:
    if startOfList:
        s = line.split()
        if len(s) >= 2:
            listX.append(float(s[0])) # Extract step
            listY.append(float(s[1])) # Extract MSD

    if "Index" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='step', ylabel='M [Å^2]', title='Mean square displacement')
plt.show()
```

## Calculate and plot velocity auto correlation function

We wish to calculate the velocity auto correlation function (VACF). We then need to specify the time between each frame in fs, which is 5 in our case (i.e., 0.5 fs time step and output every 10 steps), as well as the number of time steps to include in the plot, which we will chose to be 100. The velocity auto correlation function should use the unwrapped system, as was also the case for the mean square displacement. First select the atoms to be included by

    sel all

and then do the calculation by invoking

    calculate vacf traj_mtunwrap.dcd 800 9999 5 100 sel > vacf.dat

and then plot using the below python script:

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("vacf.dat", "r")
startOfList = False
for line in f:
    if startOfList:
        s = line.split()
        if len(s) >= 3:
            listX.append(float(s[0])) # Extract t in fs
            listY.append(float(s[2])) # Extract VACF

    if "vacf" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

ax.set(xlabel='t [fs]', ylabel='nVACF', title='Normalized auto correlation function')
plt.show()
```

## Calculate and plot vibrational density of states

The vibrarional density of states (VDOS) can be calculated both using the unwrapped and wrapped systems, where we in our case will use the unwrapped one, since we have it available. We specify the starting frame for the calculation to 800, as well as the size of the fast Fourier transform (FFT) being applied, which is also the required number of available frames after the starting frame. We specify 10, which results in a 2^10 = 1024 FFT. We start the calculation by

    sel all
    calculate vdos traj_mtunwrap.dcd 800 10 0.5 sel > vdos.dat

We can plot the result using the following Python script:

```python
import matplotlib.pyplot as plt

listX = []
listY = []
f = open("vdos.dat", "r")
startOfList = False
for line in f:
    if startOfList:
        s = line.split()
        if len(s) >= 6:
            listX.append(float(s[4])) # Extract wavelength, lambda, in nm
            listY.append(float(s[5])) # Extract VDOS

    if "inf" in line:
        startOfList = True

fig, ax = plt.subplots()
ax.plot(listX, listY)

plt.xlim([1780,8000])
plt.ylim([0,0.00014])
ax.fill_between(listX, listY, alpha=0.2)
ax.set(xlabel=r'$\lambda$ [nm]', ylabel='VDOS', title='Vibrational density of states')
plt.show()
```
