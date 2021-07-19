Commands
********

Create a reactor
################

.. code-block:: python

   ifas.Reactor(influent='default', volume=4000, VSS=3000, DO=2, Xbulk=default_Xbulk, Sbulk=default_Sbulk, area=100, Lmax=400*(10**-6), Lbf=5*(10**-6), VSmax=25000, DBL=10*(10**-6))

Creates a Reactor object.

*influent* is a pandas dataframe with information about the influent to the reactor.
The 'default' influent is a subset of the dry weather influent for the IWA benchmark simulation model (BSM1, http://iwa-mia.org/benchmarking/).
However, the user can provide her/his own influent data. The data must contain a column with the heading 'Q' having flow data in m3/day.
Optionally, columns the headings 'Time' (in days), 'S_s' (biodegradable COD, g/m3), 'S_NH4' (gN/m3), 'S_NO2' (gN/m3), 'S_NO3' (gN/m3),
'X_OHO' (heterotrophic biomass, gVSS/m3), 'X_AOB' (gVSS/m3), 'X_NOB' (gVSS/m3), 'X_AMX' (anammox, gVSS/m3), 'X_CMX' (comammox, gVSS/m3), 'X_I' (inert biomass, gVSS/m3)
can also be provided.

*volume* is the volume of the reactor in m3.

*VSS* is the volatile suspended solids concentration in the reactor (gVSS/m3). Set to 0 if it is a biofilm reactor without suspended growth.

*DO* is the dissolved oxygen concentration in the reactor.

*Xbulk* is the concentrations of biomass components at the start of the simulation. The default is equal concentrations of each component.

*Sbulk* is the concentrations of substrate components at the start of the simulation.

*area* is the surface area for biofilm growth. Set to 0 if it is a suspended growth reactor without biofilm.

*Lmax* is the maximum biofilm thickness (m).

*Lbf* is the biofilm thickness at the start of the simulation (m).

*VSmax* is the maximum biomass concentration in the biofilm (gVSS/m3).

*DBL* is the diffusion boundary layer between the bulk liquid and the biofilm (m). A smaller thickness leads to lower diffusion resistance. 

Example: The code below will create a reactor object name *r*, using the default parameter inputs.

.. code-block:: python

   import biops
   
   r = biops.ifas.Reactor()

Add a compartment
#################

Assume we have created a reactor object named *r*. Now we want to add a compartment.

.. code-block:: python

   r.add_compartment(influent=None, volume='same', DO='same', Xbulk='same', Sbulk='same', area='same', Lmax='same', Lbf='same', VSmax='same', DBL='same')

The default inputs *same* assume that everything are the same as in the preceding reactor compartment.

The default *influent* is None. This assumes that the compartment does not receive any external influent, except for the water flowing from the preceding compartment.
Influent to this compartment can be species as a pandas dataframe or dictionary object with information about flow ('Q') and other parameters.

Add a recirculating flow
########################

Assume we have added a compartment to our reactor and now we want a recirculating flow from the second compartment to the first compartment.
This is, for example, the design of pre-denitrification systems in wastewater treatment.

.. code-block:: python

   r.add_recirculation(origin=None, target=None, Q=None)

*origin* is the compartment number from which the flow is coming (for example, 2).

*target* is the compartment number to which the flow is going (for example, 1).

*Q* is the flow rate in m3/d.   

View reactor configuration
##########################

Now, we have built our reactor system and want to have a look at it. 

.. code-block:: python

   r.view_process(savename=None)

if we specify *savename*, the resulting figure will be saved as a file.

View reactor influent
#####################

What is the actual influent to our reactor? Take a look with:

.. code-block:: python

   r.view_influent(savename=None)

if we specify *savename*, the resulting figure will be saved as a file.

Run a simulation
################

Now, we want to check how the reactor content (VSS composition, biofilm, effluent quality, etc.) changes after a certain time interval.

.. code-block:: python

   r.calculate(timestep=1, dt=0.1, set_iterations=1000, set_dts='auto')

*timestep* is the time duration that we will simulate (d).

*dt* is the resolution of the timesteps taken when data about influent characteristics are read and biomass compositions are calculated.

*set_iterations* is related to the calculations of substrate concentrations, which is done for each timestep (dt). 
This value is the maximum number of iterations done before moving forward in time (dt). A larger value might result in more accurate predictions but take longer time to complete.

*set_dts* is related to the the calculations of substrate concentrations. 'auto' means that the software will automatically determine the appropriate timestep.
If the resulting substrate concentration profiles are oscillating and unrealistic, set a low value here could solve the problem. Test, e.g. 10**-7 or smaller. 

Information about the influent, bulk water, and biofilm concentrations and activities are logged at the end of the timestep using in r.calculate(). The information is stored in the following data objects.

.. code-block:: python

   r.influentlogger
   
   r.bulklogger
   
   r.biofilmlogger
   
To simulate the performance of the reactor, the *calculate* method should be run in a loop. See specific instructions for this.
