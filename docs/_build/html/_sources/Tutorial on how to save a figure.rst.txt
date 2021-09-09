Tutorial on how to save a figure
********************************

Figures of the reactor setup or the influent characteristics can be saved using the savename parameter.

.. code-block:: python

   import biops
   r = biops.ifas.Reactor()
   r.view_process(savename='path_and_name_of_figure')
   r.view_influent(savename='path_and_name_of_figure')

If you run a simulation and then plot the results using matplotlib.pyplot, you can save the figures using plt.savefig. Here is an example where we have run a biofilm simulation and the data is stored in the variable r.biofilmlogger.
The last statement in the code saves the figure.

.. code-block:: python
   
   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(figsize=(20/2.54, 17/2.54))
   ax.plot(r.biofilmlogger['1']['Conc'][20.0]['thickness']*10**6, r.biofilmlogger['1']['Conc'][20.0]['S_O2'])
   plt.savefig('path_and_name_of_figure')

