Tutorial on how to specify influent
***********************************

Constant influent
#################

Let's assume you influent has constant values over time. Then, we can specify the influent as a python dictionary.
Below, we specify an influent with the flow rate 10 m3/d, a concentration of organic compounds of 100 g/m3 COD, and an ammonium concentration of 50 g/m3 N.
All other parameters, such as S_NO2, S_NO3, X_OHO, X_AOB, X_NOB, X_AMX, X_CMX, X_I, are assumed to be zero in this case (but we can of course include a value for them as well).

.. code-block:: python

   influent = {'Q':10, 'S_s':100, 'S_NH4':50}

When we create a reactor object, we can add the influent specified above:

.. code-block:: python

   import biops
   r = biops.ifas.Reactor(influent=influent)

Influent that varies with time
##############################

Now let's assume we have a file with influent values that change over time. The influent data should be written in .csv file (comma-separated values). You can use Excel to create such a file. When you do Save as in Excel, choose CSV in the 'Save as type' field.
We need to have a column with the heading Time, which is the time in days, and a column with the heading Q, which is the flow rate in m3/d. 
Then, we can have optional columns with S_s, S_NH4, S_NO2, S_NO3, X_OHO, X_AOB, X_NOB, X_AMX, X_CMX, and X_I as headings. If a parameter is not present in the file, it is assumed to be zero.

.. list-table::
   :widths: 15 15 15 15 15
   :header-rows: 1

   * - Time
     - Q
     - S_s
     - S_NH4
     - S_NO3
   * - 0
     - 1200
     - 220
     - 28
     - 12
   * - 0.2
     - 1220
     - 180
     - 29
     - 11
   * - 0.4
     - 1100
     - 190
     - 32
     - 8
   * - 0.6
     - 1004
     - 188
     - 36
     - 9

Next, we want to load this file into python. We use pandas.

.. code-block:: python

   import pandas as pd
   influent = pd.read_csv('path_to_csv_file')

When we create a reactor object, we can add the influent specified above:

.. code-block:: python

   import biops
   r = biops.ifas.Reactor(influent=influent)

