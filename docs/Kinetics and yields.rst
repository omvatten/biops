Kinetics and yields
*******************

Currently, biops calculates the activities of ammonia-oxidizing bacteria (AOB), nitrite-oxidizing bacteria (NOB), anammox, comammox, and ordinary heterotrophic bacteria (OHO).
The growth and substrate consumption rates are described using Monod kinetics (Haldane kinetics for commamox). Default kinetic constants are set, but they can be changed with the change_constants function.

.. code-block:: python

   r.change_constants(func_group=None, **kwargs)

Assume we have a reactor objected, r, and we want to change the maximum growth rate of AOB to 3 d-1. Then, we can use the following code:

.. code-block:: python

   r.change_constants(func_group='aob', mu_max=3)

The following list have all the kinetics constants and their default values:

AOB:

- mu_max = 2.05 (maximum growth rate, per day)
- b = 0.13 (decay rate, per day)
- K_O2 = 0.6 (half-saturation constant for oxygen, g/m3)
- K_NH4 = 2.4 (half-saturation constant for ammonium, g/m3)
- fs = 0.054 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- CHON = [5,7,2,1] (assumed elemental composition of bacterial cells, this is related to the calculation of the yield coefficient)
- fI = 0.2 (fraction of biomass that becomes inert material during decay)

NOB:

- mu_max = 1.45 (maximum growth rate, per day)
- b = 0.06 (decay rate, per day)
- K_O2 = 0.4 (half-saturation constant for oxygen, g/m3)
- K_NO2 = 0.5 (half-saturation constant for nitrite, gN/m3)
- fs = 0.07 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- CHON = [5,7,2,1] (assumed elemental composition of bacterial cells, this is related to the calculation of the yield coefficient)
- fI = 0.2 (fraction of biomass that becomes inert material during decay)

anammox:

- mu_max = 0.08 (maximum growth rate, per day)
- b = 0.003 (decay rate, per day)
- K_O2 = 0.01 (half-saturation constant for oxygen, g/m3)
- K_NO2 = 0.05 (half-saturation constant for nitrite, gN/m3)
- K_NH4 = 0.07 (half-saturation constant for ammonium, gN/m3)
- fs = 0.051 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- CHON = [1,1.74,0.31,0.2] (assumed elemental composition of bacterial cells, this is related to the calculation of the yield coefficient)
- fI = 0.2 (fraction of biomass that becomes inert material during decay)

comammox:

- mu_max = 0.14 (maximum growth rate, per day; based on Kits et al. Nature 549, 269, 2017; 14.8 umolN/mgProt.h, 400 mgProt/molN)
- b = 0.003 (decay rate, per day)
- K_O2 = 0.4 (half-saturation constant for oxygen, g/m3)
- K_NH4 = 0.012  (half-saturation constant for ammonium, gN/m3; Kits et al. Nature 549, 269, 2017; 0.84 uM)
- Ki_NH4 = 3.44 (inhibition constant for ammonium; Sakoula et al. ISME 15, 1010, 2020; 246 uM)
- fs = 0.02 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- CHON = [5,7,2,1] (assumed elemental composition of bacterial cells, this is related to the calculation of the yield coefficient)
- fI = 0.2 (fraction of biomass that becomes inert material during decay)

OHO:

- mu_max_O2 = 6 (maximum growth rate, per day)
- mu_max_NOx = 4.8 (maximum growth rate, per day)
- b = 0.62 (decay rate, per day)
- K_s = 20 (half-saturation constant for COD, g/m3)
- K_O2 = 0.2 (half-saturation constant for oxygen, g/m3)
- K_NOx = 0.3 (half-saturation constant for nitrite and nitrate, gN/m3)
- fs_O2 = 0.67 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- fs_NOx = 0.67 (fraction electron diverted to biomass, this is related to the calculation of the yield coefficient)
- CHON = [5,7,2,1] (assumed elemental composition of bacterial cells, this is related to the calculation of the yield coefficient)
- fI = 0.2 (fraction of biomass that becomes inert material during decay)


