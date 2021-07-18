import numpy as np
import matplotlib.pyplot as plt
from . import ifas

aob = {}
aob['mu_max'] = 2.05 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
aob['b'] = 0.13 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
aob['K_O2'] = 0.6 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)
aob['K_NH4'] = 2.4 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)

nob = {}
nob['mu_max'] = 1.45 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
nob['b'] = 0.06 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
nob['K_O2'] = 0.4 #gN/m3 (Laureni et al. WR 154, 104-116, 2019)
nob['K_NO2'] = 0.5 #g/m3 (Laureni et al. WR 154, 104-116, 2019)

anammox = {}
anammox['mu_max'] = 0.08 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
anammox['b'] = 0.003 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_O2'] = 0.01 #gO2/m3 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_NO2'] = 0.05 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_NH4'] = 0.07 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)

comammox = {}
comammox['mu_max'] = 0.14 #d-1 (Kits et al. Nature 549, 269, 2017; 14.8 umolN/mgProt.h, 400 mgProt/molN)
comammox['b'] = 0.003 #d-1 (assumed)
comammox['K_O2'] = 0.4 #gO2/m3 (assumed)
comammox['K_NH4'] = 0.012 #gN/m3 (Kits et al. Nature 549, 269, 2017; 0.84 uM)
comammox['Ki_NH4'] = 3.44 #gN/m3 (Sakoula et al. ISME 15, 1010, 2020; 246 uM)

oho = {}
oho['mu_max_O2'] = 6 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['mu_max_NOx'] = 4.8 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['b'] = 0.62 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_s'] = 20 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_O2'] = 0.2 #gO2/m3 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_NOx'] = 0.3 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)

#Function to view growth rates of microbes included in the model
def view_rates(savename=None):
    points = 200
    xinput = np.tile(ifas.default_Xbulk, (points, 1))

    plt.rcParams.update({'font.size':10})
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(13/2.54, 22/2.54))

    #NH4
    xlist = np.linspace(0, 10, points)
    sinput = np.tile(ifas.default_Sbulk, (points, 1))
    sinput[:, 1] = xlist
    y = ifas.bio_rates(spar=sinput.T, xpar=xinput.T)
    ax[0].plot(xlist, y[0]-y[1], label='AOB')
    ax[0].plot(xlist, y[6]-y[7], label='Comammox')
    ax[0].set_xlabel('NH4-N (mg/L)')
    ax[0].legend(bbox_to_anchor=(1,1), loc=2)

    #NO2
    sinput = np.tile(ifas.default_Sbulk, (points, 1))
    sinput[:, 2] = xlist
    y = ifas.bio_rates(spar=sinput.T, xpar=xinput.T)
    ax[1].plot(xlist, y[2]-y[3], label='NOB')
    ax[1].set_xlabel('NO2-N (mg/L)')
    ax[1].legend(bbox_to_anchor=(1,1), loc=2)

    #O2
    xlist = np.linspace(0, 5, points)
    sinput = np.tile(ifas.default_Sbulk, (points, 1))
    sinput[:, 0] = xlist
    y = ifas.bio_rates(spar=sinput.T, xpar=xinput.T)
    ax[2].plot(xlist, y[0]-y[1], label='AOB')
    ax[2].plot(xlist, y[2]-y[3], label='NOB')
    ax[2].plot(xlist, y[4]-y[5], label='Anammox')
    ax[2].plot(xlist, y[6]-y[7], label='Comammox')
    ax[2].set_xlabel('O2 (mg/L)')
    ax[2].legend(bbox_to_anchor=(1,1), loc=2)

    ax[3].plot(xlist, y[8]-y[11], label='OHO on O2')
    ax[3].plot(xlist, (y[9]+y[10])/2-y[11], label='OHO on NOx')
    ax[3].set_xlabel('O2 (mg/L)')
    ax[3].legend(bbox_to_anchor=(1,1), loc=2)
    
    for i in range(4):
        ax[i].set_yticks([])
        ax[i].set_ylabel('Net growth rate')
    
    plt.tight_layout()
    if savename != None:
        plt.savefig(savename)
