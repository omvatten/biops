import numpy as np
import pandas as pd
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import pkgutil
import io
from . import yields
from . import kinetics

def load_default_influent():
    stream = pkgutil.get_data(__name__, 'data/default_influent.csv')
    return pd.read_csv(io.BytesIO(stream), encoding='utf8')

xparlist = ['X_I', 'X_AOB', 'X_NOB', 'X_AMX', 'X_CMX', 'X_OHO'] 
sparlist = ['S_O2', 'S_NH4', 'S_NO2', 'S_NO3', 'S_s']
default_Xbulk = np.array([1, 1, 1, 1, 1, 1], dtype=float)
default_Sbulk = np.array([4, 10, 10, 10, 20], dtype=float)
BF_layers = 20
Dw_NOx = 0.00017 #m2/d (Li et al. STOTEN 669, 683-691, 2019)
Dw_NH4 = 0.00019 #m2/d (Li et al. STOTEN 669, 683-691, 2019)
Dw_O2 = 0.00022 #m2/d (Li et al. STOTEN 669, 683-691, 2019)
Dw_S = 0.0001 #m2/d (Li et al. STOTEN 669, 683-691, 2019)
Dw = np.array([Dw_O2, Dw_NH4, Dw_NOx, Dw_NOx, Dw_S])
Dbf = Dw*0.8
#xparlist = [0'X_I', 1'X_AOB', 2'X_NOB', 3'X_AMX', 4'X_CMX', 5'X_OHO'] 
#sparlist = [0'S_O2', 1'S_NH4', 2'S_NO2', 3'S_NO3', 4'S_s', ]

#Function to calculate growth and decay rates
def bio_rates(xpar=default_Xbulk, spar=default_Sbulk, aob=kinetics.aob, nob=kinetics.nob, amx=kinetics.anammox, cmx=kinetics.comammox, oho=kinetics.oho):
    AOB_growth = aob['mu_max']*xpar[1]*(spar[1]/(aob['K_NH4']+spar[1]))*(spar[0]/(aob['K_O2']+spar[0]))
    AOB_decay = aob['b']*xpar[1]
    NOB_growth = nob['mu_max']*xpar[2]*(spar[2]/(nob['K_NO2']+spar[2]))*(spar[0]/(nob['K_O2']+spar[0]))
    NOB_decay = nob['b']*xpar[2]
    AMX_growth = amx['mu_max']*xpar[3]*(spar[1]/(amx['K_NH4']+spar[1]))*(spar[2]/(amx['K_NO2']+spar[2]))*(amx['K_O2']/(amx['K_O2']+spar[0]))
    AMX_decay = amx['b']*xpar[3]
    CMX_growth = cmx['mu_max']*xpar[4]*(spar[1]/(cmx['K_NH4']+spar[1]+spar[1]**2/cmx['Ki_NH4']))*(spar[0]/(cmx['K_O2']+spar[0]))
    CMX_decay = cmx['b']*xpar[4]
    OHO_growth_O2 = oho['mu_max_O2']*xpar[5]*(spar[4]/(oho['K_s']+spar[4]))*(spar[0]/(oho['K_O2']+spar[0]))
    OHO_growth_NO2 = oho['mu_max_NOx']*xpar[5]*(spar[4]/(oho['K_s']+spar[4]))*(spar[2]/(oho['K_NOx']+spar[2]))*(oho['K_O2']/(oho['K_O2']+spar[0]))
    OHO_growth_NO3 = oho['mu_max_NOx']*xpar[5]*(spar[4]/(oho['K_s']+spar[4]))*(spar[3]/(oho['K_NOx']+spar[3]))*(oho['K_O2']/(oho['K_O2']+spar[0]))
    OHO_decay = oho['b']*xpar[5]
    return np.array([AOB_growth, AOB_decay, NOB_growth, NOB_decay, AMX_growth, AMX_decay, CMX_growth, CMX_decay, OHO_growth_O2, OHO_growth_NO2, OHO_growth_NO3, OHO_decay])
    #[0 AOB_growth, 1 AOB_decay, 2 NOB_growth, 3 NOB_decay, 4 AMX_growth, 5 AMX_decay, 6 CMX_growth, 7 CMX_decay, 8 OHO_growth_O2, 9 OHO_growth_NO2, 10 OHO_growth_NO3, 11 OHO_decay])

class Reactor:
    def __init__(self, influent='default', volume=4000, VSS=3000, DO=2, Xbulk=default_Xbulk, Sbulk=default_Sbulk, area=100, Lmax=400*(10**-6), Lbf=5*(10**-6), VSmax=25000, DBL=10*(10**-6)):
        if influent == 'default':
            inf = load_default_influent()
            for p in xparlist+sparlist:
                if p not in inf.columns:
                    inf[p] = 0
            self.influents = [inf]
        elif isinstance(influent, pd.DataFrame):
            for p in xparlist+sparlist:
                if p not in influent.columns:
                    influent[p] = 0
            self.influents = [influent]
            if 'Q' not in influent.columns:
                print('Warning, Q data is not in influent')
            if 'Time' not in influent.columns:
                print('Warning, Time data is not in influent')
        elif isinstance(influent, dict) or isinstance(influent, pd.Series):
            influent = pd.Series(influent)
            for p in xparlist+sparlist:
                if p not in influent.index:
                    influent[p] = 0
            influent = influent.to_frame().transpose()
            self.influents = [influent]
            if 'Q' not in influent.columns:
                print('Warning, Q is not in influent')
        self.volumes = [volume]
        self.VSS = [VSS]
        self.DO = [DO]
        Sbulk[0] = DO
        self.Xbulk = [VSS*Xbulk/(Xbulk.sum())]
        self.Sbulk = [Sbulk]
        self.areas = [area]
        self.Lmax = [Lmax]
        self.Lbf = [Lbf]
        self.VSmax = [VSmax]
        self.DBL = [DBL]
        bf_x = VSmax*Xbulk/Xbulk.sum()
        bf_x = np.tile(bf_x, (BF_layers+1, 1))
        bf_x[-1] = np.zeros(len(bf_x[-1]))
        bf_s = np.tile(Sbulk, (BF_layers+1, 1))
        self.Xbiofilm = [bf_x]
        self.Sbiofilm = [bf_s]
        self.recircT = {}
        self.recircO = {}

        self.time = 0
        self.bulklogger = {}
        self.biofilmlogger = {}
        self.bulklogger['1'] = {'Conc':pd.DataFrame(np.concatenate([self.Xbulk[-1],self.Sbulk[-1]]).reshape(1,-1), index=[self.time], columns=xparlist+sparlist)}
        self.bulklogger['1']['Prod'] = pd.DataFrame(0, index=[round(self.time,2)], columns=xparlist+sparlist)
        if area > 0:
            self.biofilmlogger['1'] = {'Conc':{round(self.time,2):pd.DataFrame(np.concatenate([self.Xbiofilm[-1], self.Sbiofilm[-1]], axis=1), index=range(BF_layers+1), columns=xparlist+sparlist)}}
            self.biofilmlogger['1']['Conc'][round(self.time,2)]['thickness'] = np.linspace(0, self.Lbf[-1], BF_layers+1)
            self.biofilmlogger['1']['Prod'] = pd.DataFrame(0, index=[round(self.time,2)], columns=xparlist+sparlist)
        self.influentlogger = pd.DataFrame(self.influents[-1][['Q']+xparlist+sparlist].to_numpy()[0,:].reshape(1,-1), index=[round(self.time,2)], columns=['Q']+xparlist+sparlist)

    def add_compartment(self, influent=None, volume='same', DO='same', Xbulk='same', Sbulk='same', area='same', Lmax='same', Lbf='same', VSmax='same', DBL='same'):
        if isinstance(influent, pd.DataFrame):
            for p in xparlist+sparlist:
                if p not in influent.columns:
                    influent[p] = 0
            self.influents.append(influent)
            if 'Q' not in influent.columns:
                print('Warning, Q is not in influent')
        elif isinstance(influent, dict) or isinstance(influent, pd.Series):
            influent = pd.Series(influent)
            for p in xparlist+sparlist:
                if p not in influent.index:
                    influent[p] = 0
            influent = influent.to_frame().transpose()
            self.influents.append(influent)
            if 'Q' not in influent.columns:
                print('Warning, Q is not in influent')
        else:
            self.influents.append(influent)

        self.VSS.append(self.VSS[-1])
        if volume == 'same':
            self.volumes.append(self.volumes[-1])
        else:
            self.volumes.append(volume)
        if DO == 'same':
            self.DO.append(self.DO[-1])
        else:
            self.DO.append(DO)
        if Xbulk == 'same':
            self.Xbulk.append(self.Xbulk[-1])
        else:
            Xbulk = self.VSS[-1]*Xbulk/Xbulk.sum()
            self.Xbulk.append(Xbulk)
        if Sbulk == 'same':
            self.Sbulk.append(self.Sbulk[-1])
            self.Sbulk[-1][0] = self.DO[-1]
        else:
            Sbulk[0] = self.DO[-1]
            self.Sbulk.append(Sbulk)

        if area == 'same':
            self.areas.append(self.areas[-1])
        else:
            self.areas.append(area)
        if Lmax == 'same':
            self.Lmax.append(self.Lmax[-1])
        else:
            self.Lmax.append(Lmax)
        if Lbf == 'same':
            self.Lbf.append(self.Lbf[-1])
        else:
            self.Lbf.append(Lbf)
        if VSmax == 'same':
            self.VSmax.append(self.VSmax[-1])
        else:
            self.VSmax.append(VSmax)

        if DBL == 'same':
            self.DBL.append(self.DBL[-1])
        else:
            self.DBL.append(DBL)
        bf_x = self.VSmax[-1]*self.Xbulk[-1]/self.Xbulk[-1].sum()
        bf_x = np.tile(bf_x, (BF_layers+1, 1))
        bf_x[-1] = np.zeros(len(bf_x[-1]))
        self.Xbiofilm.append(bf_x)
        bf_s = np.tile(self.Sbulk[-1], (BF_layers+1, 1))
        self.Sbiofilm.append(bf_s)

        cnr = str(len(self.volumes))
        self.bulklogger[cnr] = {'Conc':pd.DataFrame(np.concatenate([self.Xbulk[-1],self.Sbulk[-1]]).reshape(1,-1), index=[self.time], columns=xparlist+sparlist)}
        self.bulklogger[cnr]['Prod'] = pd.DataFrame(0, index=[round(self.time,2)], columns=xparlist+sparlist)
        if area > 0:
            self.biofilmlogger[cnr] = {'Conc':{round(self.time,2):pd.DataFrame(np.concatenate([self.Xbiofilm[-1], self.Sbiofilm[-1]], axis=1), index=range(BF_layers+1), columns=xparlist+sparlist)}}
            self.biofilmlogger[cnr]['Conc'][round(self.time,2)]['thickness'] = np.linspace(0, self.Lbf[-1], BF_layers+1)
            self.biofilmlogger[cnr]['Prod'] = pd.DataFrame(0, index=[round(self.time,2)], columns=xparlist+sparlist)

    def add_recirculation(self, origin=None, target=None, Q=None):
        self.recircT[target] = [origin, Q]
        self.recircO[origin] = [target, Q]

    def calculate(self, timestep=1, dt=0.1, set_iterations=1000, set_dts='auto'):
        #Process matrices for S and X
        Xpm = yields.get_matrix()[xparlist].values
        Spm = yields.get_matrix()[sparlist].values

        #Check influent dataframe and biofilm calculation parameters
        influent_Vals = {}
        influent_Time = {}
        header = ['Q']+xparlist+sparlist
        BF_Xprod = {}
        BF_Sprod = {}
        Bulk_Xprod = {}
        Bulk_Sprod = {}
        influent_log = np.zeros(len(xparlist)+len(sparlist)+1)

        for i in range(len(self.volumes)):
            if isinstance(self.influents[i], pd.DataFrame) and 'Time' in self.influents[i].columns and len(self.influents[i].index) > 1:
                influent_Vals[i] = self.influents[i][header].values
                influent_Time[i] = self.influents[i]['Time'].values
            if isinstance(self.influents[i], pd.DataFrame) and len(self.influents[i].index) == 1:
                influent_Vals[i] = self.influents[i][header].values[0]
            Bulk_Xprod[i] = np.zeros(len(xparlist))
            Bulk_Sprod[i] = np.zeros(len(sparlist))
            BF_Xprod[i] = np.zeros(len(xparlist))
            BF_Sprod[i] = np.zeros(len(sparlist))

        #Set substrate calculation step
        if set_dts == 'auto':
            testS = np.array([20]*len(sparlist))
            testX = np.array([1]*len(xparlist))
            dts = 0.001
            for i in range(len(self.volumes)):
                if self.areas[i] > 0:
                    testXbf = self.VSmax[i]*testX/sum(testX)
                    br = bio_rates(xpar=testXbf, spar=testS).reshape((-1,1)) #g/m3-d
                    brS = br*Spm
                    brS = brS.sum(axis=0)
                    if min(abs(0.1/brS[brS>0])) < dts:
                        dts = min(abs(0.1/brS[brS!=0]))
                if self.VSS[i] > 0:
                    testXbulk = self.VSS[i]*testX/sum(testX)
                    br = bio_rates(xpar=testXbulk, spar=testS).reshape((-1,1)) #g/m3-d
                    brS = br*Spm
                    brS = brS.sum(axis=0)
                    if min(abs(0.1/brS[brS>0])) < dts:
                        dts = min(abs(0.1/brS[brS>0]))
        else:
            dts = set_dts

        #### Go through the timesteps ####
        VSS = self.VSS[0]
        for tx in range(round(timestep/dt)): #This is timestep for X
            self.time = self.time + dt

            ### Get influent and biofilm characteristics at this time for each volume
            Qinf = {}
            A_matrix_bf = {}
            Sbulk_prev = {}
            bulk_br_prev = {}
            Sbf_prev = {}
            bf_br_prev = {}

            for i in range(len(self.volumes)):
                #Get influent
                if i in influent_Vals.keys():
                    if i in influent_Time.keys():
                        ix = self.time % influent_Time[i][-1]
                        ix = int(len(influent_Time[i])*ix/influent_Time[i][-1])
                        Qinf[i] = influent_Vals[i][ix]
                    else:
                        Qinf[i] = influent_Vals[i]
                    if i == 0:
                        influent_log = influent_log + dt*Qinf[0]
                
                #Set variables for bulk
                Sbulk_prev[i] = self.Sbulk[i]
                bulk_br_prev[i] = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1))

                #Set variables for biofilm calculations
                if self.areas[i] > 0:
                    Sbf_prev[i] = self.Sbiofilm[i].copy()
                    bf_br_prev[i] = bio_rates(xpar=self.Xbiofilm[i][:-1].T, spar=self.Sbiofilm[i][:-1].T)

                    dx = self.Lbf[i]/(BF_layers)
                    A_matrix_bf[i] = []
                    for s in range(len(sparlist)):
                        F_dx = 2*dts*Dbf[s]/(dx**2)
                        F_BL = 2*dts*Dw[s]/(self.DBL[i]**2)
                        main = np.zeros(BF_layers+1)
                        lower = np.zeros(BF_layers)
                        upper = np.zeros(BF_layers)
                        main[:] = 3 + 2*F_dx
                        lower[:] = -F_dx
                        upper[:] = -F_dx
                        main[0] = 3 + F_dx
                        main[-1] = 3
                        main[-2] = 3 + (F_dx+F_BL)
                        lower[-1] = 0
                        upper[-1] = -F_BL
                        A_matrix_bf[i].append(scipy.sparse.diags(diagonals=[main, lower, upper], offsets=[0, -1, 1], shape=(BF_layers+1, BF_layers+1), format='csr'))
                        #print(A_matrix_bf[i][-1].todense())

            ### Update substrate concentrations ###
            run_bulk_substrate_update = True
            run_biofilm_substrate_update = True
            run_substrate_time = 0
            run_substrate_count = 0
            while ((run_bulk_substrate_update or run_biofilm_substrate_update) and (run_substrate_time < dt and run_substrate_count <= set_iterations)) or (run_substrate_count <= 10):
                run_substrate_time = run_substrate_time + dts
                run_substrate_count += 1
                #if run_substrate_count == max_iterations:
                #    print('Warning, reached max iterations of substrate field')

                Qin_total = {}
                Qeff_balance = {}

                for i in range(len(self.volumes)):
                    delta_s = np.zeros(len(sparlist))

                    #Update biofilm concentration profiles
                    if self.areas[i] > 0:
                        new_Sbiofilm = self.Sbiofilm[i].copy()
                        bf_br = bio_rates(xpar=self.Xbiofilm[i][:-1].T, spar=self.Sbiofilm[i][:-1].T)
                        Sbiofilm_change = 0 #Variable to check change in conc after each iteration
                        for s in range(len(sparlist)):
                            bf_brS = bf_br*Spm[:, s].reshape(-1,1)
                            bf_brS = bf_brS.sum(axis=0) #g/m3.d change
                            bf_brS1 = bf_br_prev[i]*Spm[:, s].reshape(-1,1)
                            bf_brS1 = bf_brS1.sum(axis=0) #g/m3.d change

                            b = self.Sbiofilm[i][:, s].copy()
                            b = 4*b - Sbf_prev[i][:, s]
                            b[:-1] = b[:-1] + 2*dts*(2*bf_brS-bf_brS1)
 
                            new_Sbiofilm[:, s] = scipy.sparse.linalg.spsolve(A_matrix_bf[i][s], b)
                            new_Sbiofilm[-1, s] = self.Sbulk[i][s]
                            new_Sbiofilm[:, s][new_Sbiofilm[:, s]<0] = 0

                            change_conc_profile = max(abs(new_Sbiofilm[:-1, s]-self.Sbiofilm[i][:-1, s]))
                            if change_conc_profile > Sbiofilm_change:
                                Sbiofilm_change = change_conc_profile

                            #Calculate biofilm contribution to bulk removal of s
                            delta_s[s] = sum(bf_brS)*self.Lbf[i]*self.areas[i]/self.volumes[i]

                        Sbf_prev[i] = self.Sbiofilm[i].copy()
                        self.Sbiofilm[i] = new_Sbiofilm.copy()
                        bf_br_prev[i] = bf_br.copy()

                        #Check if biofilm conc profile is stable
                        max_Sbiofilm_change = (1/(24*60))*(Sbiofilm_change/dts)
                        if max_Sbiofilm_change < 0.1:
                            run_biofilm_substrate_update = False
                        #if run_substrate_count == max_iterations:
                        #    print('Biofilm conc change: ', round(max_Sbiofilm_change,2), 'mg/L.min')
                
                    ## Update bulk concentrations
                    Qin_total[i] = 0
                    Qeff_balance[i] = 0

                    #Get influent
                    if i in influent_Vals.keys():
                        Qin_total[i] = Qin_total[i] + Qinf[i][0]
                        delta_s = delta_s + Qinf[i][0]*Qinf[i][len(xparlist)+1:]/self.volumes[i]
    
                    #Get return sludge flow flow (assumed 0.5xQin)
                    if i == 0 and VSS > 0:
                        Qin_total[i] = Qin_total[i] + 0.5*Qinf[i][0]
                        delta_s = delta_s + 0.5*Qinf[i][0]*self.Sbulk[-1]/self.volumes[i] #m3/d*d*g/m3 = g
    
                    #Get recirculation flow
                    if i+1 in self.recircT.keys(): #Flow coming into i
                        Q = self.recircT[i+1][1]
                        Qin_total[i] = Qin_total[i] + Q
                        origin = self.recircT[i+1][0]-1
                        delta_s = delta_s + Q*self.Sbulk[origin]/self.volumes[i]
                    if i+1 in self.recircO.keys(): #Flow going out from i
                        Qeff_balance[i] = Qeff_balance[i] - self.recircO[i+1][1]

                    #Get flow from previous
                    if i > 0:
                        Q = Qeff_balance[i-1]
                        Qin_total[i] = Qin_total[i] + Q
                        delta_s = delta_s + Q*self.Sbulk[i-1]/self.volumes[i]

                    #Calculate Q going out and into the next compartment
                    Qeff_balance[i] = Qeff_balance[i] + Qin_total[i]

                    #Get biorates bulk
                    if VSS > 0:
                        br = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1)) #g/m3-d
                        brS = br*Spm
                        brS = brS.sum(axis=0)
                        brS1 = bulk_br_prev[i]*Spm
                        brS1 = brS1.sum(axis=0)

                    #Get new concentrations in current volume
                    new_Sbulk = (4*self.Sbulk[i]-Sbulk_prev[i]+2*dts*(2*brS-brS1)+2*dts*delta_s)/(3+2*Qin_total[i]*dts/self.volumes[i])
                    new_Sbulk[0] = self.DO[i]
                    new_Sbulk[new_Sbulk<0] = 0

                    #Check if new bulk concentration is stable
                    max_Sbulk_change = (1/(24*60))*(max(abs(new_Sbulk-self.Sbulk[i])/dts))
                    if max_Sbulk_change < 0.1:
                        run_bulk_substrate_update = False
                    #if run_substrate_count == max_iterations:
                    #    print('Bulk conc change: ', round(max_Sbulk_change,2), 'mg/L.min')
                    Sbulk_prev[i] = self.Sbulk[i].copy()
                    self.Sbulk[i] = new_Sbulk
            ###End of while loop for updating substrate conc. ###

            ### Update biofilm biomass concentrations
            for i in range(len(self.volumes)):
                if self.areas[i] > 0:
                    bf_br = bio_rates(xpar=self.Xbiofilm[i][:-1].T, spar=self.Sbiofilm[i][:-1].T)
                    for x in range(len(xparlist)):
                        bf_brX = bf_br*Xpm[:, x].reshape(-1,1)
                        bf_brX = bf_brX.sum(axis=0)
                        self.Xbiofilm[i][:-1, x] = self.Xbiofilm[i][:-1, x] + dt*bf_brX
                    #Calculate new biofilm thickness
                    new_Lbf = self.Lbf[i]*self.Xbiofilm[i][:-1].sum(axis=1).mean()/self.VSmax[i]
                    if new_Lbf > self.Lmax[i]:
                        new_Lbf = self.Lmax[i]
                    self.Xbiofilm[i] = (self.Lbf[i]/new_Lbf)*self.Xbiofilm[i]
                    for lay in range(BF_layers):
                        diffX = self.Xbiofilm[i][lay].sum() - self.VSmax[i]
                        if lay == BF_layers-1 and diffX > 0:
                            self.Xbiofilm[i][lay] = self.Xbiofilm[i][lay]/self.Xbiofilm[i][lay].sum()
                            self.Xbulk[i] = self.Xbulk[i] + diffX*(self.Xbiofilm[i][lay]/self.Xbiofilm[i][lay].sum())*self.areas[i]*(new_Lbf/BF_layers)/self.volumes[i]
                        elif lay < BF_layers-1 and diffX < 0:
                            addX = abs(diffX)*self.Xbiofilm[i][lay+1]/self.Xbiofilm[i][lay+1].sum()
                            self.Xbiofilm[i][lay] = self.Xbiofilm[i][lay] + addX
                            self.Xbiofilm[i][lay+1] = self.Xbiofilm[i][lay+1] - addX
                        elif diffX > 0:
                            subX = diffX*self.Xbiofilm[i][lay]/self.Xbiofilm[i][lay].sum()
                            self.Xbiofilm[i][lay] = self.Xbiofilm[i][lay] - subX
                            self.Xbiofilm[i][lay+1] = self.Xbiofilm[i][lay+1] + subX
                    self.Lbf[i] = new_Lbf

            ### Update bulk biomass concentrations (assume biomass is mixed in whole system)
            amt_x = np.zeros(len(xparlist)) #This will be in g
            total_volume = 0
            for i in range(len(self.volumes)):
                #Get amt in volume
                total_volume = total_volume + self.volumes[i]
                amt_x = amt_x + self.Xbulk[i]*self.volumes[i]
                #Get influent
                if i in influent_Vals.keys():
                    amt_x = amt_x + dt*Qinf[i][0]*Qinf[i][1:len(xparlist)+1]
                #Get biorates bulk
                if VSS > 0:
                    br = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1)) #g/m3-d
                    brX = br*Xpm
                    brX = brX.sum(axis=0)
                    amt_x = amt_x + dt*brX*self.volumes[i]
            #Calculate new Xbulk composition and put in volumes
            new_Xbulk = amt_x/total_volume
            if new_Xbulk.sum() > VSS:
                new_Xbulk = VSS*new_Xbulk/new_Xbulk.sum()
            for i in range(len(self.volumes)):
                self.Xbulk[i] = new_Xbulk

            ###Calculate production and consumption of S and X for updating loggers
            for i in range(len(self.volumes)):
                if VSS > 0:
                    br = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1)) #g/m3-d
                    brX = br*Xpm
                    brS = br*Spm
                    Bulk_Xprod[i] = Bulk_Xprod[i] + brX.sum(axis=0)*dt
                    Bulk_Sprod[i] = Bulk_Sprod[i] + brS.sum(axis=0)*dt

                if self.areas[i] > 0:
                    bf_br = bio_rates(xpar=self.Xbiofilm[i].T, spar=self.Sbiofilm[i].T)
                    for s in range(len(sparlist)):
                        bf_brS = bf_br*Spm[:, s].reshape(-1,1)
                        bf_brS = bf_brS.sum(axis=0) #g/m3.d change
                        bf_brS = bf_brS*(self.Lbf[i]/BF_layers) #g/m2.d per layer
                        BF_Sprod[i][s] = BF_Sprod[i][s] + dt*bf_brS.sum()*self.areas[i]/self.volumes[i] #g/m3
                    for x in range(len(xparlist)):
                        bf_brX = bf_br*Xpm[:, s].reshape(-1,1)
                        bf_brX = bf_brX.sum(axis=0) #g/m3.d change
                        bf_brX = bf_brX*(self.Lbf[i]/BF_layers) #g/m2.d per layer
                        BF_Xprod[i][s] = BF_Xprod[i][s] + dt*bf_brX.sum()*self.areas[i]/self.volumes[i] #g/m3
        #### End of dt loop

        #### Update loggers in the end of calculation ####
        for i in range(len(self.volumes)):
            self.bulklogger[str(i+1)]['Conc'].loc[self.time, xparlist+sparlist] = np.concatenate([self.Xbulk[i], self.Sbulk[i]])
            if VSS > 0:
                self.bulklogger[str(i+1)]['Prod'].loc[self.time, xparlist+sparlist] = np.concatenate([Bulk_Xprod[i]/timestep, Bulk_Sprod[i]/timestep]) #Get prod in g/m3.d
            if self.areas[i] > 0:
                self.biofilmlogger[str(i+1)]['Conc'][round(self.time,2)] = pd.DataFrame(np.concatenate([self.Xbiofilm[i], self.Sbiofilm[i]], axis=1), index=range(BF_layers+1), columns=xparlist+sparlist)
                self.biofilmlogger[str(i+1)]['Conc'][round(self.time,2)]['thickness'] = np.linspace(0, self.Lbf[i], BF_layers+1)
                self.biofilmlogger[str(i+1)]['Prod'].loc[round(self.time,2), xparlist+sparlist] = np.concatenate([BF_Xprod[i]/timestep, BF_Sprod[i]/timestep]) #g/m3.d, directly comparable to bulk
        self.influentlogger.loc[round(self.time-timestep/2, 2)] = influent_log/timestep

    # Function to plot schematic of process setup
    def view_process(self, savename=None):
        plt.rcParams.update({'font.size':10})
        fig, ax = plt.subplots(figsize=(20/2.54, 10/2.54))
        
        #Make large rectangle and settler
        rect = [Rectangle((10, 10), width=100, height=20)]
        ax.add_collection(PatchCollection(rect, facecolor='white', edgecolor='black'))
        ax.plot([120, 130, 140], [30, 10, 30], color='black', lw=1)
        ax.text(130, 26, 'Settler', horizontalalignment='center', verticalalignment='center', backgroundcolor='white')
        ax.arrow(130, 10, 0, -10, head_width=1, head_length=1, length_includes_head=True, ec='grey', fc='grey')
        if sum(self.VSS) > 0:
            ax.plot([130, 10], [5, 5], color='grey', lw=1)
            ax.arrow(10, 5, 0, 5, head_width=1, head_length=1, length_includes_head=True, ec='grey', fc='grey')

        ax.arrow(0, 20, 10, 0, head_width=1, head_length=1, length_includes_head=True, ec='black', fc='black')
        ax.arrow(110, 20, 10, 0, head_width=1, head_length=1, length_includes_head=True, ec='black', fc='black')
        ax.arrow(140, 20, 10, 0, head_width=1, head_length=1, length_includes_head=True, ec='black', fc='black')
        ax.set_xlim(0,150)
        ax.set_ylim(0,40)

        total_vol = sum(self.volumes)
        xlist = [10]
        for i in range(len(self.volumes)):
            xlist.append(xlist[-1]+100*self.volumes[i]/total_vol)
            if int(xlist[-1]) != 110:
                ax.plot([xlist[-1], xlist[-1]], [10, 30], color='black', ls='--')
            text2add = 'Comp.'+str(i+1)+'\nV='+str(self.volumes[i])+'\nVSS='+str(self.VSS[i])+'\narea='+str(self.areas[i])+'\nDO='+str(self.DO[i])
            ax.text((xlist[-1]+xlist[-2])/2, 20, text2add, horizontalalignment='center', verticalalignment='center', backgroundcolor='white')
            if (isinstance(self.influents[i], pd.DataFrame) or isinstance(self.influents[i], dict)) and i > 0:
                ax.arrow((xlist[-1]+xlist[-2])/2, 37, 0, -7, head_width=1, head_length=1, length_includes_head=True, ec='black', fc='black')

        for targ, [orig, flow] in self.recircT.items():
            xkeyT = (xlist[targ] + xlist[targ-1]) / 2
            xkeyO = (xlist[orig] + xlist[orig-1]) / 2
            ax.plot([xkeyT, xkeyO, xkeyO], [35, 35, 30], color='grey', lw=1)
            ax.arrow(xkeyT, 35, 0, -5, head_width=1, head_length=1, length_includes_head=True, ec='grey', fc='grey')
            
        ax.axis('off')
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()
        if savename != None:
            plt.savefig(savename)

    def view_influent(self, savename=None):
        plt.rcParams.update({'font.size':10})
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20/2.54, 16/2.54))

        for i in range(len(self.volumes)):
            if isinstance(self.influents[i], pd.DataFrame) and 'Time' in self.influents[i].columns and len(self.influents[i].index) > 1:
                ax[0, 0].plot(self.influents[i]['Time'], self.influents[i]['Q'], label=str(i+1))
                ax[0, 1].plot(self.influents[i]['Time'], self.influents[i]['S_s'], label=str(i+1))
                ax[1, 0].plot(self.influents[i]['Time'], self.influents[i]['S_NH4'], label=str(i+1))
                ax[1, 1].plot(self.influents[i]['Time'], self.influents[i]['S_NO2'], label=str(i+1)+':NO2')
                ax[1, 1].plot(self.influents[i]['Time'], self.influents[i]['S_NO3'], label=str(i+1)+':NO3')

            if isinstance(self.influents[i], pd.DataFrame) and len(self.influents[i].index) == 1:
                xmin, xmax = ax[0, 0].get_xlim()
                ax[0, 0].plot([xmin, xmax], [self.influents[i].loc[0, 'Q']]*2, label=str(i+1))
                xmin, xmax = ax[0, 1].get_xlim()
                ax[0, 1].plot([xmin, xmax], [self.influents[i].loc[0, 'S_s']]*2, label=str(i+1))
                xmin, xmax = ax[1, 0].get_xlim()
                ax[1, 0].plot([xmin, xmax], [self.influents[i].loc[0, 'S_NH4']]*2, label=str(i+1))
                xmin, xmax = ax[1, 1].get_xlim()
                ax[1, 1].plot([xmin, xmax], [self.influents[i].loc[0, 'S_NO2']]*2, label=str(i+1)+':NO2')
                ax[1, 1].plot([xmin, xmax], [self.influents[i].loc[0, 'S_NO3']]*2, label=str(i+1)+':NO3')

        ax[0, 0].set_ylabel('Q (m3/d)')
        ax[0, 1].set_ylabel('S_s (g/m3)')
        ax[1, 0].set_ylabel('S_NH4 (g/m3)')
        ax[1, 1].set_ylabel('S_NO2, S_NO3 (g/m3)')

        for r in range(2):
            for c in range(2):
                ax[r, c].legend()
                if r == 1:
                    ax[r, c].set_xlabel('Time (d)')

        plt.tight_layout()
        if savename != None:
            plt.savefig(savename)
