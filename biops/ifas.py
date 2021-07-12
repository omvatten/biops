import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import yields
import kinetics

xparlist = ['X_I', 'X_AOB', 'X_NOB', 'X_AMX', 'X_CMX', 'X_OHO'] 
sparlist = ['S_O2', 'S_NH4', 'S_NO2', 'S_NO3', 'S_s']
default_influent = pd.read_csv('default_influent.csv')
default_Xbulk = np.array([1, 1, 1, 1, 1, 1], dtype=float)
default_Sbulk = np.array([4, 10, 10, 10, 20], dtype=float)
default_Xbiofilm = np.tile(default_Xbulk, (1, 1))
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
    def __init__(self, influent=default_influent, volume=4000, VSS=3000, DO=2, Xbulk=default_Xbulk, Sbulk=default_Sbulk, area=100, Lmax=400, VSmax=25000, DBL=10, BFL=5, Xbiofilm=default_Xbiofilm):
        if isinstance(influent, pd.DataFrame):
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
        Xbulk = VSS*Xbulk/(Xbulk.sum())
        self.Xbulk = [Xbulk]
        self.Sbulk = [Sbulk]

        self.areas = [area]
        self.Lmax = [Lmax]
        self.VSmax = [VSmax]
        self.DBL = [DBL]
        self.BFL = [BFL]
        bf_layers = int(Lmax/BFL)
        bf_x = np.zeros((bf_layers, len(xparlist)))
        bf_x[bf_layers-len(Xbiofilm[:, 0]):, :] = Xbiofilm
        bf_s = np.tile(Sbulk, (bf_layers, 1))
        self.Xbiofilm = [bf_x]
        self.Sbiofilm = [bf_s]
        self.recircT = {}
        self.recircO = {}

        self.time = 0
        self.bulklogger = {}
        self.biofilmlogger = {}
        self.bulklogger['1'] = {'Conc':pd.DataFrame(np.concatenate([self.Xbulk[-1],self.Sbulk[-1]]).reshape(1,-1), index=[self.time], columns=xparlist+sparlist)}
        self.bulklogger['1']['Prod'] = pd.DataFrame(0, index=[self.time], columns=xparlist+sparlist)
        if area > 0:
            self.biofilmlogger['1'] = {'Conc':{self.time:pd.DataFrame(np.concatenate([self.Xbiofilm[-1], self.Sbiofilm[-1]], axis=1), index=np.arange(0, self.Lmax[-1], self.BFL[-1])+self.BFL[-1]/2, columns=xparlist+sparlist)}}
            self.biofilmlogger['1']['Prod'] = pd.DataFrame(0, index=[self.time], columns=xparlist+sparlist)
        
    def add_compartment(self, influent=None, volume='same', DO='same', Xbulk='same', Sbulk='same', area='same', Lmax='same', L='same', VSmax='same', DBL='same', BFL='same', Xbiofilm='same'):
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
        if VSmax == 'same':
            self.VSmax.append(self.VSmax[-1])
        else:
            self.VSmax.append(VSmax)

        if DBL == 'same':
            self.DBL.append(self.DBL[-1])
        else:
            self.DBL.append(DBL)
        if BFL == 'same':
            self.BFL.append(self.BFL[-1])
        else:
            self.BFL.append(BFL)
        bf_layers = int(self.Lmax[-1]/self.BFL[-1])
        bf_x = np.zeros((bf_layers, len(xparlist)))
        if Xbiofilm == 'same':
            bf_x[bf_layers-len(self.Xbiofilm[-1][:, 0]):, :] = self.Xbiofilm[-1]
            self.Xbiofilm.append(bf_x)
        else:
            bf_x[bf_layers-len(Xbiofilm[:, 0]):, :] = Xbiofilm
            self.Xbiofilm.append(bf_x)
        bf_s = np.tile(self.Sbulk[-1], (bf_layers, 1))
        self.Sbiofilm.append(bf_s)

        cnr = str(len(self.volumes))
        self.bulklogger[cnr] = {'Conc':pd.DataFrame(np.concatenate([self.Xbulk[-1],self.Sbulk[-1]])[0], index=[self.time], columns=xparlist+sparlist)}
        self.bulklogger[cnr]['Prod'] = pd.DataFrame(0, index=[self.time], columns=xparlist+sparlist)
        if area > 0:
            self.biofilmlogger[cnr] = {'Conc':{self.time:pd.DataFrame(np.concatenate([self.Xbiofilm[-1], self.Sbiofilm[-1]], axis=1), index=np.arange(0, self.Lmax[-1], self.BFL[-1])+self.BFL[-1]/2, columns=xparlist+sparlist)}}
            self.biofilmlogger[cnr]['Prod'] = pd.DataFrame(0, index=[self.time], columns=xparlist+sparlist)

    def add_recirculation(self, origin=None, target=None, Q=None):
        self.recircT[target] = [origin, Q]
        self.recircO[origin] = [target, Q]

    def calculate(self, timestep=1, dtx=0.1, dts=0.001):
        #Process matrices for S and X
        Xpm = yields.get_matrix()[xparlist].values
        Spm = yields.get_matrix()[sparlist].values

        #Check influent dataframe and biofilm calculation parameters
        influent_Vals = {}
        influent_Time = {}
        header = ['Q']+xparlist+sparlist
        dtbf = {}
        BF_Xprod = {}
        BF_Sprod = {}
        Bulk_Xprod = {}
        Bulk_Sprod = {}
        
        for i in range(len(self.influents)):
            if isinstance(self.influents[i], pd.DataFrame) and 'Time' in self.influents[i].columns and len(self.influents[i].index) > 1:
                influent_Vals[i] = self.influents[i][header].values
                influent_Time[i] = self.influents[i]['Time'].values
            if isinstance(self.influents[i], pd.DataFrame) and len(self.influents[i].index) == 1:
                influent_Vals[i] = self.influents[i][header].values[0]
            dtbf[i] = 0.4*((self.BFL[i]*10**-6)**2)/Dbf.max()
            Bulk_Xprod[i] = np.zeros(len(xparlist))
            Bulk_Sprod[i] = np.zeros(len(sparlist))
            BF_Xprod[i] = np.zeros(len(xparlist))
            BF_Sprod[i] = np.zeros(len(sparlist))
        bf_tsteps = min(int(dts/dtbf[i]), 500)

        ## Go through the timesteps
        VSS = self.VSS[0]
        for tx in range(round(timestep/dtx)): #This is timestep for X
            self.time = self.time + dtx
            totalQin = {}
            BF_len = {}
            BF_s = {}
            BF_x = {}

            for i in range(len(self.volumes)):
                #Get influent
                if i in influent_Vals.keys():
                    if i in influent_Time.keys():
                        ix = self.time % influent_Time[i][-1]
                        ix = int(len(influent_Time[i])*ix/influent_Time[i][-1])
                        inf = influent_Vals[i][ix]
                    else:
                        inf = influent_Vals[i]

                #Get biofilm length and subset to part with biomass
                BF_len[i] = len(self.Xbiofilm[i][:, 0][self.Xbiofilm[i][:, 0] != 0])
                BF_s[i] = self.Sbiofilm[i][len(self.Sbiofilm[i][:,0])-BF_len[i]:]
                BF_x[i] = self.Xbiofilm[i][len(self.Xbiofilm[i][:,0])-BF_len[i]:]

            # Go through dts
            for ts in range(int(dtx/dts)):
                delta_x = np.zeros(len(xparlist))
                delta_s = np.zeros(len(sparlist))

                for i in range(len(self.volumes)):
                    totalQin[i] = 0
                    #Get influent
                    if i in influent_Vals.keys():
                        delta_x = delta_x + inf[0]*dts*inf[1:len(xparlist)+1] #m3/d*d*g/m3 = g
                        delta_s = delta_s + inf[0]*dts*inf[len(xparlist)+1:] #m3/d*d*g/m3 = g
                        totalQin[i] = totalQin[i] + inf[0]
    
                    #Get return sludge flow flow
                    if i == 0 and VSS != 0:
                        totalQin[i] = totalQin[i] + 0.5*inf[0]
                        returned_flowX = self.Xbulk[-1].copy()*3
                        delta_x = delta_x + 0.5*inf[0]*dts*returned_flowX #m3/d*d*g/m3 = g
                        delta_s = delta_s + 0.5*inf[0]*dts*self.Sbulk[-1] #m3/d*d*g/m3 = g
    
                    #Get recirculation flow
                    if i+1 in self.recircT.keys(): #Flow coming into i
                        Q = self.recircT[i+1][1]
                        totalQin[i] = totalQin[i] + Q
                        origin = self.recircT[i+1][0]-1
                        delta_x = delta_x + Q*dts*self.Xbulk[origin] #m3/d*d*g/m3 = g
                        delta_s = delta_s + Q*dts*self.Sbulk[origin] #m3/d*d*g/m3 = g
                    if i+1 in self.recircO.keys(): #Flow going out from i
                        Q = self.recircO[i+1][1]
                        totalQin[i] = totalQin[i] - Q
                        delta_x = delta_x - Q*dts*self.Xbulk[i] #m3/d*d*g/m3 = g
                        delta_s = delta_s - Q*dts*self.Sbulk[i] #m3/d*d*g/m3 = g
        
                    #Get flow from previous
                    if i > 0:
                        Q = totalQin[i-1]
                        totalQin[i] = totalQin[i] + Q
                        delta_x = delta_x + Q*dts*self.Xbulk[i-1] #m3/d*d*g/m3 = g
                        delta_s = delta_s + Q*dts*self.Sbulk[i-1] #m3/d*d*g/m3 = g
        
                    #Subtract effluent and change to concentration
                    delta_x = delta_x - totalQin[i]*dts*self.Xbulk[i] #m3/d*d*g/m3 = g
                    delta_s = delta_s - totalQin[i]*dts*self.Sbulk[i] #m3/d*d*g/m3 = g
                    delta_x = delta_x/self.volumes[i] #g/m3
                    delta_s = delta_s/self.volumes[i] #g/m3

                    #Get biorates bulk
                    br = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1)) #g/m3-d
                    brX = br*Xpm
                    brS = br*Spm
                    delta_x = delta_x + brX.sum(axis=0)*dts
                    delta_s = delta_s + brS.sum(axis=0)*dts

                    ##CALCULATE NEW BULK CONCENTRATIONS
                    self.Xbulk[i] = self.Xbulk[i] + delta_x
                    if self.Xbulk[i].sum() > VSS:
                        self.Xbulk[i] = VSS*self.Xbulk[i]/(self.Xbulk[i].sum())
                    self.Xbulk[i][self.Xbulk[i] < 0] = 0
                    self.Sbulk[i] = self.Sbulk[i] + delta_s
                    self.Sbulk[i][0] = self.DO[i]
                    self.Sbulk[i][self.Sbulk[i] < 0] = 0

                    ##BIOFILM
                    
                    #Calculate substrate profiles
                    if self.areas[i] > 0:
                        for bfti in range(bf_tsteps):
                            bf_br = bio_rates(xpar=BF_x[i].T, spar=BF_s[i].T)
                            for s in range(len(sparlist)):
                                d_DBL = dtbf[i]*Dbf[s]/(self.DBL[i]**2)
                                d_BFL = dtbf[i]*Dbf[s]/(self.BFL[i]**2)
                                bf_brS = bf_br*Spm[:, s].reshape(-1,1)
                                bf_brS = bf_brS.sum(axis=0) #g/m3.d change
                                if len(BF_s[i][:, 0] == 1):
                                    BF_s[i][0,s] = BF_s[i][0,s] + d_DBL*(self.Sbulk[i][s]-BF_s[i][0,s])+dtbf[i]*bf_brS[0]
                                elif len(BF_s[i][:, 0] == 2):
                                    BF_s[i][0,s] = BF_s[i][0,s] + d_DBL*(self.Sbulk[i][s]-BF_s[i][0,s])+d_BFL*(BF_s[i][1,s]-BF_s[i][0,s])+dtbf[i]*bf_brS[0]
                                    BF_s[i][1,s] = BF_s[i][1,s] + d_BFL*(BF_s[i][0,s]-BF_s[i][1,s])+dtbf[i]*bf_brS[1]
                                else:
                                    BF_s[i][0,s] = BF_s[i][0,s] + d_DBL*(self.Sbulk[i][s]-BF_s[i][0,s])+d_BFL*(BF_s[i][1,s]-BF_s[i][0,s])+dtbf[i]*bf_brS[0]
                                    BF_s[i][-1,s] = BF_s[i][-1,s] + d_BFL*(BF_s[i][-2,s]-BF_s[i][-1,s])+dtbf[i]*bf_brS[-1]
                                    BF_s[i][1:-1,s] = BF_s[i][1:-1,s] + d_BFL*(BF_s[i][:-2,s]+BF_s[i][2:,s]-2*BF_s[i][1:-1,s])+dtbf[i]*bf_brS[1:-1]

                        #Update S in bulk due to biofilm reactions
                        for s in range(len(sparlist)):
                            bf_brS = bf_br*Spm[:, s].reshape(-1,1)
                            bf_brS = bf_brS.sum(axis=0) #g/m3.d change
                            sconsumed = bf_brS*self.BFL[i]*(10**-6) #g/m2.d per layer
                            sconsumed = sconsumed.sum() #sum all layers
                            self.Sbulk[i][s] = self.Sbulk[i][s] + dts*sconsumed*self.areas[i]/self.volumes[i]
                        self.Sbulk[i][0] = self.DO[i]
                        self.Sbulk[i][self.Sbulk[i] < 0] = 0

            #Update biofilm biomass
            for i in range(len(self.volumes)):
                if self.areas[i] > 0:
                    self.Sbiofilm[i][len(self.Sbiofilm[i][:,0])-BF_len[i]:] = BF_s[i]
                    self.Sbiofilm[i][:len(self.Sbiofilm[i][:,0])-BF_len[i]] = np.tile(self.Sbulk[i], (len(self.Sbiofilm[i][:,0])-BF_len[i], 1))
                    bf_br = bio_rates(xpar=self.Xbiofilm[i].T, spar=self.Sbiofilm[i].T)
                    for x in range(len(xparlist)):
                        bf_brX = bf_br*Xpm[:, x].reshape(-1,1)
                        bf_brX = bf_brX.sum(axis=0)
                        self.Xbiofilm[i][:, x] = self.Xbiofilm[i][:, x] + dtx*bf_brX
                    for j in range(len(self.Xbiofilm[i][:,0])-1,0,-1):
                        excessX = self.Xbiofilm[i][j,:].sum() - self.VSmax[i]
                        if excessX > 0:
                            self.Xbiofilm[i][j-1,:] = self.Xbiofilm[i][j-1,:] + excessX*self.Xbiofilm[i][j,:]/(self.Xbiofilm[i][j,:].sum())
                            self.Xbiofilm[i][j,:] = self.VSmax[i]*self.Xbiofilm[i][j,:]/(self.Xbiofilm[i][j,:].sum())
                        elif excessX < 0 and self.Xbiofilm[i][j-1,:].sum() > 0:
                            self.Xbiofilm[i][j,:] = self.Xbiofilm[i][j,:] + abs(excessX)*self.Xbiofilm[i][j-1,:]/(self.Xbiofilm[i][j-1,:].sum())
                            self.Xbiofilm[i][j-1,:] = self.Xbiofilm[i][j-1,:] - abs(excessX)*self.Xbiofilm[i][j-1,:]/(self.Xbiofilm[i][j-1,:].sum())
                        if self.Xbiofilm[i][j-1,:].sum() == 0:
                            break
                    if self.Xbiofilm[i][0,:].sum() > self.VSmax[i]:
                        excessX = self.Xbiofilm[i][0,:].sum() - self.VSmax[i]
                        self.Xbiofilm[i][0,:] = self.VSmax[i]*self.Xbiofilm[i][0,:]/(self.Xbiofilm[i][0,:].sum())
                        self.Xbulk[i] = self.Xbulk[i] + (self.BFL[i]*10**-6)*(self.areas[i])*(excessX*self.Xbiofilm[i][0,:]/(self.Xbiofilm[i][0,:].sum()))/self.volumes[i]
                        BF_Xprod[i] = BF_Xprod[i] + (self.BFL[i]*10**-6)*(excessX*self.Xbiofilm[i][0,:]/(self.Xbiofilm[i][0,:].sum()))

                    bf_br = bio_rates(xpar=self.Xbiofilm[i].T, spar=self.Sbiofilm[i].T)
                    for s in range(len(sparlist)):
                        bf_brS = bf_br*Spm[:, s].reshape(-1,1)
                        bf_brS = bf_brS.sum(axis=0) #g/m3.d change
                        bf_brS = bf_brS*(self.BFL[i]*10**-6)
                        BF_Sprod[i][s] = BF_Sprod[i][s] + dtx*bf_brS.sum()

            #Calculate VSS prod
            for i in range(len(self.volumes)):
                if VSS > 0:
                    if self.Xbulk[i].sum() > VSS:
                        self.Xbulk[i] = VSS*self.Xbulk[i]/(self.Xbulk[i].sum())
                    br = bio_rates(xpar=self.Xbulk[i], spar=self.Sbulk[i]).reshape((-1,1)) #g/m3-d
                    brX = br*Xpm
                    brS = br*Spm
                    Bulk_Xprod[i] = Bulk_Xprod[i] + brX.sum(axis=0)*dtx
                    Bulk_Sprod[i] = Bulk_Sprod[i] + brS.sum(axis=0)*dtx

        #Update loggers
        for i in range(len(self.volumes)):
            self.bulklogger[str(i+1)]['Conc'].loc[self.time, xparlist+sparlist] = np.concatenate([self.Xbulk[i], self.Sbulk[i]])
            self.bulklogger[str(i+1)]['Prod'].loc[self.time, xparlist+sparlist] = np.concatenate([Bulk_Xprod[i]/timestep, Bulk_Sprod[i]/timestep])
            if self.areas[i] > 0:
                self.biofilmlogger[str(i+1)]['Conc'][self.time] = pd.DataFrame(np.concatenate([self.Xbiofilm[i], self.Sbiofilm[i]], axis=1), index=np.arange(0, self.Lmax[i], self.BFL[i])+self.BFL[i]/2, columns=xparlist+sparlist)
                self.biofilmlogger[str(i+1)]['Prod'].loc[self.time, xparlist+sparlist] = np.concatenate([BF_Xprod[i]/timestep, BF_Sprod[i]/timestep])

    def view(self, savename=None):
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


