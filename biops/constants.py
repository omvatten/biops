import pandas as pd

aob = {}
aob['mu_max'] = 2.05 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
aob['b'] = 0.13 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
aob['K_O2'] = 0.6 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)
aob['K_NH4'] = 2.4 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)
aob['fs'] = 0.054
aob['CHON'] = [5,7,2,1]
aob['fI'] = 0.2

nob = {}
nob['mu_max'] = 1.45 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
nob['b'] = 0.06 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
nob['K_O2'] = 0.4 #gN/m3 (Laureni et al. WR 154, 104-116, 2019)
nob['K_NO2'] = 0.5 #g/m3 (Laureni et al. WR 154, 104-116, 2019)
nob['fs'] = 0.07
nob['CHON'] = [5,7,2,1]
nob['fI'] = 0.2

anammox = {}
anammox['mu_max'] = 0.08 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
anammox['b'] = 0.003 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_O2'] = 0.01 #gO2/m3 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_NO2'] = 0.05 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)
anammox['K_NH4'] = 0.07 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)
anammox['fs'] = 0.051
anammox['CHON'] = [1,1.74,0.31,0.2]
anammox['fI'] = 0.2

comammox = {}
comammox['mu_max'] = 0.14 #d-1 (Kits et al. Nature 549, 269, 2017; 14.8 umolN/mgProt.h, 400 mgProt/molN)
comammox['b'] = 0.003 #d-1 (assumed)
comammox['K_O2'] = 0.4 #gO2/m3 (assumed)
comammox['K_NH4'] = 0.012 #gN/m3 (Kits et al. Nature 549, 269, 2017; 0.84 uM)
comammox['Ki_NH4'] = 3.44 #gN/m3 (Sakoula et al. ISME 15, 1010, 2020; 246 uM)
comammox['fs'] = 0.02
comammox['CHON'] = [5,7,2,1]
comammox['fI'] = 0.2

oho = {}
oho['mu_max_O2'] = 6 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['mu_max_NOx'] = 4.8 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['b'] = 0.62 #d-1 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_s'] = 20 #g/m3 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_O2'] = 0.2 #gO2/m3 (Li et al. STOTEN 669, 683-691, 2019)
oho['K_NOx'] = 0.3 #gN/m3 (Li et al. STOTEN 669, 683-691, 2019)
oho['fs_O2'] = 0.67
oho['fs_NOx'] = 0.67
oho['substrateCHON'] = [1,2,1,0]
oho['CHON'] = [5,7,2,1]
oho['fI'] = 0.2

## Functions for calculating yields
def Y_aob(fs=aob['fs'], CHON=aob['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]

    #Energy yielding reaction
    nh4_ox = 1/6 #mol NH4+ to NO2

    #Catabolism
    o2_red = (1-fs)/4 #mol O2 reduced

    #Anabolism    
    hco3_bm = fs/(4-Cox) #mol HCO3- to build biomass
    bm = hco3_bm/CHON[0] #mol CHON produced
    nh4_bm = bm*CHON[3] #mol NH4+ into biomass
    
    nh4 = nh4_ox+nh4_bm
    h2o = o2_red*2+hco3_bm*3-bm*CHON[2]-nh4_ox*2 #mol H2O produced
    hplus = nh4*4 + hco3_bm-bm*CHON[1]-h2o*2 #mol H+ produced

    o2p = str(round(o2_red/nh4, dec))
    hco3p = str(round(hco3_bm/nh4, dec))
    bmp = str(round(bm/nh4, dec))
    no2p = str(round(nh4_ox/nh4, dec))
    h2op = str(round(abs(h2o)/nh4, dec))
    hp = str(round(abs(hplus)/nh4, dec))

    reactants = 'NH4+ + '+o2p+' O2 + '+hco3p+' HCO3-'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+no2p+' NO2-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gNH4-N'] = (bm*MWbm)/(nh4*14) #gVSS/gNH4-N
    out['Y_gCODx_gNH4-N'] = (bm*(4-Cox)*CHON[0]*8)/(nh4*14) #gCODx/gNH4-N
    out['Y_gNO2-N_gNH4-N'] = (nh4_ox*14)/(nh4*14) #gNO2-N/gNH4-N
    out['Y_gO2_gNH4-N'] = -(o2_red*32)/(nh4*14) #gO2/gNH4-N 
    out['Y_meqAlk_gNH4-N'] = -(hco3_bm+hplus)/(nh4*14) #gO2/gNH4-N 
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def Y_nob(fs=nob['fs'], CHON=nob['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]
    
    hco3_bm = fs/(4-Cox)  #mol HCO3- to build biomass
    bm = hco3_bm/CHON[0] #mol CHON produced
    nh4_bm = bm*CHON[3]
    no2 = 1/2 #mol NO2 to NO3
    o2_red = (1-fs)/4 #mol O2 reduced
    h2o = o2_red*2+hco3_bm*3+no2*2-bm*CHON[2]-no2*3 #mol H2O produced
    hplus = hco3_bm+nh4_bm*4-bm*CHON[1]-2*h2o #mol H+ produced

    o2p = str(round(o2_red/no2, dec))
    hco3p = str(round(hco3_bm/no2, dec))
    nh4p = str(round(nh4_bm/no2, dec))
    bmp = str(round(bm/no2, dec))
    no3p = str(round(no2/no2, dec))
    h2op = str(round(abs(h2o)/no2, dec))
    hp = str(round(abs(hplus)/no2, dec))

    reactants = 'NO2- + '+o2p+' O2 + '+hco3p+' HCO3- + '+nh4p+' NH4+'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+no3p+' NO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gNO2-N'] = (bm*MWbm)/(no2*14) #gVSS/gNO2-N 
    out['Y_gCODx_gNO2-N'] = (bm*(4-Cox)*CHON[0]*8)/(no2*14) #gCODx/gNO2-N 
    out['Y_gNO3-N_gNO2-N'] = (no2*14)/(no2*14) #gNO3-N/gNO2-N 
    out['Y_gNH4-N_gNO2-N'] = -(nh4_bm*14)/(no2*14) #gNH4-N/gNO2-N
    out['Y_gO2_gNO2-N'] = -(o2_red*32)/(no2*14) #gO2/gNO2-N 
    out['Y_meqAlk_gNO2-N'] = -(hco3_bm+hplus)/(no2*14) #meq consumed/gNO2-N
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def Y_anammox(fs=anammox['fs'], CHON=anammox['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]

    hco3_bm = fs/(4-Cox) #mol HCO3- to build biomass
    bm = hco3_bm/CHON[0] #mol CHON produced
    nh4_bm = bm*CHON[3] #mol NH4+ into biomass
    no2_ox = fs/2 #mol NO2 oxidized to NO3 to generate e- for HCO3 reduction
    nh4_ox = (1-fs)/6 #mol NH4 oxidized to NO2
    no2_red = (1-fs)/3 #mol NO2 reduced to N2
    no2 = no2_ox + no2_red - nh4_ox
    nh4 = nh4_bm + nh4_ox
    n2 = no2_red/2

    h2o = no2*2+hco3_bm*3-bm*CHON[2]-no2_ox*3 #mol H2O produced
    hplus = nh4*4+hco3_bm-bm*CHON[1]-h2o*2 #mol H+ produced

    no2p = str(round(no2/nh4, dec))
    hco3p = str(round(hco3_bm/nh4, dec))
    bmp = str(round(bm/nh4, dec))
    no3p = str(round(no2_ox/nh4, dec))
    h2op = str(round(abs(h2o)/nh4, dec))
    hp = str(round(abs(hplus)/nh4, dec))
    n2p = str(round(n2/nh4, dec))

    reactants = 'NH4+ + '+no2p+' NO2- + '+hco3p+' HCO3-'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+n2p+' N2 + '+no3p+' NO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gNH4-N'] = (bm*MWbm)/(nh4*14) #gVSS/gNH4-N 
    out['Y_gCODx_gNH4-N'] = (bm*(4-Cox)*CHON[0]*8)/(nh4*14) #gCODx/gNH4-N 
    out['Y_gNO2-N_gNH4-N'] = -(no2*14)/(nh4*14) #gNO2-N/gNH4-N 
    out['Y_gNO3-N_gNH4-N'] = (no2_ox*14)/(nh4*14) #gNO3-N/gNH4-N 
    out['Y_meqAlk_gNH4-N'] = -(hco3_bm+hplus)/(nh4*14) #meq/gNH4-N
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

#(Kits et al. Nature 549, 269, 2017; 400 mg protein/mol NH3)
def Y_comammox(fs=comammox['fs'], CHON=comammox['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]

    #Energy yielding reaction
    nh4_ox = 1/8 #mol NH4+ to NO3

    #Catabolism
    o2_red = (1-fs)/4 #mol O2 reduced

    #Anabolism    
    hco3_bm = fs/(4-Cox) #mol HCO3- to build biomass
    bm = hco3_bm/CHON[0] #mol CHON produced
    nh4_bm = bm*CHON[3] #mol NH4+ into biomass
    
    nh4 = nh4_ox+nh4_bm
    h2o = o2_red*2+hco3_bm*3-bm*CHON[2]-nh4_ox*3 #mol H2O produced
    hplus = nh4*4 + hco3_bm-bm*CHON[1]-h2o*2 #mol H+ produced

    o2p = str(round(o2_red/nh4, dec))
    hco3p = str(round(hco3_bm/nh4, dec))
    bmp = str(round(bm/nh4, dec))
    no3p = str(round(nh4_ox/nh4, dec))
    h2op = str(round(abs(h2o)/nh4, dec))
    hp = str(round(abs(hplus)/nh4, dec))

    reactants = 'NH4+ + '+o2p+' O2 + '+hco3p+' HCO3-'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+no3p+' NO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gNH4-N'] = (bm*MWbm)/(nh4*14) #gVSS/gNH4-N
    out['Y_gCODx_gNH4-N'] = (bm*(4-Cox)*CHON[0]*8)/(nh4*14) #gCODx/gNH4-N
    out['Y_gNO3-N_gNH4-N'] = (nh4_ox*14)/(nh4*14) #gNO3-N/gNH4-N
    out['Y_gO2_gNH4-N'] = -(o2_red*32)/(nh4*14) #gO2/gNH4-N 
    out['Y_meqAlk_gNH4-N'] = -(hco3_bm+hplus)/(nh4*14) #gO2/gNH4-N 
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def Y_oho_O2(fs=oho['fs_O2'], substrateCHON=oho['substrateCHON'], CHON=oho['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]
    substrateCox = (substrateCHON[1]*(-1)+substrateCHON[2]*(2)+substrateCHON[3]*(3))/substrateCHON[0]
    
    hco3_red = fs/(4-Cox)  #mol HCO3- to build biomass
    hco3_prod = 1/(4-substrateCox)  #mol HCO3- produced from subst oxidation
    hco3 = hco3_prod - hco3_red

    bm = hco3_red/CHON[0] #mol CHON produced
    sub = hco3_prod/substrateCHON[0] #mol substrateCHON oxidized
    nh4_bm = bm*CHON[3]
    o2_red = (1-fs)/4 #mol O2 reduced
    h2o = o2_red*2+sub*substrateCHON[2]-bm*CHON[2]-hco3*3 #mol H2O produced
    hplus = sub*substrateCHON[1]+nh4_bm*4-bm*CHON[1]-2*h2o-hco3 #mol H+ produced

    o2p = str(round(o2_red/sub, dec))
    nh4p = str(round(nh4_bm/sub, dec))
    bmp = str(round(bm/sub, dec))
    hco3p = str(round(hco3/sub, dec))
    h2op = str(round(abs(h2o)/sub, dec))
    hp = str(round(abs(hplus)/sub, dec))

    subcomp = 'C'+str(substrateCHON[0])+'H'+str(substrateCHON[1])+'O'+str(substrateCHON[2])+'N'+str(substrateCHON[3])
    reactants = subcomp+' + '+o2p+' O2 + '+nh4p+'NH4+'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+hco3p+' HCO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14
    COD = sub*substrateCHON[0]*(4-substrateCox)*8

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gCOD'] = (bm*MWbm)/COD #gVSS/gCOD 
    out['Y_gCODx_gCOD'] = (bm*(4-Cox)*CHON[0]*8)/COD #gCODx/gCOD 
    out['Y_gNH4-N_gCOD'] = -(nh4_bm*14)/COD #gNH4-N/gCOD 
    out['Y_gO2_gCOD'] = -(o2_red*32)/COD #gO2/gCOD 
    out['Y_meqAlk_gCOD'] = (hco3-hplus)/COD #meq/gCOD 
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def Y_oho_NO2(fs=oho['fs_NOx'], substrateCHON=oho['substrateCHON'], CHON=oho['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]
    substrateCox = (substrateCHON[1]*(-1)+substrateCHON[2]*(2)+substrateCHON[3]*(3))/substrateCHON[0]
    
    hco3_red = fs/(4-Cox)  #mol HCO3- to build biomass
    hco3_prod = 1/(4-substrateCox)  #mol HCO3- produced from subst oxidation
    hco3 = hco3_prod - hco3_red

    bm = hco3_red/CHON[0] #mol CHON produced
    sub = hco3_prod/substrateCHON[0] #mol substrateCHON oxidized
    nh4_bm = bm*CHON[3]
    no2_red = (1-fs)/3 #mol NO2- reduced
    h2o = no2_red*2+sub*substrateCHON[2]-bm*CHON[2]-hco3*3 #mol H2O produced
    hplus = sub*substrateCHON[1]+nh4_bm*4-bm*CHON[1]-2*h2o-hco3 #mol H+ produced

    no2p = str(round(no2_red/sub, dec))
    nh4p = str(round(nh4_bm/sub, dec))
    bmp = str(round(bm/sub, dec))
    hco3p = str(round(hco3/sub, dec))
    n2p = str(round((no2_red/2)/sub, dec))
    h2op = str(round(abs(h2o)/sub, dec))
    hp = str(round(abs(hplus)/sub, dec))

    subcomp = 'C'+str(substrateCHON[0])+'H'+str(substrateCHON[1])+'O'+str(substrateCHON[2])+'N'+str(substrateCHON[3])
    reactants = subcomp+' + '+no2p+' NO2- + '+nh4p+'NH4+'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+n2p+' N2 + '+hco3p+' HCO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14
    COD = sub*substrateCHON[0]*(4-substrateCox)*8

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gCOD'] = (bm*MWbm)/COD #gVSS/gCOD 
    out['Y_gCODx_gCOD'] = (bm*(4-Cox)*CHON[0]*8)/COD #gCODx/gCOD
    out['Y_gNH4-N_gCOD'] = -(nh4_bm*14)/COD #gNH4-N/gCOD 
    out['Y_gNO2-N_gCOD'] = -(no2_red*14)/COD #gNO2-N/gCOD 
    out['Y_meqAlk_gCOD'] = (hco3-hplus)/COD #meq/gCOD
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def Y_oho_NO3(fs=oho['fs_NOx'], substrateCHON=oho['substrateCHON'], CHON=oho['CHON'], dec=2):
    Cox = (CHON[1]*(-1)+CHON[2]*(2)+CHON[3]*(3))/CHON[0]
    substrateCox = (substrateCHON[1]*(-1)+substrateCHON[2]*(2)+substrateCHON[3]*(3))/substrateCHON[0]
    
    hco3_red = fs/(4-Cox)  #mol HCO3- to build biomass
    hco3_prod = 1/(4-substrateCox)  #mol HCO3- produced from subst oxidation
    hco3 = hco3_prod - hco3_red

    bm = hco3_red/CHON[0] #mol CHON produced
    sub = hco3_prod/substrateCHON[0] #mol substrateCHON oxidized
    nh4_bm = bm*CHON[3]
    no3_red = (1-fs)/5 #mol NO3- reduced
    h2o = no3_red*3+sub*substrateCHON[2]-bm*CHON[2]-hco3*3 #mol H2O produced
    hplus = sub*substrateCHON[1]+nh4_bm*4-bm*CHON[1]-2*h2o-hco3 #mol H+ produced

    no3p = str(round(no3_red/sub, dec))
    nh4p = str(round(nh4_bm/sub, dec))
    bmp = str(round(bm/sub, dec))
    hco3p = str(round(hco3/sub, dec))
    n2p = str(round((no3_red/2)/sub, dec))
    h2op = str(round(abs(h2o)/sub, dec))
    hp = str(round(abs(hplus)/sub, dec))

    subcomp = 'C'+str(substrateCHON[0])+'H'+str(substrateCHON[1])+'O'+str(substrateCHON[2])+'N'+str(substrateCHON[3])
    reactants = subcomp+' + '+no3p+' NO3- + '+nh4p+'NH4+'

    bmcomp = 'C'+str(CHON[0])+'H'+str(CHON[1])+'O'+str(CHON[2])+'N'+str(CHON[3])
    products = bmp+' '+bmcomp+' + '+n2p+' N2 + '+hco3p+' HCO3-'
    if h2o < 0:
        reactants = reactants + ' + '  + h2op + ' H2O'
    elif h2o > 0:
        products = products + ' + '  + h2op + ' H2O'
    if hplus < 0:
        reactants = reactants + ' + '  + hp + ' H+'
    elif hplus > 0:
        products = products + ' + '  + hp + ' H+'
    
    eq = reactants+' --> '+products
    MWbm = CHON[0]*12+CHON[1]+CHON[2]*16+CHON[3]*14
    COD = sub*substrateCHON[0]*(4-substrateCox)*8

    out = {}
    out['eq'] = eq
    out['Y_gVSS_gCOD'] = (bm*MWbm)/COD #gVSS/gCOD
    out['Y_gCODx_gCOD'] = (bm*(4-Cox)*CHON[0]*8)/COD #gCODx/gCOD 
    out['Y_gNH4-N_gCOD'] = -(nh4_bm*14)/COD #gNH4-N/gCOD 
    out['Y_gNO3-N_gCOD'] = -(no3_red*14)/COD #gNO3-N/gCOD 
    out['Y_meqAlk_gCOD'] = (hco3-hplus)/COD #meq/gCOD 
    out['i_gN_gVSS'] = CHON[3]*14/MWbm
    out['i_gN_gCODx'] = CHON[3]*14/((4-Cox)*CHON[0]*8)
    return out

def get_matrix(values=None):
    processlist = ['AOB_growth', 'AOB_decay', 'NOB_growth', 'NOB_decay', 'AMX_growth', 'AMX_decay', 'CMX_growth', 'CMX_decay', 'OHO_growth_O2', 'OHO_growth_NO2', 'OHO_growth_NO3', 'OHO_decay']
    cols = ['X_AOB', 'X_NOB', 'X_AMX', 'X_CMX', 'X_OHO', 'X_I', 'S_NH4', 'S_NO2', 'S_NO3', 'S_s', 'S_O2', 'Equation']
    proc_matrix = pd.DataFrame(0, index=processlist, columns=cols)

    #AOB
    if values == None:
        y = Y_aob()
        y['fI'] = aob['fI']
    else:
        y = Y_aob(fs=values['aob']['fs'], CHON=values['aob']['CHON'])
        y['fI'] = values['aob']['fI']
    proc_matrix.loc['AOB_growth', 'X_AOB'] = 1
    proc_matrix.loc['AOB_growth', 'S_NH4'] = -1/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AOB_growth', 'S_NO2'] = y['Y_gNO2-N_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AOB_growth', 'S_O2'] = y['Y_gO2_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AOB_growth', 'Equation'] = 'X_AOB*mu_max*(S_NH4/(K_NH4+S_NH4))*(S_O2/(K_O2+S_O2))'
    proc_matrix.loc['AOB_decay', 'X_AOB'] = -1
    proc_matrix.loc['AOB_decay', 'X_I'] = y['fI']
    proc_matrix.loc['AOB_decay', 'S_s'] = 1-y['fI']
    proc_matrix.loc['AOB_decay', 'Equation'] = 'X_AOB*b'

    #NOB
    if values == None:
        y = Y_nob()
        y['fI'] = nob['fI']
    else:
        y = Y_nob(fs=values['nob']['fs'], CHON=values['nob']['CHON'])
        y['fI'] = values['nob']['fI']
    proc_matrix.loc['NOB_growth', 'X_NOB'] = 1
    proc_matrix.loc['NOB_growth', 'S_NO2'] = -1/y['Y_gVSS_gNO2-N']
    proc_matrix.loc['NOB_growth', 'S_NO3'] = y['Y_gNO3-N_gNO2-N']/y['Y_gVSS_gNO2-N']
    proc_matrix.loc['NOB_growth', 'S_NH4'] = y['Y_gNH4-N_gNO2-N']/y['Y_gVSS_gNO2-N']
    proc_matrix.loc['NOB_growth', 'S_O2'] = y['Y_gO2_gNO2-N']/y['Y_gVSS_gNO2-N']
    proc_matrix.loc['NOB_growth', 'Equation'] = 'X_NOB*mu_max*(S_NO2/(K_NO2+S_NO2))*(S_O2/(K_O2+S_O2))'
    proc_matrix.loc['NOB_decay', 'X_NOB'] = -1
    proc_matrix.loc['NOB_decay', 'X_I'] = y['fI']
    proc_matrix.loc['NOB_decay', 'S_s'] = 1-y['fI']
    proc_matrix.loc['NOB_decay', 'Equation'] = 'X_NOB*b'

    #AMX
    if values == None:
        y = Y_anammox()
        y['fI'] = anammox['fI']
    else:
        y = Y_anammox(fs=values['anammox']['fs'], CHON=values['anammox']['CHON'])
        y['fI'] = values['anammox']['fI']
    proc_matrix.loc['AMX_growth', 'X_AMX'] = 1
    proc_matrix.loc['AMX_growth', 'S_NH4'] = -1/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AMX_growth', 'S_NO2'] = y['Y_gNO2-N_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AMX_growth', 'S_NO3'] = y['Y_gNO3-N_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['AMX_growth', 'Equation'] = 'X_AMX*mu_max*(S_NH4/(K_NH4+S_NH4))*(S_NO2/(K_NO2+S_NO2))*(K_O2/(K_O2+S_O2))'
    proc_matrix.loc['AMX_decay', 'X_AMX'] = -1
    proc_matrix.loc['AMX_decay', 'X_I'] = y['fI']
    proc_matrix.loc['AMX_decay', 'S_s'] = 1-y['fI']
    proc_matrix.loc['AMX_decay', 'Equation'] = 'X_AMX*b'

    #CMX
    if values == None:
        y = Y_comammox()
        y['fI'] = comammox['fI']
    else:
        y = Y_comammox(fs=values['comammox']['fs'], CHON=values['comammox']['CHON'])
        y['fI'] = values['comammox']['fI']
    proc_matrix.loc['CMX_growth', 'X_AOB'] = 1
    proc_matrix.loc['CMX_growth', 'S_NH4'] = -1/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['CMX_growth', 'S_NO3'] = y['Y_gNO3-N_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['CMX_growth', 'S_O2'] = y['Y_gO2_gNH4-N']/y['Y_gVSS_gNH4-N']
    proc_matrix.loc['CMX_growth', 'Equation'] = 'X_CMX*mu_max*(S_NH4/(K_NH4+S_NH4+S_NH4/Ki_NH4^2))*(S_O2/(K_O2+S_O2))'
    proc_matrix.loc['CMX_decay', 'X_CMX'] = -1
    proc_matrix.loc['CMX_decay', 'X_I'] = y['fI']
    proc_matrix.loc['CMX_decay', 'S_s'] = 1-y['fI']
    proc_matrix.loc['CMX_decay', 'Equation'] = 'X_CMX*b'

    #OHO
    if values == None:
        y = Y_oho_O2()
        y['fI'] = oho['fI']
    else:
        y = Y_oho_O2(fs=values['oho']['fs_O2'], substrateCHON=values['oho']['substrateCHON'], CHON=values['oho']['CHON'])
        y['fI'] = values['oho']['fI']
    proc_matrix.loc['OHO_growth_O2', 'X_OHO'] = 1
    proc_matrix.loc['OHO_growth_O2', 'S_s'] = -1/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_O2', 'S_O2'] = y['Y_gO2_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_O2', 'S_NH4'] = y['Y_gNH4-N_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_O2', 'Equation'] = 'X_OHO*mu_max*(S_s/(K_s+S_s))*(S_O2/(K_O2+S_O2))'

    if values == None:
        y = Y_oho_NO2()
        y['fI'] = oho['fI']
    else:
        y = Y_oho_NO2(fs=values['oho']['fs_NOx'], substrateCHON=values['oho']['substrateCHON'], CHON=values['oho']['CHON'])
        y['fI'] = values['oho']['fI']
    proc_matrix.loc['OHO_growth_NO2', 'X_OHO'] = 1
    proc_matrix.loc['OHO_growth_NO2', 'S_s'] = -1/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO2', 'S_NO2'] = y['Y_gNO2-N_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO2', 'S_NH4'] = y['Y_gNH4-N_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO2', 'Equation'] = 'X_OHO*mu_max*(S_s/(K_s+S_s))*(S_NO2/(K_NOx+S_NO2))*(K_O2/(K_O2+S_O2))'

    if values == None:
        y = Y_oho_NO3()
        y['fI'] = oho['fI']
    else:
        y = Y_oho_NO3(fs=values['oho']['fs_NOx'], substrateCHON=values['oho']['substrateCHON'], CHON=values['oho']['CHON'])
        y['fI'] = values['oho']['fI']
    proc_matrix.loc['OHO_growth_NO3', 'X_OHO'] = 1
    proc_matrix.loc['OHO_growth_NO3', 'S_s'] = -1/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO3', 'S_NO3'] = y['Y_gNO3-N_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO3', 'S_NH4'] = y['Y_gNH4-N_gCOD']/y['Y_gVSS_gCOD']
    proc_matrix.loc['OHO_growth_NO3', 'Equation'] = 'X_OHO*mu_max*(S_s(K_s+S_s))*(S_NO3/(K_NOx+S_NO3))*(K_O2/(K_O2+S_O2))'

    proc_matrix.loc['OHO_decay', 'X_OHO'] = -1
    proc_matrix.loc['OHO_decay', 'X_I'] = y['fI']
    proc_matrix.loc['OHO_decay', 'S_s'] = 1-y['fI']
    proc_matrix.loc['OHO_decay', 'Equation'] = 'X_OHO*b'
    
    return proc_matrix
