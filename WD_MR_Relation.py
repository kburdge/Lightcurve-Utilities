#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:59:42 2020

@author: kburdge
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def MR_Relation(m,T_eff):
    if m<0.448:
        mu=0.48
    elif m<0.503:
        mu=4.2*m-1.4
    else:
        mu=0.78*m+0.32
    
    si=0.984-0.021*m
    
    if m>0.0:
        mu_core=2.0
        R_ch=2.45354/mu_core*(m/(5.816/mu_core**2))**(-1/3)*(1-(m/(5.816/mu_core**2))**(4/3))**(1/2)
        R_WD=si*R_ch*(1-1/mu*(T_eff/588862)*m**(-1)*(si*R_ch))**(-1)

    return R_WD*0.0091577

def R_Relation(R,T_eff):
    m=np.linspace(0.16,1.4,10000)
    radii=[]
    for i in m:
        radii.append(MR_Relation(i,T_eff))
    f=interpolate.interp1d(radii,m)
    print(f(R))
    plt.plot(radii,m)
    
    return f(R)

