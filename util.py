from settings import params, constants, earth_field

import numpy as np

def c2k(magnet, current):
    '''
    convert a current to a quad strength k
    '''
    if magnet == 'QF':
        peak_grad = params['QF'][0]*params['QF'][2]
    elif magnet == 'QD':
        peak_grad = params['QD'][0]*params['QD'][2]
    elif magnet == 'QR1':
        peak_grad = params['QR1'][0]*params['QR1'][2]
    elif magnet == 'YQ':
        peak_grad = params['YQ'][0]*params['YQ'][2]
    else:
        peak_grad = 0
        
    return current * peak_grad / constants['ridg']
    
def c2angle(magnet, current, efield=1):
    '''
    convert a current to an angle for a dipole and returns field as well
    efield is a boolean whether to include efield or not
    
    A current = 2.9858e02 A corresponds to 10 degree bend w/ no earth field
    '''
    
    peak_field = params['D'][0]*params['D'][2]
    l_eff = params['D'][1]*params['D'][3]*1e2 # in cm
    if magnet == 'PD':
        peak_field = params['PD'][0]*params['PD'][2]
        l_eff = params['PD'][1]*params['PD'][3]*1e2 # in cm
    
    angleR = l_eff * peak_field * current / (constants['ridg']*1e2) # cm
    by = peak_field * current # field in tesla
    
    if efield:
        angleR += -1*earth_field[magnet][3]*1e-3*32 / (constants['ridg']*1e6) # average across 32 cm
        by += -1*earth_field[magnet][3]*1e-3*32*1e-4
        
    return angleR
        
def angle2c(magnet, angle, efield=1):
    '''
    convert an angle to a current, see c2angle
    angle should be in radians
    '''
    peak_field = params['D'][0]*params['D'][2]
    l_eff = params['D'][1]*params['D'][3]*1e2 # in cm
    if magnet == 'PD':
        peak_field = params['PD'][0]*params['PD'][2]
        l_eff = params['PD'][1]*params['PD'][3]*1e2 # in cm

    current = angle * constants['ridg']*1e2 / (l_eff * peak_field)
    if efield:
        current -= -1*earth_field[magnet][3]*1e-7*32 / (l_eff * peak_field)

    return current







    
