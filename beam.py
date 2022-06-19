from warp import *

UB = {
    'pencil': {
        'ibeam': -0.6e-3,
        'emitx': 7.6e-6,
        'emity': 7.6e-6,
        'xcent': 0.0,
        'xpcent': 0.0,
        'aper': [0.00025, 0.00025, -0.00127, -0.00127],
        'ring': {'100%': [0.001692, 0.001716, -0.006587, 0.007104],
                 '94%':  [0.001659, 0.001482, -0.008613, 0.007757],
                 '83%':  [0.001717, 0.001663, -0.007678, 0.007346],
                 '69%':  [0.001845, 0.001931, -0.006684, 0.006747],
                }
        },
    '6mA': {
        'ibeam': -6.0e-3,
        'emitx': 25.5e-6,
        'emity': 25.5e-6,
        'xcent': 0.0005,
        'xpcent': 0.0,
        'aper': [0.000875, 0.000875, -0.0043, -0.0043],
        'ring': {'100%': [0.003218, 0.003085, -0.016961, 0.015836],
                 '94%':  [0.003308, 0.003281, -0.016084, 0.015544],
                 '83%':  [0.003533, 0.003679, -0.014735, 0.014916],
                 '69%':  [0.003953, 0.004356, -0.013312, 0.014167],
                }
        },
    '23mA': {
        'ibeam': -21.0e-3,
        'emitx': 30.0e-6,
        'emity': 30.0e-6,
        'xcent': 0.0005,
        'xpcent': 0.0,
        'aper': [0.0015, 0.0015, -0.0067, -0.0067],
        'ring': {'100%': [0.004578, 0.004775, -0.021763, 0.022386],
                 '94%':  [0.004825, 0.005107, -0.021347, 0.022279],
                 '83%':  [0.005364, 0.005803, -0.020576, 0.022015],
                 '69%':  [0.006286, 0.007065, -0.019748, 0.021923],
                }
        },
    '80mA': {
        'ibeam': -78.0e-3,
        'emitx': 58.9e-6,
        'emity': 58.9e-6,
        'xcent': 0.0005,
        'xpcent': 0.0,
        'aper': [0.00285, 0.00285, -0.01273, -0.01273],
        'ring': {'100%': [0.008366, 0.008814, -0.038605, 0.040407],
                 '94%':  [0.008856, 0.009405, -0.038036, 0.040322],
                 '83%':  [0.009916, 0.010764, -0.037323, 0.040331],
                 '69%':  [0.011680, 0.013153, -0.036182, 0.040585],
                }
        },
    '100mA': {
        'ibeam': -104.0e-3,
        'emitx': 64.0e-6,
        'emity': 64.0e-6,
        'xcent': 0.0005,
        'xpcent': 0.0,
        'aper': [0.0032, 0.0032, -0.0143, -0.0143],
        'ring': {'100%': [0.009570, 0.010084, -0.043915, 0.046068],
                 '94%':  [0.010131, 0.010769, -0.043315, 0.046025],
                 '83%':  [0.011369, 0.012309, -0.042482, 0.045980],
                 '69%':  [0.013379, 0.015061, -0.041343, 0.046110],
                }
        },
}


def createBeam(aperture='6mA'):
    # --- Create the beam species
    beam = Species(type=Electron,name="beam species",charge_state=-1)

    b = UB[aperture]
    # --- Set input parameters describing the beam,
    beam.emit     = b['emitx'] # emittance
    beam.a0       = b['ring']['100%'][0] # beam edge in x-x'
    beam.b0       = b['ring']['100%'][1] # beam edge in y-y'
    beam.ap0      = b['ring']['100%'][2]
    beam.bp0      = b['ring']['100%'][3]
    beam.ibeam    = b['ibeam'] 
    beam.ekin     = 10000.0
    beam.aion     = emass/amu
    top.lrelativ  = True
    derivqty()
    #beam.ibeam    = beam.ibeam/(top.gammabar**2)
    e_spread = 100
    beam.vthz = beam.vbeam*e_spread/(2.*beam.ekin)
    
    return beam
