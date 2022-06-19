mpaths = {
    'u8.umd.edu':     '/home/WARP/Magnets/',
    'ebte.umd.edu':   '/home/WARP/Magnets/',
    'accel.umd.edu':  '/home/WARP/Magnets/',
    'marvin.umd.edu': '/work/shared/Magnets/',
    'umer.umd.edu':   '/ebte/pywarp/rscripts/',
    'scooby.localdomain':   '/ebte/pywarp/rscripts/',
}

# --- Locate Paths to pdb arrays
import os
#f mpaths.has_key(os.environ['HOST']):  Bbase = mpaths[os.environ['HOST']]
#lse:
#   Bbase = '/ebte/pywarp/rscripts/'
#   print """ !!! WARNING !!!!
#       It is much better to copy the magnet pdb files to your local machine
#       Be sure to edit this line to supply the correct path once done """

Bbase = '/work/haber/ring/respmat/Magnet'

Ebase = "/humer/shared/Universal/BEarth-data/"  # Location of Earth field tables

# --- Primary Magnets
#mRQ  = os.path.join(Bbase, 'RQUAD3.MAGLI.pdb')    # Main Ring Quad (3rd-gen)
mRQ  = os.path.join(Bbase, 'RQUAD3.MAGLI.pkl')    # Main Ring Quad (3rd-gen)
#mRQs = os.path.join(Bbase, 'RQUAD3S.MAGLI.pdb')   # Main Ring Quad (3rd-gen, Spiral)
mRQs = os.path.join(Bbase, 'RQUAD3S.MAGLI.pkl')   # Main Ring Quad (3rd-gen, Spiral)
#mRQg = os.path.join(Bbase, 'RQUAD3G.MAGLI.pdb')   # Main Ring Quad (3rd-gen, Spiral with glitch)
mRQg = os.path.join(Bbase, 'RQUAD3G.MAGLI.pkl')   # Main Ring Quad (3rd-gen, Spiral with glitch)
#mBD  = os.path.join(Bbase, 'NBDDIPOS.MAGLI.pdb')  # Main bending dipoles
mBD  = os.path.join(Bbase, 'NBDDIPOS.MAGLI.pkl')  # Main bending dipoles
#mBDh = os.path.join(Bbase, 'CRDIPOH.MAGLI.pdb')   # Main bending dipoles (Hi-Res)
mBDh = os.path.join(Bbase, 'CRDIPOH.MAGLI.pkl')   # Main bending dipoles (Hi-Res)
#mD5  = os.path.join(Bbase, 'SD5DIPOS.MAGLI.pdb')  # SD5 - regular dipole over straight pipe
mD5  = os.path.join(Bbase, 'SD5DIPOS.MAGLI.pkl')  # SD5 - regular dipole over straight pipe

# --- Y Magnets (injection/recirculation)
mYQi = os.path.join(Bbase, 'YQi.MAGLI.pdb')    # Y-Quad (tilted)
mPDi = os.path.join(Bbase, 'NPDi.MAGLI.pdb')   # Pulsed dipole
mQ1i = os.path.join(Bbase, 'RQi.MAGLI.pdb')    # QR1

mYQr = os.path.join(Bbase, 'YQr.MAGLI.pdb')   # Y-Quad (tilted)
mPDr = os.path.join(Bbase, 'NPDr.MAGLI.pdb')  # Pulsed dipole
mQ1r = os.path.join(Bbase, 'RQr.MAGLI.pdb')   # QR1

# --- Old Magnet Models (Wrong Transform)
mBDo = os.path.join(Bbase, 'CRDIPO.MAGLI.pdb')  # Main bending dipoles (wrongtransform-loop)
mBDs = os.path.join(Bbase, 'CRDIPOS.MAGLI.pdb') # Main bending dipoles (wrongtransform-Spiral)

#mYQio = os.path.join(Bbase, 'NYQ_A10_D8.pdb')    # Y-Quad (tilted)
mYQio = os.path.join(Bbase, 'NYQ_A10_D8.pkl')    # Y-Quad (tilted)
#mPDio = os.path.join(Bbase, 'NPD_A10_D8.pdb')    # Pulsed dipole
mPDio = os.path.join(Bbase, 'NPD_A10_D8.pkl')    # Pulsed dipole
#mQ1io = os.path.join(Bbase, 'NRQ_A10_D8.pdb')    # QR1
mQ1io = os.path.join(Bbase, 'NRQ_A10_D8.pkl')    # QR1

#mYQro = os.path.join(Bbase, 'NYQr_A10_D8.pdb')   # Y-Quad (tilted)
mYQro = os.path.join(Bbase, 'NYQr_A10_D8.pkl')   # Y-Quad (tilted)
#mPDro = os.path.join(Bbase, 'NPDr_A10_D8.pdb')   # Pulsed dipole
mPDro = os.path.join(Bbase, 'NPDr_A10_D8.pkl')   # Pulsed dipole
#mQ1ro = os.path.join(Bbase, 'NRQr_A10_D8.pdb')   # QR1
mQ1ro = os.path.join(Bbase, 'NRQr_A10_D8.pkl')   # QR1

# --- Steering Dipoles: Not yet implemented
#mD4 = os.path.join(Bbase, 'DIPO4.MAGLI.pdb')    # 4.5"  -flange steering dipoles
mD4 = os.path.join(Bbase, 'DIPO4.MAGLI.pkl')    # 4.5"  -flange steering dipoles
#mD3 = os.path.join(Bbase, 'DIPO3.MAGLI.pdb')    # 3.375"-flange steering dipoles
mD3 = os.path.join(Bbase, 'DIPO3.MAGLI.pkl')    # 3.375"-flange steering dipoles
#mSD = os.path.join(Bbase, 'SDIPO.MAGLI.pdb')    # Short steering dipoles
mSD = os.path.join(Bbase, 'SDIPO.MAGLI.pkl')    # Short steering dipoles



