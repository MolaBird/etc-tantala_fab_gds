import gdstk
import numpy as np


def _taperHolePeriod(d):
    dmax = HoleTaper_Period
    return 1 - dmax*np.abs(2*d**3 - 3*d**2 + 1)


def _taperHoleRadius_a(d):
    dmax = HoleTaper_Ra
    return 1 - dmax*np.abs(2*d**3 - 3*d**2 + 1)
    

def _taperHoleRadius_b(d):
    dmax = HoleTaper_Rb
    return 1 - dmax*np.abs(2*d**3 - 3*d**2 + 1)


def initHoleParams():
    global Tapering_Period, Tapering_Ra, Tapering_Rb
    global CavityRange, MirrorRange, CavitySize, MirrorSize
    global HolePosArray, HoleRaArray, HoleRbArray

    ii_cavityHole = np.arange(1, HoleNum_Cavity+1)
    pos_normalized = ii_cavityHole / HoleNum_Cavity

    Tapering_Period = _taperHolePeriod(pos_normalized)
    tapered_period = Tapering_Period * HolePeriod

    pos_array = np.array([sum(tapered_period[:i+1]) for i in range(len(tapered_period))]) - tapered_period[0]/2

    Tapering_Ra = _taperHoleRadius_a(pos_normalized)
    ra_array = np.ones_like(pos_array) * HoleRadius_a * Tapering_Ra

    Tapering_Rb = _taperHoleRadius_b(pos_normalized)
    rb_array = np.ones_like(pos_array) * HoleRadius_b * Tapering_Rb

    CavityRange = (0, pos_array.max() + 0.5*HolePeriod)
    CavitySize = 2*(CavityRange[1] - CavityRange[0])

    ii_mirrorHole = np.arange(1,HoleNum_Mirror+1)
    pos_array = np.append(
        pos_array, 
        pos_array.max() + ii_mirrorHole*HolePeriod
    )

    ra_array = np.append(
        ra_array, 
        np.ones_like(ii_mirrorHole)*HoleRadius_a
    )

    rb_array = np.append(
        rb_array, 
        np.ones_like(ii_mirrorHole)*HoleRadius_b
    )

    MirrorRange = (CavityRange[1], CavityRange[1] + HoleNum_Mirror*HolePeriod)
    MirrorSize = 2*(MirrorRange[1] - MirrorRange[0])

    HolePosArray, HoleRaArray, HoleRbArray = pos_array, ra_array, rb_array
    

def genNanobeam():
    nanobeam = gdstk.FlexPath(
        [(-Crystal_L/2, 0), (Crystal_L/2, 0)],
        width=Crystal_W,
        layer=1,
        datatype=0
    )
    return nanobeam


def genHoleArray():
    holes = []
    for pos, ra, rb in zip(HolePosArray, HoleRaArray, HoleRbArray):
        hole_r = gdstk.ellipse(
            (pos, 0), 
            (ra, rb),
            tolerance=1e-5,
            layer=1,
            datatype=0
        )
        hole_l = gdstk.ellipse(
            (-pos, 0), 
            (ra, rb),
            tolerance=1e-5,
            layer=1,
            datatype=0
        )
        holes.append(hole_l)
        holes.append(hole_r)

    return holes


def save2GDS(fp):
    from pathlib import Path
    fp = Path(fp)
    if not fp.parent.exists():
        fp.parent.mkdir(parents=True, exist_ok=True)

    gds_lib = gdstk.Library()
    gds_cell = gds_lib.new_cell("tantala")
    
    # holes = genHoleArray()
    structure = gdstk.boolean(
        operand1=genNanobeam(),
        operand2=genHoleArray(),
        operation="not",
        layer=1,
        datatype=0
    )

    gds_cell.add(*structure)
    gds_lib.write_gds(fp)





#======================================================================
ORIGIN = (0, 0, 0)

HolePeriod = 1.2294394              # untapered lattice constant
HoleRadius_a = 0.3428               # longitudinal radius of the hole
HoleRadius_b = 0.3428               # transverse radius of the hole
HoleNum_Cavity = 12                 # oneside
HoleNum_Mirror = 27                 # oneside
HoleNum_Couple = 0                  # oneside
HoleTaper_Period = 0.2074775552     #    
HoleTaper_Ra = 0.1647557            # 
HoleTaper_Rb = 0.1647557

HolePosArray = np.array([])
HoleRaArray = np.array([])
HoleRbArray =np.array([])

CavityRange = (0,0)
MirrorRange = (0,0)

initHoleParams()

#======================================================================
Crystal_L = 200
Crystal_W = 1.4739
Crystal_H = 0.8                     # um


fp = "tantala_design.gds"
# save2GDS(fp)

