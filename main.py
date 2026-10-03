

import gdstk
import numpy as np

def _iter_units():
    for i in range(NCOL):
        for j in range(NROW):
            yield i, j

def _iter_nanobeams():
    for i in range(NCOL):
        for j in range(NROW):
            for k in range(NANOBEAM_NUM_PER_UNIT):
                yield i, j, k


def init_gds_lib(lib_name="tantala_pcc_array", main_cell_name="MAIN"):
    global GDS_lib, GDS_cell
    GDS_lib = gdstk.Library(lib_name)
    GDS_cell = GDS_lib.new_cell(main_cell_name)
    return GDS_lib, GDS_cell


def load_design_config():
    import json 
    with open(DesignConfig_fp, "r") as f:
        return json.load(f)['designs']


def init_pos_of_units():
    global PosOfUnits
    step_x = WIDTH / NCOL
    step_y = HEIGHT / NROW

    for i, j in _iter_units():
        x = i * step_x + step_x / 2 - WIDTH/2
        y = j * step_y + step_y / 2 - HEIGHT/2
        PosOfUnits[i, j] = (x, y)


def init_pos_of_nanobeams():
    global RelPosYOfNanobeams
    idx = np.arange(NANOBEAM_NUM_PER_UNIT)
    ytop = UNITREGION_H/2
    RelPosYOfNanobeams = ytop - (idx + 0.5) * NANOBEAM_SPACING


def init_pos_of_couplers():
    global RelPosXOfCouplers
    w = COUPLER_REGION_W
    o = COUPLER_REGION_OFFSET_FROM_WALL
    g = SIDEWALL_GAP
    xl = -g/2 + o + w/2
    xr = g/2 - o - w/2
    RelPosXOfCouplers = (xl, xr)





def AuxLine_mark_origin(layer, datatype=0):
    circle = gdstk.ellipse(
        (0, 0),
        20,
        tolerance=TOLERANCE,
        layer=layer,
        datatype=datatype
    )
    GDS_cell.add(circle)


def AuxLine_mark_main_region(layer, datatype=0):
    rect = gdstk.rectangle(
        (-WIDTH/2, -HEIGHT/2),
        (WIDTH/2, HEIGHT/2),
        layer=layer,
        datatype=datatype
    )
    GDS_cell.add(rect)


def AuxLine_draw_main_grid(layer, datatype=0, line_width=0.1):
    step_x = WIDTH / NCOL
    step_y = HEIGHT / NROW
    for i, j in _iter_units():
        x, y = PosOfUnits[i, j]
        x -= step_x / 2
        y -= step_y / 2
        v_line = gdstk.FlexPath([(x, -HEIGHT/2), (x, HEIGHT/2)], line_width, layer=layer, datatype=datatype)
        GDS_cell.add(v_line)
        h_line = gdstk.FlexPath([(-WIDTH/2, y), (WIDTH/2, y)], line_width, layer=layer, datatype=datatype)
        GDS_cell.add(h_line)
    x += step_x
    y += step_y
    v_line = gdstk.FlexPath([(x, -HEIGHT/2), (x, HEIGHT/2)], line_width, layer=layer, datatype=datatype)
    GDS_cell.add(v_line)
    h_line = gdstk.FlexPath([(-WIDTH/2, y), (WIDTH/2, y)], line_width, layer=layer, datatype=datatype)
    GDS_cell.add(h_line)
    return GDS_cell


def AuxLine_mark_unit_regions(layer, datatype=0):
    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        rect = gdstk.rectangle(
            (xc-UNITREGION_W/2, yc-UNITREGION_H/2), 
            (xc+UNITREGION_W/2, yc + UNITREGION_H/2), 
            layer=layer, datatype=datatype
        )
        GDS_cell.add(rect)
    return GDS_cell


def AuxLine_mark_coupler_regions(layer, datatype=0):
    w = COUPLER_REGION_W
    h = COUPLER_REGION_H
    def _mark_1_pair_of_coupler_regions(i, j, k):
        xc, yc = PosOfUnits[i, j]
        xl = xc + RelPosXOfCouplers[0]
        xr = xc + RelPosXOfCouplers[1]
        y = RelPosYOfNanobeams[k] + yc
        rect = gdstk.rectangle(
            (xl-w/2, y-h/2),
            (xl+w/2, y+h/2),
            layer=layer, datatype=datatype
        )
        GDS_cell.add(rect)
        rect = gdstk.rectangle(
            (xr-w/2, y-h/2),
            (xr+w/2, y+h/2),
            layer=layer, datatype=datatype
        )
        GDS_cell.add(rect)  
    for i, j, k in _iter_nanobeams():
        _mark_1_pair_of_coupler_regions(i, j, k)


def Main_place_sidewalls(layer, datatype=0):
    gds_sidewalls = []
    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        coord1 = (xc-UNITREGION_W/2, yc-UNITREGION_H/2)
        coord2 = (xc-UNITREGION_W/2+SIDEWALL_W, yc+UNITREGION_H/2)
        sidewall_left = gdstk.rectangle(
            coord1, 
            coord2, 
            layer=layer, datatype=datatype
        )
        sidewall_left.fillet(SIDEWALL_CORNER_RAD, tolerance=TOLERANCE)
        gds_sidewalls.append(sidewall_left)
        coord1 = (xc+UNITREGION_W/2-SIDEWALL_W, yc-UNITREGION_H/2)
        coord2 = (xc+UNITREGION_W/2, yc+UNITREGION_H/2)
        sidewall_right = gdstk.rectangle(
            coord1, 
            coord2, 
            layer=layer, datatype=datatype
        )
        sidewall_right.fillet(SIDEWALL_CORNER_RAD, tolerance=TOLERANCE)
        gds_sidewalls.append(sidewall_right)
    GDS_cell.add(*gds_sidewalls)


def Main_gen_nanobeams(layer, datatype=0):
    global GDS_nanobeams
    l = NANOBEAM_L
    w = NANOBEAM_W
    def place_1_nanobeam(i, j, k):
        xc, yc = PosOfUnits[i, j]
        xbeam = xc
        ybeam = RelPosYOfNanobeams[k] + yc
        nanobeam = gdstk.rectangle(
            (xbeam - l/2, ybeam - w/2),
            (xbeam + l/2, ybeam + w/2),
            layer=layer, datatype=datatype
        )
        GDS_nanobeams.append(nanobeam)

    for i, j, k in _iter_nanobeams():
        place_1_nanobeam(i, j, k)


def Main_gen_holes(layer, datatype=0):
    global GDS_pcc
    def _tapering(distance, taper_factor):
        dmax = taper_factor  
        d = distance 
        taper = 1 - dmax*np.abs(2*d**3 - 3*d**2 + 1)
        return taper

    def _place_1_device(design, scaling, single_side_mirror_hole_num, xc, yc):
        rx = design['hole_rad_x'] * scaling
        ry = design['hole_rad_y'] * scaling
        period = design['period'] * scaling
        n_cav_hole = design['single_side_cavity_hole_num']
        arr_cavityhole_idx = np.arange(1, n_cav_hole+1)
        arr_normalized_pos = arr_cavityhole_idx / n_cav_hole
        arr_period = period * _tapering(arr_normalized_pos, design['period_taper'])

        arr_pos = np.array([sum(arr_period[:i+1]) for i in range(len(arr_period))]) - arr_period[0]/2
        arr_rx = rx * _tapering(arr_normalized_pos, design['hole_rad_x_taper'])
        arr_ry = ry * _tapering(arr_normalized_pos, design['hole_rad_y_taper'])

        arr_mirrorhole_idx = np.arange(1, single_side_mirror_hole_num+1)
        arr_pos = np.append(
            arr_pos,
            arr_pos.max() + arr_mirrorhole_idx * arr_period[0]
        )

        arr_rx = np.append(
            arr_rx,
            rx * np.ones_like(arr_mirrorhole_idx)
        )

        arr_ry = np.append(
            arr_ry,
            ry * np.ones_like(arr_mirrorhole_idx)
        )

        GDS_holes = []
        for pos, rx, ry in zip(arr_pos, arr_rx, arr_ry):
            hole_r = gdstk.ellipse(
                (pos+xc, yc),
                (rx, ry),
                tolerance=TOLERANCE,
                layer=layer,
                datatype=datatype
            )
            GDS_holes.append(hole_r)
            hole_l = gdstk.ellipse(
                (-pos+xc, yc),
                (rx, ry),
                tolerance=TOLERANCE,
                layer=layer,
                datatype=datatype
            )
            GDS_holes.append(hole_l)

        return GDS_holes

    for i, j, k in _iter_nanobeams():
        if k == 0:
            GDS_pcc.append([GDS_nanobeams.pop(0)])
            continue
        xc, yc = PosOfUnits[i, j]
        design = i // UNITREPLICA_NUM
        scaling = SCALING_LIST[j]
        xbeam = xc
        ybeam = RelPosYOfNanobeams[k] + yc
        GDS_holes = _place_1_device(
            DesignConfig[design],
            scaling,
            NUM_MIRROR_HOLE_LIST[k-1],
            xbeam,
            ybeam
        )
        GDS_pcc.append(
            gdstk.boolean(
                GDS_nanobeams.pop(0),
                GDS_holes,
                "not",
                layer=layer
            )
        )


def Main_place_pccs():
    for pcc in GDS_pcc:
        GDS_cell.add(*pcc)


def Main_place_ids(layer, datatype=0):
    def _center_x_text(polygons, xc):
        min_x = min(p.points[:, 0].min() for p in polygons)
        max_x = max(p.points[:, 0].max() for p in polygons)
        xc_ = (min_x + max_x) / 2.0        
        for p in polygons:
            p.translate(xc-xc_, 0)
        return polygons

    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        unit_id = i % UNITREPLICA_NUM
        design = i // UNITREPLICA_NUM
        scaling = SCALING_LIST[j]
        id_text = gdstk.text(
            f"S{scaling:.2f}-D{design+1}-U{unit_id+1}",
            size=ID_TEXT_SIZE,
            position=(xc, yc + UNITREGION_H/2 + ID_TEXT_OFFSET_Y), # align to center
            layer=layer,
            datatype=datatype
        )
        id_text = _center_x_text(id_text, xc)
        GDS_cell.add(*id_text)


def Main_place_positioners(layer, datatype=0):
    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        xl = xc + RelPosXOfCouplers[0]
        xr = xc + RelPosXOfCouplers[1]
        y = yc + UNITREGION_H/2 + POSITIONER_OFFSET_Y
        circ_l = gdstk.ellipse(
            (xl, y),
            POSITIONER_RADIUS,
            layer=layer,
            datatype=datatype
        )
        GDS_cell.add(circ_l)
        circ_r = gdstk.ellipse(
            (xr, y),
            POSITIONER_RADIUS,
            layer=layer,
            datatype=datatype
        )
        GDS_cell.add(circ_r)


# def 


def export_gds(lib, gds_filename="output.gds"):
    lib.write_gds(gds_filename)
    print(f"GDS successfully saved to: {gds_filename}")


TOLERANCE = 1e-4

WIDTH = 4000.0
HEIGHT = 6000.0
NROW = 20
NCOL = 20

NANOBEAM_L = 48
NANOBEAM_W = 1.7
NANOBEAM_SPACING = 30
NANOBEAM_NUM_PER_UNIT = 6

SIDEWALL_W = 20
SIDEWALL_H = NANOBEAM_SPACING * NANOBEAM_NUM_PER_UNIT
SIDEWALL_GAP = 70
SIDEWALL_CORNER_RAD = 0.8

UNITREGION_W = SIDEWALL_W + SIDEWALL_GAP + SIDEWALL_W
UNITREGION_H = SIDEWALL_H

UNITREPLICA_NUM = 4

COUPLER_REGION_W = 9
COUPLER_REGION_H = 9
COUPLER_REGION_OFFSET_FROM_WALL = 2

POSITIONER_RADIUS = 2
POSITIONER_OFFSET_Y = 10

ID_TEXT_SIZE = 12
ID_TEXT_OFFSET_Y =  20

SCALING_LIST = np.arange(-10, 10, 1) * 0.01 + 1
NUM_MIRROR_HOLE_LIST = (0, 2, 4, 6, 8)


DesignConfig_fp = "./tantala_pcc_designs.json"
DesignConfig = load_design_config()

PosOfUnits = np.full((NROW, NCOL, 2), np.nan)
RelPosYOfNanobeams = np.full(NANOBEAM_NUM_PER_UNIT, np.nan)
RelPosXOfCouplers = ()

GDS_lib = None
GDS_cell = None

GDS_nanobeams = []
GDS_pcc = []


if __name__ == "__main__":
    init_gds_lib()
    init_pos_of_units()
    init_pos_of_nanobeams()
    init_pos_of_couplers()

    AuxLine_draw_main_grid(layer=1)
    # AuxLine_mark_main_region(layer=1)
    AuxLine_mark_origin(layer=1)
    # AuxLine_mark_unit_regions(layer=1)
    AuxLine_mark_coupler_regions(layer=1)

    Main_place_sidewalls(layer=2)
    Main_gen_nanobeams(layer=2)
    Main_gen_holes(layer=2)
    Main_place_ids(layer=2)
    Main_place_positioners(layer=2)


    # GDS_cell.add(*GDS_nanobeams)
    for pcc in GDS_pcc:
        GDS_cell.add(*pcc)

    export_gds(GDS_lib, gds_filename="cavity_array.gds")



