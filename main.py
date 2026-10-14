

import gdstk
import numpy as np
LAYER_MAIN = 0
LAYER_UNDERCUT = 1
# LAYER_UNDERCUT = 1
# LAYER_UNDERCUT = 2
LAYER_AUX = 2


def _iter_units():
    for i in range(NCOL):
        for j in range(NROW):
            yield i, j


def _iter_devices():
    n = 0
    for i in range(NCOL):
        for j in range(NROW):
            for k in range(NANOBEAM_NUM_PER_UNIT):
                yield i, j, k, n
                n += 1


def init_gds_lib(lib_name="tantala_pcc_array", main_cell_name="MAIN"):
    global GDS_lib, GDS_cell
    GDS_lib = gdstk.Library(lib_name)
    GDS_cell = GDS_lib.new_cell(main_cell_name)
    return GDS_lib, GDS_cell


def load_design_config():
    import json 
    with open(DesignConfig_fp, "r") as f:
        return json.load(f)['designs']


def load_coupler(fp_coupler, arm_spacing, arm_width, arm_taper_radius, layer=LAYER_MAIN, datatype=0):
    r = arm_taper_radius
    def _create_arm_taper(origin, n_arc_points=64):
        x0, y0 = origin
        polygons = []
        for sign in [1, -1]:
            xc, yc = x0 + r, y0 + sign*(arm_width/2 + r)
            x1, y1 = x0, y0 + sign*arm_width/2
            points = [(x1, y1)]
            angles = np.linspace(0, sign*0.5*np.pi, n_arc_points + 1) + np.pi
            xx, yy = xc + r * np.cos(angles), yc + r * np.sin(angles)
            points.extend(zip(xx, yy))
            polygons.append(gdstk.Polygon(points, layer=layer, datatype=datatype))
        return polygons

    lib = gdstk.read_gds(fp_coupler)
    cell = lib.cells[0]
    coupler = [p.copy() for p in cell.polygons]
    coupler = []
    extents = []
    xmin, xmax, ymin, ymax = None, None, None, None
    for p in cell.polygons:
        p_copy = p.copy()
        p_copy.layer = layer
        p_copy.datatype = datatype
        coupler.append(p_copy)
        bbox = p.bounding_box()
        xmin = bbox[0][0] if xmin is None else min(xmin, bbox[0][0])
        xmax = bbox[1][0] if xmax is None else max(xmax, bbox[1][0])
        ymin = bbox[0][1] if ymin is None else min(ymin, bbox[0][1])
        ymax = bbox[1][1] if ymax is None else max(ymax, bbox[1][1])
    x_extent = xmax - xmin
    underfill = COUPLER_REGION_W - x_extent
    if underfill < 0:
        raise ValueError("Coupler x-extent exceeds the defined coupler region width.")
    gap_from_wall = COUPLER_REGION_OFFSET_FROM_WALL + underfill / 2
    gap_from_nanobeam = underfill / 2
    arms = [
        gdstk.rectangle(
            (xmin - gap_from_wall, arm_spacing/2),
            (xmin, arm_spacing/2+arm_width),
            layer=layer,
            datatype=datatype
        ),
        gdstk.rectangle(
            (xmin - gap_from_wall, -arm_spacing/2 - arm_width),
            (xmin, -arm_spacing/2),
            layer=layer,
            datatype=datatype
        )
    ]
    extents.append(
        gdstk.rectangle(
            (xmax, NANOBEAM_W/2),
            (xmax + gap_from_nanobeam, -NANOBEAM_W/2),
            layer=layer,
            datatype=datatype
        )
    )
    extents.extend(arms)
    extents.extend(_create_arm_taper((xmin - gap_from_wall, arm_spacing/2+arm_width/2)))
    extents.extend(_create_arm_taper((xmin - gap_from_wall, -arm_spacing/2-arm_width/2)))
    return Coupler(coupler=coupler, extents=extents, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)


def init_pos_of_units():
    global PosOfUnits
    step_x = WIDTH / NCOL
    step_y = HEIGHT / NROW

    for i, j in _iter_units():
        x = i * step_x + step_x / 2 - WIDTH/2
        y = j * step_y + step_y / 2 - HEIGHT/2
        PosOfUnits[i, j] = (x, y)


def init_pos_of_nanobeams():
    global RelPosYOfDevices
    idx = np.arange(NANOBEAM_NUM_PER_UNIT)
    ytop = UNITREGION_H/2
    RelPosYOfDevices = ytop - (idx + 0.5) * NANOBEAM_SPACING


def init_pos_of_couplers():
    global RelPosXOfCouplers
    w = COUPLER_REGION_W
    o = COUPLER_REGION_OFFSET_FROM_WALL
    g = SIDEWALL_GAP
    xl = -g/2 + o + w/2
    xr = g/2 - o - w/2
    RelPosXOfCouplers = (xl, xr)


def AuxLine_mark_origin(layer=LAYER_AUX, datatype=0):
    circle = gdstk.ellipse(
        (0, 0),
        20,
        tolerance=TOLERANCE,
        layer=layer,
        datatype=datatype
    )
    GDS_cell.add(circle)


def AuxLine_mark_main_region(layer=LAYER_AUX, datatype=0):
    rect = gdstk.rectangle(
        (-WIDTH/2, -HEIGHT/2),
        (WIDTH/2, HEIGHT/2),
        layer=layer,
        datatype=datatype
    )
    GDS_cell.add(rect)


def AuxLine_draw_main_grid(layer=LAYER_AUX, datatype=0, line_width=0.1):
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


def AuxLine_mark_unit_regions(layer=LAYER_AUX, datatype=0):
    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        rect = gdstk.rectangle(
            (xc-UNITREGION_W/2, yc-UNITREGION_H/2), 
            (xc+UNITREGION_W/2, yc + UNITREGION_H/2), 
            layer=layer, datatype=datatype
        )
        GDS_cell.add(rect)
    return GDS_cell


def AuxLine_mark_coupler_regions(layer=LAYER_AUX, datatype=0):
    w = COUPLER_REGION_W
    h = COUPLER_REGION_H
    for i, j, k, _ in _iter_devices():
        xc, yc = PosOfUnits[i, j]
        y = RelPosYOfDevices[k] + yc
        xl = xc + RelPosXOfCouplers[0]
        xr = xc + RelPosXOfCouplers[1]
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


def Main_place_sidewalls(layer=LAYER_MAIN, datatype=0):
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


def Main_gen_nanobeams(layer=LAYER_MAIN, datatype=0):
    global GDS_devices, GDS_devices_mask
    l = NANOBEAM_L
    w = NANOBEAM_W
    def place_1_nanobeam(i, j, k):
        xc, yc = PosOfUnits[i, j]
        xbeam = xc
        ybeam = RelPosYOfDevices[k] + yc
        nanobeam = gdstk.rectangle(
            (xbeam - l/2, ybeam - w/2),
            (xbeam + l/2, ybeam + w/2),
            layer=layer, datatype=datatype
        )
        # GDS_devices.append([nanobeam])
        pad = UNDERCUT_PAD
        mask = gdstk.rectangle(
            (xbeam - l/2 - pad, ybeam - w/2 - pad),
            (xbeam + l/2 + pad, ybeam + w/2 + pad),
            layer=LAYER_UNDERCUT, datatype=datatype
        )
        return nanobeam, mask

    for i, j, k, _ in _iter_devices():
        nanobeam, mask = place_1_nanobeam(i, j, k)
        GDS_devices.append([nanobeam])
        GDS_devices_mask.append([mask])
        # GDS_devices_mask.append([])




def _place_1_device(design, scaling, single_side_mirror_hole_num, xc, yc, layer, datatype=0):
    def _tapering(distance, taper_factor):
        dmax = taper_factor  
        d = distance 
        taper = 1 - dmax*np.abs(2*d**3 - 3*d**2 + 1)
        return taper
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
        arr_pos.max() + arr_mirrorhole_idx * period
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


def Main_gen_holes(layer=LAYER_MAIN, datatype=0):
    global GDS_devices
    for i, j, k, n in _iter_devices():
        if k == 0: continue

        xc, yc = PosOfUnits[i, j]
        design = i // UNITREPLICA_NUM
        scaling = SCALING_LIST[j]
        xbeam = xc
        ybeam = RelPosYOfDevices[k] + yc
        GDS_holes = _place_1_device(
            DesignConfig[design],
            scaling,
            NUM_MIRROR_HOLE_LIST[k-1],
            xbeam,
            ybeam,
            layer,
            datatype
        )
        GDS_devices[n] = gdstk.boolean(
            GDS_devices[n],
            GDS_holes,
            "not",
            precision=TOLERANCE,
            layer=layer
        )
    # GDS_devices.extend(temp_GDS_devices)


def Main_place_couplers(layer=LAYER_MAIN, datatype=0):
    global GDS_devices, GDS_devices_mask
    # temp_GDS_devices = []
    # temp_GDS_devices_mask = []
    for i, j, k, n in _iter_devices():
        if i % UNITREPLICA_NUM == 0:
            coupler = Coupler_1
        elif i % UNITREPLICA_NUM == 1:
            coupler = Coupler_2
        else:
            continue
        xc, yc = PosOfUnits[i, j]
        y = yc + RelPosYOfDevices[k]
        xl = xc + RelPosXOfCouplers[0]
        xr = xc + RelPosXOfCouplers[1]
        polys = []
        for p_ in coupler.coupler + coupler.extents:
            pl = p_.copy()
            pl.layer, pl.datatype = layer, datatype
            pl.translate(xl, y)
            pr = p_.copy()
            pr.layer, pr.datatype = layer, datatype
            pr.mirror((0, 0), (0, 1)).translate(xr, y)
            polys.extend([pl, pr])

        GDS_devices[n].extend(polys)

        pad = UNDERCUT_PAD
        # mask = GDS_devices_mask.pop(0)
        GDS_devices_mask[n].append(
            gdstk.rectangle(
                (xl + coupler.xmin - pad, y + coupler.ymin - pad),
                (xl + coupler.xmax + pad, y + coupler.ymax + pad),
                layer=LAYER_UNDERCUT,
                datatype=datatype
            )
        )
        GDS_devices_mask[n].append(
            gdstk.rectangle(
                (xr + coupler.xmin - pad, y + coupler.ymin - pad),
                (xr + coupler.xmax + pad, y + coupler.ymax + pad),
                layer=LAYER_UNDERCUT,
                datatype=datatype
            )
        )
        # temp_GDS_devices_mask.append(mask)
        # temp_GDS_devices.append(
        #     gdstk.boolean(
        #         device,
        #         polys,
        #         "or",
        #         precision=TOLERANCE,
        #         layer=layer
        #     )
        # )
    # GDS_devices.extend(temp_GDS_devices)
    # GDS_devices_mask.extend(temp_GDS_devices_mask)


def Main_place_ids(layer=LAYER_MAIN, datatype=0):
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
            f"S{scaling:.3f}-D{design+1}-U{unit_id+1}",
            size=ID_TEXT_SIZE,
            position=(xc, yc + UNITREGION_H/2 + ID_TEXT_OFFSET_Y), # align to center
            layer=layer,
            datatype=datatype
        )
        id_text = _center_x_text(id_text, xc)
        GDS_cell.add(*id_text)


def Main_place_positioners(layer=LAYER_MAIN, datatype=0):
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


def Main_place_undercuts(layer=LAYER_UNDERCUT, datatype=0):
    global GDS_devices_mask
    for i, j in _iter_units():
        xc, yc = PosOfUnits[i, j]
        pad = UNDERCUT_PAD
        pos_masks = gdstk.rectangle(
            (xc - SIDEWALL_GAP/2 + pad, yc - UNITREGION_H/2),
            (xc + SIDEWALL_GAP/2 - pad, yc + UNITREGION_H/2),
            layer=layer,
            datatype=datatype
        )
        neg_masks = []
        for _ in range(NANOBEAM_NUM_PER_UNIT):
            neg_masks.extend(GDS_devices_mask.pop(0))
        GDS_cell.add(
            *gdstk.boolean(
                pos_masks,
                neg_masks,
                "not",
                precision=TOLERANCE,
                layer=layer
            )
        )



def export_gds(lib, gds_filename="output.gds"):
    lib.write_gds(gds_filename)
    print(f"GDS successfully saved to: {gds_filename}")


def gen_individual_gds(folder, layer=LAYER_MAIN, datatype=0):
    from pathlib import Path
    folder_path = Path(folder)
    # if folder_path.exists():
    #     folder_path.rmdir()
    folder_path.mkdir(parents=True, exist_ok=True)

    scaling = 1.0
    xc, yc = 0, 0
    nanobeam = gdstk.rectangle(
        (-NANOBEAM_L/2, -NANOBEAM_W/2),
        (NANOBEAM_L/2, NANOBEAM_W/2),
        layer=layer, datatype=datatype
    )
    for i, design in enumerate(DesignConfig):
        for mirror_hole_num in NUM_MIRROR_HOLE_LIST:
            gds_lib = gdstk.Library(str(folder_path/f'S1-D{i+1}-M{mirror_hole_num}'))
            gds_cell = gds_lib.new_cell('MAIN')
            holes = _place_1_device(design, scaling, mirror_hole_num, xc, yc, layer=layer, datatype=datatype)
            beam = nanobeam.copy()
            gds_cell.add(
                *gdstk.boolean(
                    beam,
                    holes,
                    "not",
                    precision=TOLERANCE,
                    layer=layer
                )
            )
            gds_lib.write_gds(folder_path/f'S1-D{i+1}-M{mirror_hole_num}.gds')
    

TOLERANCE = 1e-5

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
COUPLER_REGION_OFFSET_FROM_NANOBEAM = SIDEWALL_GAP/2 - COUPLER_REGION_W - COUPLER_REGION_OFFSET_FROM_WALL

UNDERCUT_PAD = 2

POSITIONER_RADIUS = 2
POSITIONER_OFFSET_Y = 15

ID_TEXT_SIZE = 12
ID_TEXT_OFFSET_Y =  30

SCALING_LIST = np.arange(-10, 11, 1) * 0.005 + 1
NUM_MIRROR_HOLE_LIST = (0, 2, 4, 6, 8)


DesignConfig_fp = "./tantala_pcc_designs.json"
DesignConfig = load_design_config()
# print(DesignConfig)
PosOfUnits = np.full((NROW, NCOL, 2), np.nan)
RelPosYOfDevices = np.full(NANOBEAM_NUM_PER_UNIT, np.nan)
RelPosXOfCouplers = ()
# RelPosYOfUndercuts = np.full(NANOBEAM_NUM_PER_UNIT, np.nan)

GDS_lib = None
GDS_cell = None

GDS_devices = []
GDS_devices_mask = []

fp_coupler_1 = "./coupler_old.gds"
fp_coupler_2 = "./coupler_new_Austin.gds"
# fp_coupler_3 = "./coupler_new_Henry.gds"

from collections import namedtuple
Coupler = namedtuple("Coupler", ["coupler", "extents", "xmin", "xmax", "ymin", "ymax"])
Coupler_1 :Coupler= load_coupler(fp_coupler_1, arm_spacing=2, arm_width=0.5, arm_taper_radius=0.8)
Coupler_2 :Coupler= load_coupler(fp_coupler_2, arm_spacing=2, arm_width=1.0, arm_taper_radius=1)
# GDS_coupler_3 = load_coupler(fp_coupler_3)

if __name__ == "__main__":
    init_pos_of_units()
    init_pos_of_nanobeams()
    init_pos_of_couplers()

    # gen_individual_gds("individual_designs")
    # import sys
    # sys.exit(0)

    init_gds_lib()

    # AuxLine_draw_main_grid()
    # AuxLine_mark_main_region()
    # AuxLine_mark_origin()
    # AuxLine_mark_unit_regions()
    # AuxLine_mark_coupler_regions()

    Main_place_sidewalls()
    Main_gen_nanobeams()
    Main_gen_holes()
    Main_place_ids()
    Main_place_positioners()
    Main_place_couplers()
    Main_place_undercuts()


    for dev in GDS_devices:
        GDS_cell.add(*dev)

    export_gds(GDS_lib, gds_filename="cavity_array.gds")



