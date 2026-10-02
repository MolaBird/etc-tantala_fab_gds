

import gdstk
import numpy as np

def init_gds_lib(lib_name="cavity_library", main_cell_name="MAIN"):
    """
    Initializes a new GDS library and creates the top-level cell.
    
    Args:
        lib_name (str): Name of the GDS library.
        main_cell_name (str): Name of the primary cell.
        
    Returns:
        tuple: (gdstk.Library, gdstk.Cell)
    """
    lib = gdstk.Library(lib_name)
    main_cell = lib.new_cell(main_cell_name)
    return lib, main_cell


def load_design_config():
    import json 
    with open(DesignConfig_fp, "r") as f:
        return json.load(f)['designs']


def init_pos_of_units():
    global PosOfUnits
    step_x = WIDTH / NCOL
    step_y = HEIGHT / NROW

    for i in range(NCOL):
        for j in range(NROW):
            x = i * step_x + step_x / 2
            y = j * step_y + step_y / 2
            PosOfUnits[j, i] = (x+step_x/2, y+step_y/2)


def init_pos_of_nanobeams():
    global RelPosYOfNanobeams
    idx = np.arange(NANOBEAM_NUM_PER_UNIT)
    ytop = UNITREGION_H/2
    RelPosYOfNanobeams = ytop - (idx + 0.5) * NANOBEAM_SPACING


def AuxLine_draw_main_grid(cell, layer, datatype=0, line_width=0.1):
    """
    Draws a grid pattern on a specified layer using FlexPaths.
    
    Args:
        cell (gdstk.Cell): The cell to add the grid to.
        layer (int): GDS layer number for the grid.
        datatype (int): GDS datatype (default 0).
        line_width (float): Width of the grid lines.
    """
    step_x = WIDTH / NCOL
    step_y = HEIGHT / NROW
    for i in range(NCOL):
        for j in range(NROW):
            x, y = PosOfUnits[j, i]
            x -= step_x / 2
            y -= step_y / 2
            v_line = gdstk.FlexPath([(x, 0), (x, HEIGHT)], line_width, layer=layer, datatype=datatype)
            cell.add(v_line)
            h_line = gdstk.FlexPath([(0, y), (WIDTH, y)], line_width, layer=layer, datatype=datatype)
            cell.add(h_line)
    return cell


def AuxLine_mark_unit_regions(cell, layer, datatype=0):
    """
    Marks the region of each unit in the grid.
    
    Args:
        cell (gdstk.Cell): The cell to add the unit regions to.
        layer (int): GDS layer number for the unit regions.
        datatype (int): GDS datatype (default 0).
        line_width (float): Width of the unit region borders.
    """
    for i in range(NCOL):
        for j in range(NROW):
            xc, yc = PosOfUnits[j, i]
            rect = gdstk.rectangle(
                (xc-UNITREGION_W/2, yc-UNITREGION_H/2), 
                (xc+UNITREGION_W/2, yc + UNITREGION_H/2), 
                layer=layer, datatype=datatype
            )
            cell.add(rect)
    return cell


def AuxLine_mark_coupler_regions(cell, layer, datatype=0):
    w = COUPLER_REGION_W
    h = COUPLER_REGION_H
    o = COUPLER_REGION_OFFSET_FROM_WALL
    d = NANOBEAM_SPACING
    def _mark_1_pair_of_coupler_regions(i_col, j_row, n_idx):
        xc, yc = PosOfUnits[j_row, i_col]
        x0 = xc - SIDEWALL_GAP/2 + o
        x1 = xc + SIDEWALL_GAP/2 - o
        y = RelPosYOfNanobeams[n_idx] + yc
        rect = gdstk.rectangle(
            (x0, y-h/2),
            (x0+w, y+h/2),
            layer=layer, datatype=datatype
        )
        cell.add(rect)
        rect = gdstk.rectangle(
            (x1, y-h/2),
            (x1-w, y+h/2),
            layer=layer, datatype=datatype
        )
        cell.add(rect)  
    for i in range(NCOL):
        for j in range(NROW):
            for k in range(NANOBEAM_NUM_PER_UNIT):
                _mark_1_pair_of_coupler_regions(i, j, k)


def Main_place_sidewalls(layer, datatype=0):
    global GDS_sidewalls
    for i in range(NROW):
        for j in range(NCOL):
            xc, yc = PosOfUnits[i, j]
            coord1 = (xc-UNITREGION_W/2, yc-UNITREGION_H/2)
            coord2 = (xc-UNITREGION_W/2+SIDEWALL_W, yc+UNITREGION_H/2)
            sidewall_left = gdstk.rectangle(
                coord1, 
                coord2, 
                layer=layer, datatype=datatype
            )
            sidewall_left.fillet(SIDEWALL_CORNER_RAD, tolerance=TOLERANCE)
            GDS_sidewalls.append(sidewall_left)
            coord1 = (xc+UNITREGION_W/2-SIDEWALL_W, yc-UNITREGION_H/2)
            coord2 = (xc+UNITREGION_W/2, yc+UNITREGION_H/2)
            sidewall_right = gdstk.rectangle(
                coord1, 
                coord2, 
                layer=layer, datatype=datatype
            )
            sidewall_right.fillet(SIDEWALL_CORNER_RAD, tolerance=TOLERANCE)
            GDS_sidewalls.append(sidewall_right)


def Main_place_nanobeams(layer, datatype=0):
    global GDS_nanobeams
    d = NANOBEAM_SPACING
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

    for i in range(NROW):
        for j in range(NCOL):
            for k in range(NANOBEAM_NUM_PER_UNIT):
                place_1_nanobeam(i, j, k)


def Main_place_holes(layer, datatype=0):
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

    for i in range(NROW):
        for j in range(NCOL):
            xc, yc = PosOfUnits[i, j]
            design = j % len(DesignConfig)
            scaling = SCALING_LIST[i]
            for k in range(NANOBEAM_NUM_PER_UNIT):
                if k == 0:
                    GDS_pcc.append([GDS_nanobeams.pop(0)])
                    continue
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


def save_and_export_preview(lib, gds_filename="output.gds", svg_filename="preview.svg"):
    # Save the standard GDSII file
    lib.write_gds(gds_filename)
    print(f"GDS successfully saved to: {gds_filename}")
    
    # Export an SVG representation of the top-level cell for quick previewing
    # cell.write_svg(svg_filename)
    # print(f"Preview SVG successfully saved to: {svg_filename}")



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

COUPLER_REGION_W = 9
COUPLER_REGION_H = 9
COUPLER_REGION_OFFSET_FROM_WALL = 2


SCALING_LIST = np.arange(-10, 10, 1) * 0.01 + 1
NUM_MIRROR_HOLE_LIST = (0, 2, 4, 6, 8)


DesignConfig_fp = "./tantala_pcc_designs.json"
DesignConfig = load_design_config()

PosOfUnits = np.full((NROW, NCOL, 2), np.nan)
RelPosYOfNanobeams = np.full(NANOBEAM_NUM_PER_UNIT, np.nan)

GDS_nanobeams = []
GDS_pcc = []
GDS_sidewalls = []


if __name__ == "__main__":
    my_lib, my_cell = init_gds_lib(lib_name="PhC_Arrays", main_cell_name="OP")
    init_pos_of_units()
    init_pos_of_nanobeams()

    # AuxLine_mark_unit_regions(my_cell, layer=1)
    AuxLine_mark_coupler_regions(my_cell, layer=1)

    Main_place_sidewalls(layer=2)
    Main_place_nanobeams(layer=2)
    Main_place_holes(layer=2)

    my_cell.add(*GDS_sidewalls)
    for pcc in GDS_pcc:
        my_cell.add(*pcc)

    save_and_export_preview(my_lib, gds_filename="cavity_array.gds", svg_filename="cavity_array.svg")



