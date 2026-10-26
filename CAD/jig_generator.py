# -*- coding: utf-8 -*-
"""Generate the light-tight FLIM measuring jig inside CAD/FLIMLab-meas-setup.FCStd.

Parts built here:
  jig-cup         cup that sits on the leaf: detection port (-X, 45 deg) made from the
                  filter-holder Body in the document, excitation tube (-Y, 45 deg), AMBIT cradle on top
  jig-lid         lid that closes the AMBIT cradle (4 x M3 into heat-set inserts)
  jig-sma-insert  printed SMA905 receptacle that slides in the excitation tube

Coordinates: origin = laser spot on the leaf, Z up, mm. Detection axis points to -X/+Z,
excitation axis to -Y/+Z, the AMBIT stick lies along Y with its head toward -Y.

Run from the repository root, with the file closed in the FreeCAD GUI:
    "%LOCALAPPDATA%\\Programs\\FreeCAD 1.1\\bin\\freecadcmd.exe" CAD\\jig_generator.py

Environment overrides:
    AMBIT_STEP       AMBIT casing STEP (A280.0101.000 ... Assembly_Rev06.STEP)
    JIG_FCSTD        document to update (default CAD/FLIMLab-meas-setup.FCStd)
    JIG_EXPORT_DIR   STL/STEP output folder (default CAD/export)
    JIG_PREVIEW_DIR  if set, tessellations and sections are written there for preview plots
    JIG_SKIP_LIGHTCHECK=1  skip the slow light-tightness boolean check

Edit PARAMS below and rerun; previously generated objects are replaced, the holder Body is kept.
"""
import os
import sys
import math
import time
import json

import FreeCAD as App
import Part
import MeshPart
import Import
import numpy as np

V = App.Vector
R2 = math.sqrt(2.0)
HERE = os.path.dirname(os.path.abspath(__file__))

AMBIT_STEP = os.environ.get(
    "AMBIT_STEP",
    os.path.join(os.path.expanduser("~"), "Downloads",
                 "A280.0101.000_3d_JII Asens._Assembly_Rev06.STEP"))
FCSTD = os.environ.get("JIG_FCSTD", os.path.join(HERE, "FLIMLab-meas-setup.FCStd"))
EXPORT_DIR = os.environ.get("JIG_EXPORT_DIR", os.path.join(HERE, "export"))
PREVIEW_DIR = os.environ.get("JIG_PREVIEW_DIR", "")
SKIP_LIGHTCHECK = os.environ.get("JIG_SKIP_LIGHTCHECK", "") == "1"

PARAMS = dict(
    # ---- cup and cavity -------------------------------------------------------------
    cav_half=22.0,          # cavity half width at the leaf (opening 44 x 44)
    cav_ceiling=27.0,       # flat ceiling height inside the cavity
    roof_top=30.0,          # outer roof height
    box_half=25.0,          # outer half width (wall 3 mm)
    px_ceil_edge=12.0,      # +X inner 45 deg slope meets the ceiling at this x
    py_ceil_edge=18.0,      # +Y inner 45 deg slope meets the ceiling at this y
    flange_half=31.0,       # bottom flange half width (6 mm land for foam tape)
    flange_t=3.0,
    # ---- detection port (-X, 45 deg), plate = filter-holder Body ------------------------
    D_face=38.0,            # FiberPort mating face distance from the spot along the axis
    D_boss_R=21.5,          # boss radius behind the plate (holder outline is R 19.05)
    D_boss_start=18.0,      # boss starts this far from the spot (inside the wall)
    D_light_tunnel_d=7.0,   # light path from the cavity to the plate aperture
    D_head_tunnel_d=6.5,    # access for the FiberPort screw heads
    D_tunnel_start=10.0,
    filter_t=3.0,           # display only: filter thickness
    fiberport_len=26.0,     # display only: approximate FiberPort body length
    # ---- excitation port (-Y, 45 deg) -------------------------------------------------
    E_slope_in=24.0,        # inner -Y slope distance from the spot
    E_slope_out=27.0,       # outer -Y slope distance
    E_tube_R=10.0,
    E_tube_start=20.0,
    E_tube_end=51.0,
    E_bore_d=12.3,          # sliding fit for the insert
    E_stop=26.0,            # insert stop shoulder distance from the spot
    E_stop_d=10.0,          # beam aperture at the stop
    E_cone_half=15.0,       # opening cone from the stop into the cavity, deg
    E_setscrew_boss_h=5.0,  # half diagonal of the diamond boss for the grub-screw insert
    E_setscrew_boss_end=12.5,
    # ---- SMA905 insert ------------------------------------------------------------
    ins_d=12.0,
    ins_len=14.0,           # body length (inner face to stud shoulder)
    stud_d=6.2,             # the 1/4-36 coupling nut forms its own thread on this
    stud_len=6.5,
    ferrule_bore_d=3.3,     # SMA905 ferrule is 3.175
    ferrule_len=9.8,        # ferrule length from the connector shoulder (measure yours)
    ferrule_bore_extra=0.5,
    ins_cone_half=15.0,     # exit cone after the fibre tip, deg
    ins_min_engagement=10.0,
    fibre_NA=0.22,
    L_E_display=42.0,       # fibre tip distance used for the displayed assembly
    # ---- AMBIT cradle and lid ---------------------------------------------------------
    ambit_gap=30.0,         # leaf to lowest point of the AMBIT casing
    ambit_centre_stepz=3.0, # STEP z of the casing put over the spot (between LED and PD)
    cradle_half=12.5,
    cradle_y0=-10.0,
    cradle_y1=38.0,
    cradle_z0=27.0,
    parting_z=36.0,         # cradle/lid split, close to the AMBIT shell seam
    gusset_y=25.0,          # tongue starts at the +Y wall, 45 deg gusset below it
    clear_xy=0.3,           # lateral clearance around the casing
    clear_floor=0.1,        # under the casing
    clear_lid=0.2,          # above the casing
    lid_wall=2.0,
    hm_dx=0.15,             # height-map slice spacing (x)
    hm_dy=0.1,              # height-map bin size (y)
    hm_simplify=0.03,       # profile simplification tolerance (added to the clearance)
    hm_merge=0.05,          # strips whose profiles differ less than this are merged
    window_half_x=4.0,      # optical window under the AMBIT head
    window_y0=-6.0,
    window_y1=15.0,
    screw_y_head=-6.0,
    screw_y_tail=32.0,
    heatset_d=4.0,          # M3 x 5.7 heat-set insert hole
    heatset_depth=6.0,
    screw_clear_d=3.4,
    screw_wall_min=1.0,
)

UD = V(-1, 0, 1) * (1 / R2)   # detection axis: spot -> FiberPort
VD = V(1, 0, 1) * (1 / R2)    # up-slope direction in the detection plane
UE = V(0, -1, 1) * (1 / R2)   # excitation axis: spot -> fibre
X = V(1, 0, 0)

CASING_KEYS = ("Shell", "Cable Cover", "Mask", "Vent")
OPTIC_KEYS = ("PD15", "MLX90632", "AS7341", "TEMD7100", "OSLON", "FFC")

T0 = time.time()


def log(msg):
    print("[jig %6.1fs] %s" % (time.time() - T0, msg))
    sys.stdout.flush()


# ----------------------------------------------------------------------------- helpers
def halfspace_ge(n, d, L=600.0):
    """Big box covering {p : p.n >= d}."""
    n = V(n)
    n.normalize()
    b = Part.makeBox(L, L, L, V(-L / 2, -L / 2, 0))
    b.Placement = App.Placement(n * d, App.Rotation(V(0, 0, 1), n))
    return b


def single_solid(shape, name):
    sols = sorted(shape.Solids, key=lambda s: -s.Volume)
    if not sols:
        raise RuntimeError("%s: boolean produced no solid" % name)
    if len(sols) > 1:
        log("WARNING %s: %d solids, keeping the largest (others %s mm3)"
            % (name, len(sols), ", ".join("%.3f" % s.Volume for s in sols[1:6])))
    return sols[0]


def transformed(shape, mat):
    s = shape.copy()
    s.transformShape(mat)
    return s


def matrix_from_axes(cx, cy, cz, t):
    return App.Matrix(cx.x, cy.x, cz.x, t.x,
                      cx.y, cy.y, cz.y, t.y,
                      cx.z, cy.z, cz.z, t.z,
                      0, 0, 0, 1)


def holder_outline_prism(R, flat, z0, z1):
    """Holder outline (circle R cut by flats at |x| = flat) extruded from z0 to z1 (local)."""
    cyl = Part.makeCylinder(R, z1 - z0, V(0, 0, z0), V(0, 0, 1))
    box = Part.makeBox(2 * flat, 2 * R + 2, z1 - z0, V(-flat, -R - 1, z0))
    return cyl.common(box)


def slice_x(shape, x):
    for d in (0.0, 0.0131, -0.0173, 0.0291):
        try:
            return shape.slice(X, x + d)
        except Exception:
            continue
    log("WARNING slice at x=%.3f failed" % x)
    return []


def dp_simplify(pts, tol):
    pts = np.asarray(pts, float)
    if len(pts) < 3:
        return pts
    keep = np.zeros(len(pts), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        p, q = pts[a], pts[b]
        seg = q - p
        L = math.hypot(seg[0], seg[1])
        rel = pts[a + 1:b] - p
        if L > 0:
            d = np.abs(seg[0] * rel[:, 1] - seg[1] * rel[:, 0]) / L
        else:
            d = np.hypot(rel[:, 0], rel[:, 1])
        i = int(np.argmax(d))
        if d[i] > tol:
            m = a + 1 + i
            keep[m] = True
            stack += [(a, m), (m, b)]
    return pts[keep]


def common_volume(a, b):
    try:
        return a.common(b).Volume
    except Exception as e:
        log("common() failed: %s" % e)
        return float("nan")


def min_dist(a, b):
    try:
        return a.distToShape(b)[0]
    except Exception as e:
        log("distToShape failed: %s" % e)
        return float("nan")


# ----------------------------------------------------------------------------- inputs
def load_holder(doc):
    bodies = [o for o in doc.Objects if o.TypeId == "PartDesign::Body"]
    if not bodies:
        raise RuntimeError("no PartDesign Body (filter holder) in %s" % FCSTD)
    body = bodies[0]
    s = body.Shape.copy()
    s.Placement = App.Placement()
    bb = s.BoundBox
    info = dict(flat=bb.XMax, R=bb.YMax, t=bb.ZLength, hole_r=None)
    holes, aperture_r, pocket = [], None, None
    for f in s.Faces:
        srf = f.Surface
        if "Cylinder" not in srf.TypeId or abs(abs(srf.Axis.z) - 1) > 1e-6:
            continue
        r = srf.Radius
        c = (round(srf.Center.x, 3), round(srf.Center.y, 3))
        if r < 2.0:
            info["hole_r"] = r
            if c not in holes:
                holes.append(c)
        elif r < 4.0 and abs(c[0]) < 1e-3 and abs(c[1]) < 1e-3:
            aperture_r = r
        elif r > 10.0 and f.BoundBox.ZMax < bb.ZMax - 1e-3:
            pocket = (r, f.BoundBox.ZMin, f.BoundBox.ZMax)
    info.update(holes=holes, aperture_r=aperture_r, pocket=pocket)
    log("holder Body '%s': flats +-%.2f, R %.2f, t %.2f, %d holes r=%s, aperture r=%s, pocket %s"
        % (body.Label, info["flat"], info["R"], info["t"], len(holes), info["hole_r"],
           aperture_r, pocket))
    if abs(bb.ZMin) > 1e-6 or len(holes) != 4 or aperture_r is None or pocket is None:
        raise RuntimeError("holder Body does not look like the filter holder (z from 0, 4 holes, aperture, pocket)")
    return body, s, info


def load_ambit(P):
    if not os.path.isfile(AMBIT_STEP):
        raise RuntimeError("AMBIT STEP not found: %s (set AMBIT_STEP)" % AMBIT_STEP)
    sdoc = App.newDocument("ambit_step_tmp")
    Import.insert(AMBIT_STEP, sdoc.Name)
    casing, optics = [], []
    for o in sdoc.Objects:
        if o.TypeId != "Part::Feature" or o.Shape.isNull():
            continue
        s = o.Shape.copy()
        s.Placement = o.getGlobalPlacement()
        if any(k in o.Label for k in CASING_KEYS):
            casing.append(s)
        elif any(k in o.Label for k in OPTIC_KEYS):
            optics.append(s)
    App.closeDocument(sdoc.Name)
    if len(casing) < 4:
        raise RuntimeError("expected the AMBIT shells in the STEP, found %d casing parts" % len(casing))
    ymin = min(s.optimalBoundingBox().YMin for s in casing)
    # STEP leaf side (-Y) points down, head (+Z) toward jig -Y
    pl = App.Placement(V(0, P["ambit_centre_stepz"], P["ambit_gap"] - ymin),
                       App.Rotation(V(1, 0, 0), 90))
    mat = pl.toMatrix()
    casing = [transformed(s, mat) for s in casing]
    optics = [transformed(s, mat) for s in optics]
    bb = Part.makeCompound(casing).BoundBox
    log("AMBIT placed: %d casing parts, %d optics; casing x[%.2f,%.2f] y[%.2f,%.2f] z[%.2f,%.2f]"
        % (len(casing), len(optics), bb.XMin, bb.XMax, bb.YMin, bb.YMax, bb.ZMin, bb.ZMax))
    return casing, optics


# ----------------------------------------------------------------------------- AMBIT height map
def heightmap(casing, P):
    dx, dy = P["hm_dx"], P["hm_dy"]
    half = P["cradle_half"]
    xs = np.arange(-half, half + 1e-9, dx)
    y0 = P["cradle_y0"] - 1.0
    y1 = P["cradle_y1"] + 1.0
    ny = int(round((y1 - y0) / dy))
    low = np.full((len(xs), ny), np.inf)
    high = np.full((len(xs), ny), -np.inf)
    for i, x in enumerate(xs):
        for sh in casing:
            bb = sh.BoundBox
            if x < bb.XMin - 1e-6 or x > bb.XMax + 1e-6:
                continue
            for w in slice_x(sh, float(x)):
                try:
                    pts = w.discretize(Distance=dy * 0.4)
                except Exception:
                    continue
                if not pts:
                    continue
                a = np.array([(p.y, p.z) for p in pts])
                j = np.floor((a[:, 0] - y0) / dy).astype(int)
                ok = (j >= 0) & (j < ny)
                np.minimum.at(low[i], j[ok], a[ok, 1])
                np.maximum.at(high[i], j[ok], a[ok, 1])
    return dict(xs=xs, y0=y0, dy=dy, low=low, high=high)


HM_KEYS = ("hm_dx", "hm_dy", "cradle_half", "cradle_y0", "cradle_y1", "ambit_gap", "ambit_centre_stepz")


def cached_heightmap(casing, P):
    """Slicing the casing takes minutes; keep the result next to the script, keyed by inputs."""
    import hashlib
    st = os.stat(AMBIT_STEP)
    key = json.dumps([os.path.basename(AMBIT_STEP), st.st_size, int(st.st_mtime)] +
                     [P[k] for k in HM_KEYS])
    tag = hashlib.md5(key.encode()).hexdigest()[:12]
    cdir = os.path.join(HERE, ".jig_cache")
    path = os.path.join(cdir, "ambit_heightmap_%s.npz" % tag)
    if os.path.isfile(path):
        d = np.load(path)
        log("height map from cache %s" % path)
        return dict(xs=d["xs"], y0=float(d["y0"]), dy=float(d["dy"]), low=d["low"], high=d["high"])
    log("slicing the AMBIT casing for the height map (several minutes, cached afterwards)")
    hm = heightmap(casing, P)
    os.makedirs(cdir, exist_ok=True)
    np.savez_compressed(path, **hm)
    return hm


def build_profiles(hm, P, mode):
    """Dilated, merged column profiles: list of [x_left, x_right, profile[ny]]."""
    xs, dy = hm["xs"], hm["dy"]
    arr = hm["low"] if mode == "low" else hm["high"]
    nx, ny = arr.shape
    c = P["clear_xy"]
    mx = int(math.ceil(c / P["hm_dx"] - 1e-9))
    my = int(math.ceil(c / dy - 1e-9))
    red = np.min if mode == "low" else np.max
    fill = np.inf if mode == "low" else -np.inf
    strips = []
    for k in range(nx - 1):
        sub = arr[max(0, k - mx):min(nx, k + 2 + mx)]
        prof = red(sub, axis=0)
        ext = np.concatenate([np.full(my, fill), prof, np.full(my, fill)])
        prof = red(np.stack([ext[s:s + ny] for s in range(2 * my + 1)]), axis=0)
        strips.append([float(xs[k]), float(xs[k + 1]), prof])
    merged = []
    for s in strips:
        if not np.isfinite(s[2]).any():
            continue
        if merged:
            m = merged[-1]
            same_mask = np.array_equal(np.isfinite(m[2]), np.isfinite(s[2]))
            touching = abs(m[1] - s[0]) < 1e-9
            if same_mask and touching:
                fin = np.isfinite(s[2])
                if np.max(np.abs(m[2][fin] - s[2][fin])) < P["hm_merge"]:
                    m[1] = s[1]
                    m[2] = red(np.stack([m[2], s[2]]), axis=0)
                    continue
        merged.append([s[0], s[1], s[2].copy()])
    return merged


def profile_prisms(profiles, hm, P, mode, z_cap):
    """One prism per contiguous run: columns from the casing surface to z_cap."""
    y0, dy = hm["y0"], hm["dy"]
    tol = P["hm_simplify"]
    off = -(P["clear_floor"] + tol) if mode == "low" else (P["clear_lid"] + tol)
    pick = min if mode == "low" else max
    prisms = []
    for xl, xr, prof in profiles:
        fin = np.isfinite(prof)
        n = len(prof)
        j = 0
        while j < n:
            if not fin[j]:
                j += 1
                continue
            ja = j
            while j < n and fin[j]:
                j += 1
            jb = j - 1
            pts = []
            for b in range(ja, jb + 2):
                if b == ja:
                    v = prof[ja]
                elif b == jb + 1:
                    v = prof[jb]
                else:
                    v = pick(prof[b - 1], prof[b])
                pts.append((y0 + b * dy, v))
            pts = dp_simplify(pts, tol)
            # heights on a 0.01 mm grid, rounded toward more clearance, so neighbouring strips
            # either match or step by a printable amount (no sliver faces in the STL)
            q = (lambda z: math.floor(z * 100 + 1e-6) / 100) if mode == "low" else \
                (lambda z: math.ceil(z * 100 - 1e-6) / 100)
            poly = [V(xl, pts[0][0], z_cap)]
            poly += [V(xl, round(y, 4), q(z + off)) for y, z in pts]
            poly += [V(xl, pts[-1][0], z_cap), V(xl, pts[0][0], z_cap)]
            clean = [poly[0]]
            for p in poly[1:]:
                if (p - clean[-1]).Length > 1e-6:
                    clean.append(p)
            face = Part.Face(Part.makePolygon(clean))
            prisms.append(face.extrude(V(xr - xl, 0, 0)))
    return prisms


def footprint_half_x(profiles, hm, ya, yb):
    """Largest |x| of the dilated pocket between ya and yb."""
    y0, dy = hm["y0"], hm["dy"]
    ja = max(0, int(math.floor((ya - y0) / dy)))
    jb = int(math.ceil((yb - y0) / dy))
    ext = 0.0
    for xl, xr, prof in profiles:
        if np.isfinite(prof[ja:jb]).any():
            ext = max(ext, abs(xl), abs(xr))
    return ext


# ----------------------------------------------------------------------------- parts
def excitation_numbers(P):
    tan_t = math.tan(math.asin(P["fibre_NA"]))
    tip_from_face = P["ins_len"] + P["stud_len"] - P["ferrule_len"]
    L_min = P["E_stop"] + tip_from_face
    L_clip = P["E_stop"] + (P["E_stop_d"] / 2) / tan_t
    L_eng = P["E_tube_end"] - P["ins_min_engagement"] + tip_from_face
    L_max = min(L_clip, L_eng)
    face_max = L_max - tip_from_face
    if face_max >= P["E_stop"] + P["ins_len"]:
        raise RuntimeError("insert travel longer than its body: no common grub-screw position")
    s_set = 0.5 * (face_max + P["E_stop"] + P["ins_len"])
    return dict(tan_t=tan_t, tip_from_face=tip_from_face, L_min=L_min, L_max=L_max, s_set=s_set)


def build_insert(P):
    body = Part.makeCylinder(P["ins_d"] / 2, P["ins_len"])
    stud = Part.makeCylinder(P["stud_d"] / 2, P["stud_len"], V(0, 0, P["ins_len"]))
    top = P["ins_len"] + P["stud_len"]
    depth = P["ferrule_len"] + P["ferrule_bore_extra"]
    z_end = top - depth
    if z_end <= 0.5:
        raise RuntimeError("insert too short for the ferrule")
    bore = Part.makeCylinder(P["ferrule_bore_d"] / 2, depth + 1, V(0, 0, z_end))
    r_big = P["ferrule_bore_d"] / 2 + (z_end + 0.01) * math.tan(math.radians(P["ins_cone_half"]))
    cone = Part.makeCone(P["ferrule_bore_d"] / 2, r_big, z_end + 0.01, V(0, 0, z_end + 0.01), V(0, 0, -1))
    ins = single_solid(body.fuse(stud).cut([bore, cone]), "insert").removeSplitter()
    ins = single_solid(ins, "insert")
    log("insert: tip %.1f mm from the inner face, exit cone opens to d=%.2f (wall %.2f)"
        % (top - P["ferrule_len"], 2 * r_big, P["ins_d"] / 2 - r_big))
    return ins


def insert_matrix(face_s):
    return App.Placement(UE * face_s, App.Rotation(V(0, 0, 1), UE)).toMatrix()


def sma_envelope(P):
    sma = Part.makeCylinder(4.0, 9.0, V(0, 0, P["ins_len"] + 0.3))
    return sma.fuse([Part.makeCylinder(3.2, 17.0, V(0, 0, P["ins_len"] + 9.3)),
                     Part.makeCylinder(2.5, 20.0, V(0, 0, P["ins_len"] + 26.3))])


def build_cup(P, holder, hinfo, MD, cradle_prisms, heatset_xy, ex):
    H, h = P["box_half"], P["cav_half"]
    plate_t = hinfo["t"]
    D_inner = P["D_face"] - plate_t

    outer = Part.makeBox(2 * H, 2 * H, P["roof_top"], V(-H, -H, 0))
    dpx = (P["px_ceil_edge"] + P["cav_ceiling"]) / R2
    dpy = (P["py_ceil_edge"] + P["cav_ceiling"]) / R2
    outer = single_solid(outer.cut([halfspace_ge(UD, P["D_face"]),
                                    halfspace_ge(UE, P["E_slope_out"]),
                                    halfspace_ge(V(1, 0, 1), dpx + (H - h))]), "outer box")

    F = P["flange_half"]
    flange = Part.makeBox(2 * F, 2 * F, P["flange_t"], V(-F, -F, 0))
    # boss ring from the FiberPort face inward: the plate sits in it and the screw-head tunnels
    # (which reach past the plate outline) stay closed
    boss = transformed(holder_outline_prism(P["D_boss_R"], hinfo["flat"], 0.0,
                                            plate_t + D_inner - P["D_boss_start"]), MD)
    tube = Part.makeCylinder(P["E_tube_R"], P["E_tube_end"] - P["E_tube_start"],
                             UE * P["E_tube_start"], UE)
    ps = UE * ex["s_set"]
    hb = P["E_setscrew_boss_h"]
    diamond = Part.Face(Part.makePolygon([V(0, ps.y + hb, ps.z), V(0, ps.y, ps.z + hb),
                                          V(0, ps.y - hb, ps.z), V(0, ps.y, ps.z - hb),
                                          V(0, ps.y + hb, ps.z)]))
    ss_boss = diamond.extrude(V(P["E_setscrew_boss_end"], 0, 0))
    cw = P["cradle_half"]
    block = Part.makeBox(2 * cw, P["cradle_y1"] - P["cradle_y0"], P["parting_z"] - P["cradle_z0"],
                         V(-cw, P["cradle_y0"], P["cradle_z0"]))
    gy = P["gusset_y"]
    gz = P["cradle_z0"] - (P["cradle_y1"] - gy)
    gusset = Part.Face(Part.makePolygon([V(-cw, gy - 0.01, gz), V(-cw, gy - 0.01, P["cradle_z0"] + 0.01),
                                         V(-cw, P["cradle_y1"], P["cradle_z0"] + 0.01),
                                         V(-cw, gy - 0.01, gz)])).extrude(V(2 * cw, 0, 0))

    body = single_solid(outer.fuse([flange, boss, tube, ss_boss, block, gusset]), "cup outer")
    body = single_solid(body.removeSplitter(), "cup outer")
    log("cup outer body fused")

    plate_region = transformed(holder_outline_prism(hinfo["R"], hinfo["flat"], -1.0, plate_t), MD)
    plate = transformed(holder, MD)
    body = single_solid(body.cut(plate_region).fuse(plate), "cup with plate")
    log("holder plate fused at the detection port")

    cav = Part.makeBox(2 * h, 2 * h, P["cav_ceiling"] + 1, V(-h, -h, -1))
    cav = single_solid(cav.cut([halfspace_ge(UD, D_inner), halfspace_ge(UE, P["E_slope_in"]),
                                halfspace_ge(V(1, 0, 1), dpx), halfspace_ge(V(0, 1, 1), dpy)]), "cavity")
    tl = D_inner - P["D_tunnel_start"]
    tools = [cav]
    tools.append(transformed(Part.makeCylinder(P["D_light_tunnel_d"] / 2, tl, V(0, 0, plate_t)), MD))
    for hx, hy in hinfo["holes"]:
        tools.append(transformed(Part.makeCylinder(P["D_head_tunnel_d"] / 2, tl, V(hx, hy, plate_t)), MD))
    tools.append(Part.makeCylinder(P["E_bore_d"] / 2, P["E_tube_end"] - P["E_stop"] + 5,
                                   UE * P["E_stop"], UE))
    r0 = P["E_stop_d"] / 2
    tools.append(Part.makeCone(r0, r0 + 10 * math.tan(math.radians(P["E_cone_half"])), 10,
                               UE * P["E_stop"], UE * -1))
    x_in = P["E_bore_d"] / 2 - 0.5
    tools.append(Part.makeCylinder(P["heatset_d"] / 2, P["E_setscrew_boss_end"] + 1 - x_in,
                                   V(x_in, ps.y, ps.z), X))
    wx = P["window_half_x"]
    tools.append(Part.makeBox(2 * wx, P["window_y1"] - P["window_y0"], P["parting_z"] - 25.0,
                              V(-wx, P["window_y0"], 25.0)))
    for x, y in heatset_xy:
        tools.append(Part.makeCylinder(P["heatset_d"] / 2, P["heatset_depth"] + 1,
                                       V(x, y, P["parting_z"] - P["heatset_depth"])))
    tools += cradle_prisms
    log("cutting %d void tools (%d cradle strips)" % (len(tools), len(cradle_prisms)))
    cup = single_solid(body.cut(tools), "cup")
    log("cup cut done")
    try:
        cup = single_solid(cup.removeSplitter(), "cup refined")
    except Exception as e:
        log("removeSplitter skipped: %s" % e)
    return cup


def build_lid(P, lid_prisms, lid_top, heatset_xy):
    cw = P["cradle_half"]
    zp = P["parting_z"]
    lid = Part.makeBox(2 * cw, P["cradle_y1"] - P["cradle_y0"], lid_top - zp, V(-cw, P["cradle_y0"], zp))
    holes = [Part.makeCylinder(P["screw_clear_d"] / 2, lid_top - zp + 2, V(x, y, zp - 1))
             for x, y in heatset_xy]
    lid = single_solid(lid.cut(lid_prisms + holes), "lid")
    try:
        lid = single_solid(lid.removeSplitter(), "lid refined")
    except Exception as e:
        log("lid removeSplitter skipped: %s" % e)
    return lid


# ----------------------------------------------------------------------------- references
def excitation_cone(ex, L):
    length = 90.0
    cone = Part.makeCone(0.025, 0.025 + length * ex["tan_t"], length, UE * L, UE * -1)
    return cone.common(Part.makeBox(200, 200, 100, V(-100, -100, 0.0)))


def references(P, hinfo, MD, ex, insert):
    refs = {}
    fp = holder_outline_prism(hinfo["R"], hinfo["flat"], -P["fiberport_len"], 0.0)
    fc = Part.makeCylinder(6.0, 40.0, V(0, 0, -P["fiberport_len"] - 40.0))
    refs["ref_FiberPort_envelope"] = transformed(fp.fuse(fc), MD)
    refs["ref_filter"] = transformed(Part.makeCylinder(12.5, P["filter_t"]), MD)
    mi = insert_matrix(P["L_E_display"] - ex["tip_from_face"])
    refs["insert_at_display"] = transformed(insert, mi)
    refs["ref_SMA_connector_envelope"] = transformed(sma_envelope(P), mi)
    pocket_depth = hinfo["pocket"][2] - hinfo["pocket"][1]
    refs["ref_detection_view"] = Part.makeCylinder(2.0, P["D_face"] - pocket_depth - 0.2, V(0, 0, 0), UD)
    refs["ref_excitation_cone"] = excitation_cone(ex, P["L_E_display"])
    refs["ref_leaf"] = Part.makeBox(90, 110, 0.3, V(-45, -45, -0.3))
    return refs


# ----------------------------------------------------------------------------- document
def write_document(doc, body, MD, parts, refs, P, extra):
    old = doc.getObject("JigGenerated")
    if old is not None:
        for o in list(old.Group):
            doc.removeObject(o.Name)
        doc.removeObject(old.Name)
    grp = doc.addObject("App::DocumentObjectGroup", "JigGenerated")
    grp.Label = "jig (generated by jig_generator.py)"
    objs = {}
    for name, shape in list(parts.items()) + list(refs.items()):
        o = doc.addObject("Part::Feature", "Jig_" + name.replace("-", "_"))
        o.Label = name
        o.Shape = shape
        grp.addObject(o)
        objs[name] = o
    shown = ("ref_AMBIT", "ref_FiberPort_envelope", "insert_at_display",
             "ref_SMA_connector_envelope", "ref_excitation_cone")
    for name in refs:
        objs[name].Visibility = name in shown
    objs["jig-sma-insert"].Visibility = False   # the positioned copy is shown instead
    sheet = doc.addObject("Spreadsheet::Sheet", "JigParams")
    sheet.Label = "jig parameters (edit in jig_generator.py)"
    grp.addObject(sheet)
    sheet.set("A1", "parameter")
    sheet.set("B1", "value")
    row = 2
    for k, v in list(P.items()) + list(extra.items()):
        sheet.set("A%d" % row, k)
        sheet.set("B%d" % row, "%g" % v)
        try:
            sheet.setAlias("B%d" % row, k)
        except Exception:
            pass
        row += 1
    # the holder Body is placed where it sits in the cup and hidden (the cup contains it)
    body.Placement = App.Placement(MD)
    body.Visibility = False
    doc.recompute()
    return objs


COLOURS = {"jig-cup": (0.22, 0.22, 0.25), "jig-lid": (0.40, 0.40, 0.45), "jig-sma-insert": (0.85, 0.55, 0.10),
           "insert_at_display": (0.85, 0.55, 0.10), "ref_AMBIT": (0.20, 0.50, 0.85),
           "ref_FiberPort_envelope": (0.75, 0.75, 0.75), "ref_filter": (0.80, 0.20, 0.20),
           "ref_SMA_connector_envelope": (0.90, 0.80, 0.30), "ref_excitation_cone": (0.20, 0.80, 0.30),
           "ref_detection_view": (0.90, 0.20, 0.20), "ref_leaf": (0.30, 0.60, 0.20)}
TRANSPARENT = ("ref_FiberPort_envelope", "ref_excitation_cone", "ref_detection_view", "ref_leaf")


def style_in_gui(objs):
    """Colours only exist with a GUI (run as a macro from FreeCAD); freecadcmd skips this."""
    if not App.GuiUp:
        return
    for name, o in objs.items():
        vo = getattr(o, "ViewObject", None)
        if vo is None:
            continue
        if name in COLOURS:
            vo.ShapeColor = COLOURS[name]
        if name in TRANSPARENT:
            vo.Transparency = 60


def read_gui_entries(path):
    """freecadcmd saves without GuiDocument.xml; keep the previous view data to put back."""
    import zipfile
    import re
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if "GuiDocument.xml" not in names:
            return None
        xml = z.read("GuiDocument.xml").decode("utf-8")
        ent = [("GuiDocument.xml", xml)]
        if "thumbnails/Thumbnail.png" in names:
            ent.append(("thumbnails/Thumbnail.png", z.read("thumbnails/Thumbnail.png")))
        seen = set()
        for f in re.findall(r'file="([^"]+)"', xml):
            if f in names and f not in seen:
                seen.add(f)
                ent.append((f, z.read(f)))
    return ent


def restore_gui_entries(path, ent, hidden=("Body",)):
    import zipfile
    import re
    with zipfile.ZipFile(path) as z:
        existing = set(z.namelist())
    if ent is None or "GuiDocument.xml" in existing:
        return False
    xml = ent[0][1]
    for name in hidden:
        a = xml.find('<ViewProvider name="%s"' % name)
        if a < 0:
            continue
        b = xml.find("</ViewProvider>", a)
        block = re.sub(r'(<Property name="Visibility" type="App::PropertyBool"[^>]*>\s*<Bool value=")true(")',
                       r"\1false\2", xml[a:b], count=1)
        xml = xml[:a] + block + xml[b:]
    with zipfile.ZipFile(path, "a", zipfile.ZIP_DEFLATED) as z:
        z.writestr("GuiDocument.xml", xml)
        for name, data in ent[1:]:
            if name not in existing:
                z.writestr(name, data)
    return True


def export(objs, parts):
    """Write the STLs and a STEP of the printed parts; return mesh problems per part."""
    os.makedirs(EXPORT_DIR, exist_ok=True)
    problems = {}
    for name in parts:
        mesh = MeshPart.meshFromShape(Shape=parts[name], LinearDeflection=0.02,
                                      AngularDeflection=math.radians(8), Relative=False)
        path = os.path.join(EXPORT_DIR, name + ".stl")
        mesh.write(path)
        problems[name] = (mesh.isSolid(), mesh.hasNonManifolds(), len(mesh.getSelfIntersections()))
        log("wrote %s (%d facets)" % (path, mesh.CountFacets))
    step = os.path.join(EXPORT_DIR, "jig-parts.step")
    Import.export([objs[n] for n in parts], step)
    log("wrote %s" % step)
    return problems


def dump_preview(shapes, planes):
    os.makedirs(PREVIEW_DIR, exist_ok=True)
    for name, s in shapes.items():
        try:
            pts, tris = s.tessellate(0.15)
            np.savez_compressed(os.path.join(PREVIEW_DIR, "mesh_%s.npz" % name),
                                v=np.array([(p.x, p.y, p.z) for p in pts]), t=np.array(tris))
        except Exception as e:
            log("tessellate %s failed: %s" % (name, e))
    sections = {}
    for axis, val in planes:
        n = {"x": V(1, 0, 0), "y": V(0, 1, 0), "z": V(0, 0, 1)}[axis]
        key = "%s=%g" % (axis, val)
        sections[key] = {}
        for name, s in shapes.items():
            wires = []
            for d in (0.0, 0.011, -0.017):
                try:
                    wires = s.slice(n, val + d)
                    break
                except Exception:
                    wires = []
            polys = []
            for w in wires:
                try:
                    polys.append([[p.x, p.y, p.z] for p in w.discretize(Deflection=0.05)])
                except Exception:
                    pass
            sections[key][name] = polys
    with open(os.path.join(PREVIEW_DIR, "sections.json"), "w") as f:
        json.dump(sections, f)
    log("preview data written to %s" % PREVIEW_DIR)


# ----------------------------------------------------------------------------- checks
def light_check(P, cup, lid, hinfo, MD, ex, lid_top):
    """Fill every intended opening and check that the cavity air does not reach the outside."""
    cw = P["cradle_half"]
    fills = [lid,
             Part.makeBox(2 * cw, P["cradle_y1"] - P["cradle_y0"], lid_top - P["cradle_z0"],
                          V(-cw, P["cradle_y0"], P["cradle_z0"])),           # AMBIT in its nest
             transformed(holder_outline_prism(hinfo["R"] + 0.5, hinfo["flat"] + 0.5,
                                              -10.0, 0.5), MD),               # filter + FiberPort
             Part.makeCylinder(P["E_bore_d"] / 2 + 0.2, 60, UE * P["E_stop"], UE),  # insert + SMA
             Part.makeCylinder(P["heatset_d"] / 2 + 0.2, 20, UE * ex["s_set"], X),  # grub screw
             Part.makeBox(160, 180, 0.8, V(-80, -80, -0.5))]                   # leaf + foam
    for hx, hy in hinfo["holes"]:                                              # FiberPort screws
        fills.append(transformed(Part.makeCylinder(hinfo["hole_r"] + 0.2, hinfo["t"] + 1, V(hx, hy, -0.5)), MD))
    dom = Part.makeBox(150, 170, 90, V(-75, -75, -0.4))
    air = dom.cut([cup] + fills)
    probe = V(0, 0, 8)
    for s in air.Solids:
        if s.isInside(probe, 1e-6, True):
            bb = s.BoundBox
            inside = (bb.XMin > -40 and bb.XMax < 40 and bb.YMin > -40 and bb.YMax < 40 and bb.ZMax < 45)
            return inside, bb, len(air.Solids)
    return False, None, len(air.Solids)


def overhangs(shape, limit_nz=-0.72, skip_z=(0.0,)):
    bad = []
    for f in shape.Faces:
        try:
            u0, u1, v0, v1 = f.ParameterRange
            n = f.normalAt(0.5 * (u0 + u1), 0.5 * (v0 + v1))
        except Exception:
            continue
        if n.z < limit_nz and f.Area > 2.0:
            c = f.CenterOfMass
            if "Plane" in f.Surface.TypeId and any(abs(c.z - z) < 1e-3 for z in skip_z):
                continue
            bad.append((round(c.x, 1), round(c.y, 1), round(c.z, 1), round(f.Area, 1), round(n.z, 2)))
    return bad


# ----------------------------------------------------------------------------- main
def main():
    P = PARAMS
    log("document %s" % FCSTD)
    gui_entries = None if App.GuiUp else read_gui_entries(FCSTD)
    doc = App.openDocument(FCSTD)
    body, holder, hinfo = load_holder(doc)
    # holder local frame: x -> -VD (across the flats), y -> Y, z -> -UD (inward); z=0 face at D_face
    MD = matrix_from_axes(VD * -1, V(0, 1, 0), UD * -1, UD * P["D_face"])
    ex = excitation_numbers(P)
    log("excitation: fibre tip %.1f..%.1f mm from the spot, grub screw at %.1f mm"
        % (ex["L_min"], ex["L_max"], ex["s_set"]))

    casing, optics = load_ambit(P)
    hm = cached_heightmap(casing, P)
    log("height map: %d slices x %d bins" % hm["low"].shape)

    xs, y0, dy = hm["xs"], hm["y0"], hm["dy"]
    ix = np.where(np.abs(xs) <= P["window_half_x"] + 1e-9)[0]
    jy = np.arange(int(round((P["window_y0"] - y0) / dy)), int(round((P["window_y1"] - y0) / dy)))
    window_covered = bool(np.isfinite(hm["low"][np.ix_(ix, jy)]).all())

    low_prof = build_profiles(hm, P, "low")
    high_prof = build_profiles(hm, P, "high")
    zp = P["parting_z"]
    cradle_prisms = profile_prisms(low_prof, hm, P, "low", zp + 5.0)
    lid_prisms = profile_prisms(high_prof, hm, P, "high", zp - 5.0)
    yl = np.arange(hm["high"].shape[1]) * dy + y0
    inlid = (yl >= P["cradle_y0"]) & (yl <= P["cradle_y1"])
    sub = hm["high"][:, inlid]
    zmax = float(np.max(sub[np.isfinite(sub)]))
    lid_top = round(zmax + P["clear_lid"] + P["hm_simplify"] + P["lid_wall"], 2)
    log("cradle: %d strips, lid: %d strips, casing top under the lid %.2f -> lid top %.2f"
        % (len(cradle_prisms), len(lid_prisms), zmax, lid_top))

    heatset_xy = []
    r_ins = P["heatset_d"] / 2
    for ys in (P["screw_y_head"], P["screw_y_tail"]):
        half = max(footprint_half_x(low_prof, hm, ys - r_ins - 1, ys + r_ins + 1),
                   footprint_half_x(high_prof, hm, ys - r_ins - 1, ys + r_ins + 1))
        xsc = round(half + r_ins + P["screw_wall_min"], 2)
        if xsc + r_ins + 0.8 > P["cradle_half"]:
            raise RuntimeError("no room for the lid screw at y=%.1f (pocket half width %.2f)" % (ys, half))
        heatset_xy += [(-xsc, ys), (xsc, ys)]
    log("lid screws at %s" % heatset_xy)

    insert = build_insert(P)
    cup = build_cup(P, holder, hinfo, MD, cradle_prisms, heatset_xy, ex)
    log("cup: volume %.1f cm3, valid %s, faces %d" % (cup.Volume / 1000, cup.isValid(), len(cup.Faces)))
    lid = build_lid(P, lid_prisms, lid_top, heatset_xy)
    log("lid: volume %.1f cm3, valid %s" % (lid.Volume / 1000, lid.isValid()))

    parts = {"jig-cup": cup, "jig-lid": lid, "jig-sma-insert": insert}
    refs = references(P, hinfo, MD, ex, insert)
    refs["ref_AMBIT"] = Part.makeCompound(casing + optics)

    # ------------------------------------------------------------------ checks
    results = []

    def check(name, ok, detail):
        results.append((name, ok, detail))
        log("%s %-58s %s" % ("PASS" if ok else "FAIL", name, detail))

    check("cup, lid, insert are valid single solids",
          all(s.isValid() and len(s.Solids) == 1 for s in parts.values()), "")
    check("window under the AMBIT covered by the casing", window_covered,
          "x +-%.1f, y %.1f..%.1f" % (P["window_half_x"], P["window_y0"], P["window_y1"]))
    case_c = Part.makeCompound(casing)
    for pname in ("jig-cup", "jig-lid"):
        d = min_dist(case_c, parts[pname])
        v = 0.0 if d > 1e-4 else common_volume(case_c, parts[pname])
        check("AMBIT casing clear of %s" % pname, d > 1e-4 or v < 0.01,
              "min gap %.3f mm, overlap %.4f mm3" % (d, v))
    d = min_dist(refs["ref_detection_view"], cup)
    check("detection view (d 4 mm) clear of the cup", d > 0.05, "min gap %.2f mm" % d)
    for L in (P["L_E_display"], ex["L_max"]):
        cone = excitation_cone(ex, L)
        v = common_volume(cone, cup)
        check("excitation cone at L=%.1f clear of the cup" % L, v < 0.5, "overlap %.3f mm3" % v)
        v = common_volume(cone, case_c)
        check("excitation cone at L=%.1f clear of the AMBIT" % L, v < 0.01, "overlap %.3f mm3" % v)
    fp = refs["ref_FiberPort_envelope"]
    for other, s in (("cup", cup), ("lid", lid), ("AMBIT", case_c)):
        v = common_volume(fp, s)
        check("FiberPort envelope clear of the %s" % other, v < 0.01, "overlap %.4f mm3" % v)
    for L in (ex["L_min"], ex["L_max"]):
        mi = insert_matrix(L - ex["tip_from_face"])
        v = common_volume(transformed(insert, mi), cup)
        check("insert at L=%.1f slides free in the tube" % L, v < 0.01, "overlap %.4f mm3" % v)
        v = common_volume(transformed(sma_envelope(P), mi), cup)
        check("SMA connector at L=%.1f clear of the cup" % L, v < 0.01, "overlap %.4f mm3" % v)
    cradle_c = Part.makeCompound(cradle_prisms)
    for x, y in heatset_xy:
        hole = Part.makeCylinder(P["heatset_d"] / 2, P["heatset_depth"], V(x, y, zp - P["heatset_depth"]))
        d = min_dist(hole, cradle_c)
        check("heat-set hole (%.1f, %.1f) wall to the AMBIT pocket" % (x, y), d >= 0.8, "%.2f mm" % d)
    if not SKIP_LIGHTCHECK:
        try:
            ok, bb, n = light_check(P, cup, lid, hinfo, MD, ex, lid_top)
            det = ("cavity air x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f], %d air volumes"
                   % (bb.XMin, bb.XMax, bb.YMin, bb.YMax, bb.ZMin, bb.ZMax, n)) if bb else "probe not found"
            check("cavity sealed with FiberPort, insert, AMBIT and leaf in place", ok, det)
        except Exception as e:
            check("cavity sealed (light check)", False, "check failed to run: %s" % e)
    oh = overhangs(cup, skip_z=(0.0, P["cav_ceiling"]))
    log("cup faces overhanging more than 44 deg when printed rim-down (x, y, z, area, nz): %s" % oh[:15])

    extra = dict(lid_top=lid_top, L_E_min=round(ex["L_min"], 2), L_E_max=round(ex["L_max"], 2),
                 setscrew_s=round(ex["s_set"], 2), lid_screw_x=heatset_xy[1][0])
    objs = write_document(doc, body, MD, parts, refs, P, extra)
    style_in_gui(objs)
    doc.save()
    log("saved %s" % FCSTD)
    if not App.GuiUp and restore_gui_entries(FCSTD, gui_entries):
        log("kept the previous view settings (GuiDocument.xml) in the saved file")
    for name, (closed, nonman, nsi) in export(objs, parts).items():
        check("%s.stl closed, manifold, no self-intersections" % name, closed and not nonman and nsi == 0,
              "closed %s, non-manifold %s, self-intersecting pairs %d" % (closed, nonman, nsi))
    if PREVIEW_DIR:
        shapes = {"jig-cup": cup, "jig-lid": lid, "jig-sma-insert": refs["insert_at_display"]}
        for k in ("ref_AMBIT", "ref_FiberPort_envelope", "ref_SMA_connector_envelope",
                  "ref_excitation_cone", "ref_detection_view"):
            shapes[k] = refs[k]
        dump_preview(shapes, [("x", 0.0), ("y", 0.0), ("y", 3.0), ("z", 33.0), ("z", 28.5), ("z", 12.0)])
    fails = [r for r in results if not r[1]]
    log("%d checks, %d failed" % (len(results), len(fails)))
    for r in fails:
        log("FAILED: %s %s" % (r[0], r[2]))
    if not App.GuiUp:
        App.closeDocument(doc.Name)
    return 0 if not fails else 1


main()
