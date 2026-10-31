#!/usr/bin/env python3
"""
Case generator for the Turek-Hron benchmark (FSI2, FSI3): pimpleFoam with the
rigidBodyMotionSolver moving-mesh solver. The deformable flag is modelled as a chain
of N rigid segments connected by elastic linearJointSpring hinges.

Usage:
    python3 scripts/make_case.py --case FSI2 --n-segments 20 --cells-per-segment 6
then ./Allrun
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent.parent

# Geometry (Turek, Hron, 2006), m
LX, LY = 2.5, 0.41
XC, YC, R = 0.2, 0.2, 0.05
H_FLAG = 0.02
X_TIP = 0.6
Y_F1, Y_F2 = YC - H_FLAG/2, YC + H_FLAG/2
X_ROOT = XC + math.sqrt(R**2 - (H_FLAG/2)**2)
L_FLAG = X_TIP - X_ROOT
XA_BOX, YB_BOX, YT_BOX = 0.063, 0.063, 0.337
Z0, Z1 = -0.5, 0.5

# Case parameters (Turek, Hron, 2006), SI units
CASES = {
    "FSI2": dict(rho_s=10000.0, E=1.4e6, Um=1.0, rho_f=1000.0, g=0.0, T=15.0, dt=1e-3),
    "FSI3": dict(rho_s=1000.0,  E=5.6e6, Um=2.0, rho_f=1000.0, g=0.0, T=10.0, dt=5e-4),
}
NU_S = 0.4
NU_F = 1e-3

REFERENCE = {
    "FSI2": "ux = -14.58 +- 12.44, uy = 1.23 +- 80.6 (x1e-3 m);  drag = 208.83 +- 73.75, lift = 0.88 +- 234.2 N/m; "
            "frequency of uy ~1.93-2.0 Hz, of ux and drag ~3.86 Hz",
    "FSI3": "ux = -2.69 +- 2.53, uy = 1.48 +- 34.38 (x1e-3 m);  drag = 457.3 +- 22.66, lift = 2.22 +- 149.78 N/m; "
            "frequency of uy ~5.3-5.5 Hz, of ux and drag ~10.9 Hz",
}

FOAM_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2112                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    {loc}object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
FOOTER = "\n// ************************************************************************* //\n"


def header(obj, cls="dictionary", loc=None):
    return FOAM_HEADER.format(obj=obj, cls=cls, loc=f'location    "{loc}";\n    ' if loc else "")


def geometric_ratio(length, n, first):
    if first*n >= length:
        return 1.0
    lo, hi = 1.0 + 1e-12, 2.0
    f = lambda q: first*(q**n - 1)/(q - 1) - length
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    q = 0.5*(lo + hi)
    return q**(n - 1)

class BlockMesh:
    def __init__(self):
        self.pts2d = {}
        self.order = []
        self.blocks = []
        self.arcs = []
        self.patches = {}

    def p(self, name, x, y):
        self.pts2d[name] = (float(x), float(y))
        self.order.append(name)

    def block(self, name, v, n, g=(1.0, 1.0)):
        self.blocks.append(dict(name=name, v=v, n=n, g=g))

    def arc(self, a, b, mid):
        self.arcs.append((a, b, mid))

    def patch(self, name, ptype, faces):
        self.patches.setdefault(name, (ptype, []))[1].extend(faces)

    def vid(self, name, back=False):
        return 2*self.order.index(name) + (1 if back else 0)

    SIDES = {"x-": (0, 3), "x+": (1, 2), "y-": (0, 1), "y+": (3, 2)}

    def face_names(self, blk, side):
        a, b = self.SIDES[side]
        return blk["v"][a], blk["v"][b]

    def face_ids(self, blk, side):
        a, b = self.face_names(blk, side)
        return (self.vid(a), self.vid(b), self.vid(b, True), self.vid(a, True))

    @staticmethod
    def nodes(n, g):
        if abs(g - 1) < 1e-12:
            return np.linspace(0, 1, n + 1)
        q = g**(1.0/(n - 1))
        w = q**np.arange(n)
        return np.concatenate([[0], np.cumsum(w)/w.sum()])

    def edge_dist(self, blk, a_idx, b_idx):
        pair = {(0, 1): 0, (3, 2): 0, (0, 3): 1, (1, 2): 1}
        rev = {(1, 0): 0, (2, 3): 0, (3, 0): 1, (2, 1): 1}
        if (a_idx, b_idx) in pair:
            d = pair[(a_idx, b_idx)]
            return self.nodes(blk["n"][d], blk["g"][d])
        d = rev[(a_idx, b_idx)]
        t = self.nodes(blk["n"][d], blk["g"][d])
        return 1 - t[::-1]

    # geometry and topology checks
    def validate(self):
        errs = []
        P = self.pts2d
        for b in self.blocks:
            q = [P[n] for n in b["v"]]
            area = 0.5*sum(q[i][0]*q[(i+1) % 4][1] - q[(i+1) % 4][0]*q[i][1] for i in range(4))
            if area <= 0:
                errs.append(f"block {b['name']}: clockwise vertex order (area {area:.3e})")
            for i in range(4):
                a, c, d = q[i], q[(i+1) % 4], q[(i+2) % 4]
                cross = (c[0]-a[0])*(d[1]-c[1]) - (c[1]-a[1])*(d[0]-c[0])
                if cross <= 0:
                    errs.append(f"block {b['name']}: non-convex corner at vertex {b['v'][(i+1) % 4]}")
        edges = {}
        for b in self.blocks:
            for a, c in ((0, 1), (1, 2), (2, 3), (3, 0)):
                key = frozenset((b["v"][a], b["v"][c]))
                na, nc = b["v"][a], b["v"][c]
                if na < nc:
                    dist = self.edge_dist(b, a, c)
                else:
                    dist = 1 - self.edge_dist(b, a, c)[::-1]
                edges.setdefault(key, []).append((b["name"], dist))
        for key, lst in edges.items():
            for name, dist in lst[1:]:
                if len(dist) != len(lst[0][1]) or np.max(np.abs(dist - lst[0][1])) > 1e-9:
                    errs.append(f"edge {sorted(key)}: node mismatch between {lst[0][0]} and {name}")
        faces = {}
        for b in self.blocks:
            for side in self.SIDES:
                key = frozenset(self.face_names(b, side))
                faces.setdefault(key, []).append((b["name"], side))
        in_patch = {}
        for pname, (_, flist) in self.patches.items():
            for bname, side in flist:
                b = next(x for x in self.blocks if x["name"] == bname)
                key = frozenset(self.face_names(b, side))
                if key in in_patch:
                    errs.append(f"face {sorted(key)} belongs to two patches: {in_patch[key]} and {pname}")
                in_patch[key] = pname
        for key, lst in faces.items():
            if len(lst) > 2:
                errs.append(f"face {sorted(key)} is shared by {len(lst)} blocks")
            elif len(lst) == 2 and key in in_patch:
                errs.append(f"internal face {sorted(key)} is assigned to patch {in_patch[key]}")
            elif len(lst) == 1 and key not in in_patch:
                errs.append(f"boundary face {sorted(key)} of block {lst[0][0]} is not assigned to any patch")
        used = {n for b in self.blocks for n in b["v"]}
        for n in self.order:
            if n not in used:
                errs.append(f"vertex {n} is not used by any block")
        for a, c, m in self.arcs:
            for pt in (P[a], P[c], m):
                if abs(math.hypot(pt[0]-XC, pt[1]-YC) - R) > 1e-9:
                    errs.append(f"arc {a}-{c}: point {pt} is not on the circle")
        return errs

    def n_cells(self):
        return sum(b["n"][0]*b["n"][1] for b in self.blocks)

    def write(self, path):
        s = header("blockMeshDict") + "scale   1;\n\nvertices\n(\n"
        for name in self.order:
            x, y = self.pts2d[name]
            s += f"    ({x:.10f} {y:.10f} {Z0})   // {2*self.order.index(name)} {name}\n"
            s += f"    ({x:.10f} {y:.10f} {Z1})   // {2*self.order.index(name)+1}\n"
        s += ");\n\nblocks\n(\n"
        for b in self.blocks:
            v = b["v"]
            ids = [self.vid(n) for n in v] + [self.vid(n, True) for n in v]
            s += (f"    hex ({' '.join(map(str, ids))}) {b['name']} ({b['n'][0]} {b['n'][1]} 1) "
                  f"simpleGrading ({b['g'][0]:.8g} {b['g'][1]:.8g} 1)\n")
        s += ");\n\nedges\n(\n"
        for a, c, (xm, ym) in self.arcs:
            s += f"    arc {self.vid(a)} {self.vid(c)} ({xm:.10f} {ym:.10f} {Z0})\n"
            s += f"    arc {self.vid(a, True)} {self.vid(c, True)} ({xm:.10f} {ym:.10f} {Z1})\n"
        s += ");\n\nboundary\n(\n"
        for pname, (ptype, flist) in self.patches.items():
            s += f"    {pname}\n    {{\n        type {ptype};\n        faces\n        (\n"
            for bname, side in flist:
                b = next(x for x in self.blocks if x["name"] == bname)
                s += f"            ({' '.join(map(str, self.face_ids(b, side)))})   // {bname} {side}\n"
            s += "        );\n    }\n"
        s += "    frontAndBack\n    {\n        type empty;\n        faces\n        (\n"
        for b in self.blocks:
            v = b["v"]
            s += f"            ({' '.join(str(self.vid(n)) for n in v[::-1])})\n"
            s += f"            ({' '.join(str(self.vid(n, True)) for n in v)})\n"
        s += "        );\n    }\n);\n\nmergePatchPairs\n(\n);\n" + FOOTER
        Path(path).write_text(s)


def build_mesh(n_seg, cps, m_box):
    dx = L_FLAG/n_seg
    hf = dx/cps
    xb = X_ROOT + m_box*dx
    if not (XC + 0.75*R < xb - 0.04 or xb > XC + 1.6*R):
        raise ValueError("right boundary of the block around the cylinder is too close to the cylinder; increase m_box")

    th_f = math.asin((H_FLAG/2)/R)
    corner = lambda x, y: math.atan2(y - YC, x - XC)
    th_TR, th_TL = corner(xb, YT_BOX), corner(XA_BOX, YT_BOX)
    th_BL, th_BR = corner(XA_BOX, YB_BOX) + 2*math.pi, corner(xb, YB_BOX) + 2*math.pi
    cpt = lambda th: (XC + R*math.cos(th), YC + R*math.sin(th))
    mid = lambda a, b: cpt(0.5*(a + b))

    bm = BlockMesh()
    P = bm.p
    P("cRT", *cpt(th_f)); P("cTR", *cpt(th_TR)); P("cTL", *cpt(th_TL))
    P("cBL", *cpt(th_BL)); P("cBR", *cpt(th_BR)); P("cRB", *cpt(-th_f))
    xs = {"0": 0.0, "a": XA_BOX, "b": xb, "t": X_TIP, "e": LX}
    ys = {"0": 0.0, "b": YB_BOX, "1": Y_F1, "2": Y_F2, "t": YT_BOX, "e": LY}
    grid = [("a", "b"), ("a", "t"), ("b", "b"), ("b", "1"), ("b", "2"), ("b", "t"),
            ("0", "0"), ("0", "b"), ("0", "t"), ("0", "e"), ("a", "0"), ("a", "e"),
            ("b", "0"), ("b", "e"), ("t", "0"), ("t", "b"), ("t", "1"), ("t", "2"), ("t", "t"), ("t", "e"),
            ("e", "0"), ("e", "b"), ("e", "1"), ("e", "2"), ("e", "t"), ("e", "e")]
    for gx, gy in grid:
        P(f"g{gx}{gy}", xs[gx], ys[gy])

    n_r = m_box*cps
    n_fx = (n_seg - m_box)*cps
    n_q = max(16, round((0.5*math.pi*R)/(1.3*hf)))
    n_ur = max(12, round(0.6*(YT_BOX - Y_F2)/hf))
    n_tip = 2*max(3, round(H_FLAG/(2*hf)))
    far = 2.2*hf
    n_top = max(10, round((LY - YT_BOX)/far))
    n_bot = max(10, round(YB_BOX/far))
    n_left = max(10, round(XA_BOX/far))
    n_w = max(60, round(1.25*(LX - X_TIP)/(8*hf)))
    g_w = geometric_ratio(LX - X_TIP, n_w, hf)
    g_f = geometric_ratio(YT_BOX - Y_F2, n_ur, 1.2*hf)

    bm.block("Bur", ["cRT", "gb2", "gbt", "cTR"], (n_r, n_ur), (1, g_f))
    bm.block("Btop", ["cTR", "gbt", "gat", "cTL"], (n_r, n_q))
    bm.block("Bleft", ["cTL", "gat", "gab", "cBL"], (n_r, n_q))
    bm.block("Bbot", ["cBL", "gab", "gbb", "cBR"], (n_r, n_q))
    bm.block("Blr", ["cBR", "gbb", "gb1", "cRB"], (n_r, n_ur), (1, 1/g_f))
    bm.arc("cRT", "cTR", mid(th_f, th_TR))
    bm.arc("cTR", "cTL", mid(th_TR, th_TL))
    bm.arc("cTL", "cBL", mid(th_TL, th_BL))
    bm.arc("cBL", "cBR", mid(th_BL, th_BR))
    bm.arc("cBR", "cRB", mid(th_BR, 2*math.pi - th_f))
    
    bm.block("Ftop", ["gb2", "gt2", "gtt", "gbt"], (n_fx, n_ur), (1, g_f))
    bm.block("Fbot", ["gbb", "gtb", "gt1", "gb1"], (n_fx, n_ur), (1, 1/g_f))
    
    bm.block("Lbot", ["g00", "ga0", "gab", "g0b"], (n_left, n_bot))
    bm.block("Lmid", ["g0b", "gab", "gat", "g0t"], (n_left, n_q))
    bm.block("Ltop", ["g0t", "gat", "gae", "g0e"], (n_left, n_top))
    bm.block("Obot", ["ga0", "gb0", "gbb", "gab"], (n_q, n_bot))
    bm.block("Otop", ["gat", "gbt", "gbe", "gae"], (n_q, n_top))
    bm.block("FbotB", ["gb0", "gt0", "gtb", "gbb"], (n_fx, n_bot))
    bm.block("FtopB", ["gbt", "gtt", "gte", "gbe"], (n_fx, n_top))
    
    bm.block("W0", ["gt0", "ge0", "geb", "gtb"], (n_w, n_bot), (g_w, 1))
    bm.block("W1", ["gtb", "geb", "ge1", "gt1"], (n_w, n_ur), (g_w, 1/g_f))
    bm.block("W2", ["gt1", "ge1", "ge2", "gt2"], (n_w, n_tip), (g_w, 1))
    bm.block("W3", ["gt2", "ge2", "get", "gtt"], (n_w, n_ur), (g_w, g_f))
    bm.block("W4", ["gtt", "get", "gee", "gte"], (n_w, n_top), (g_w, 1))

    bm.patch("inlet", "patch", [("Lbot", "x-"), ("Lmid", "x-"), ("Ltop", "x-")])
    bm.patch("outlet", "patch", [(f"W{i}", "x+") for i in range(5)])
    bm.patch("walls", "wall", [("Lbot", "y-"), ("Obot", "y-"), ("FbotB", "y-"), ("W0", "y-"),
                                ("Ltop", "y+"), ("Otop", "y+"), ("FtopB", "y+"), ("W4", "y+")])
    bm.patch("cylinder", "wall", [(b, "x-") for b in ("Bur", "Btop", "Bleft", "Bbot", "Blr")])
    bm.patch("flag", "wall", [("Bur", "y-"), ("Ftop", "y-"), ("Blr", "y+"), ("Fbot", "y+"), ("W2", "x-")])

    errs = bm.validate()
    info = dict(dx=dx, hf=hf, xb=xb, n_r=n_r, n_fx=n_fx, n_q=n_q, n_ur=n_ur, n_tip=n_tip,
                n_top=n_top, n_bot=n_bot, n_left=n_left, n_w=n_w, g_w=g_w, g_f=g_f,
                cells=bm.n_cells())
    return bm, errs, info


def write_topo_and_patch(case, n_seg, dx):
    t = header("topoSetDict") + "actions\n(\n"
    cp = header("createPatchDict") + "pointSync false;\n\npatches\n(\n"
    for i in range(n_seg):
        x1 = X_ROOT + i*dx - (0.01 if i == 0 else 0.0)
        x2 = X_ROOT + (i + 1)*dx + (0.01 if i == n_seg - 1 else 0.0)
        t += f"""    {{ name seg{i}; type faceSet; action new;    source patchToFace; patch flag; }}
    {{ name seg{i}; type faceSet; action subset; source boxToFace;
      box ({x1:.10f} {Y_F1 - 0.005:.4f} {Z0 - 1}) ({x2:.10f} {Y_F2 + 0.005:.4f} {Z1 + 1}); }}
"""
        cp += f"""    {{
        name segment{i}Wall;
        patchInfo {{ type wall; }}
        constructFrom set;
        set seg{i};
    }}
"""
    t += ");\n" + FOOTER
    cp += ");\n" + FOOTER
    (case/"system"/"topoSetDict").write_text(t)
    (case/"system"/"createPatchDict").write_text(cp)


def write_dynamic_mesh_dict(case, n_seg, prm, eta, n_outer_relax):
    dx = L_FLAG/n_seg
    EI = prm["E"]/(1 - NU_S**2)*H_FLAG**3/12
    m = prm["rho_s"]*H_FLAG*dx*(Z1 - Z0)
    W = Z1 - Z0
    Ixx, Iyy, Izz = m*(H_FLAG**2 + W**2)/12, m*(dx**2 + W**2)/12, m*(dx**2 + H_FLAG**2)/12
    patches = " ".join(f"segment{i}Wall" for i in range(n_seg))

    s = header("dynamicMeshDict")
    s += f"""
dynamicFvMesh   dynamicMotionSolverFvMesh;

// librigidBodyDynamicsSpring.so: library with the linearJointSpring restraint (src/linearJointSpring)
motionSolverLibs ("librigidBodyMeshMotion.so" "libfvMotionSolvers.so" "librigidBodyDynamicsSpring.so");

motionSolver    rigidBodyMotionSolver;

meshSolver
{{
    solver          displacementLaplacian;
    diffusivity     inverseDistance (cylinder {patches});
}}

report          on;

solver
{{
    type            Newmark;
    gamma           0.5;
    beta            0.25;
}}

accelerationRelaxation  {n_outer_relax};
accelerationDamping     1;
nIter                   1;

g               (0 {(-prm['g'] if prm['g'] else 0):.6g} 0);
rho             rhoInf;
rhoInf          {prm['rho_f']:.6g};

bodies
{{
"""
    for i in range(n_seg):
        parent = "root" if i == 0 else f"segment{i-1}"
        org = f"({X_ROOT:.10f} {YC} 0)" if i == 0 else f"({dx:.10f} 0 0)"
        s += f"""    segment{i}
    {{
        type            rigidBody;
        parent          {parent};
        centreOfMass    ({dx/2:.10f} 0 0);
        mass            {m:.10g};
        inertia         ({Ixx:.10g} 0 0 {Iyy:.10g} 0 {Izz:.10g});
        transform       (1 0 0 0 1 0 0 0 1) {org};
        joint           {{ type Rz; }}
        patches         (segment{i}Wall);
    }}
"""
    s += "}\n\nrestraints\n{\n"
    for i in range(n_seg):
        k = 2*EI/dx if i == 0 else EI/dx      # k = 2EI/(dx_{i-1} + dx_i), dx_{-1} = 0
        s += f"""    hinge{i}{'   // clamp hinge, k = 2EI/dx' if i == 0 else ''}
    {{
        type            linearJointSpring;
        body            segment{i};
        stiffness       {k:.10g};
        damping         {eta*k:.10g};
    }}
"""
    s += "}\n" + FOOTER
    (case/"constant"/"dynamicMeshDict").write_text(s)
    return EI, m

def write_U(case, Um):
    s = header("U", "volVectorField", "0") + f"""dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{{
    inlet
    {{
        type            codedFixedValue;
        value           uniform (0 0 0);
        name            turekInlet;
        codeInclude
        #{{
            #include "mathematicalConstants.H"
        #}};
        code
        #{{
            const scalar Um = {Um:.10g};
            const scalar H  = {LY};
            const scalar t  = this->db().time().value();
            const scalar ramp =
                (t < 2.0) ? 0.5*(1.0 - cos(constant::mathematical::pi*t/2.0)) : 1.0;

            const vectorField& Cf = patch().Cf();
            vectorField Uin(Cf.size(), Zero);
            forAll(Cf, i)
            {{
                const scalar y = Cf[i].y();
                Uin[i] = vector(1.5*Um*y*(H - y)/sqr(0.5*H)*ramp, 0, 0);
            }}
            operator==(Uin);
        #}};
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);
    }}

    "(walls|cylinder)"
    {{
        type            noSlip;
    }}

    "(flag|segment.*Wall)"
    {{
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }}

    frontAndBack
    {{
        type            empty;
    }}
}}
""" + FOOTER
    (case/"0.orig"/"U").write_text(s)


def write_transport(case):
    s = header("transportProperties") + f"""transportModel  Newtonian;

nu              {NU_F};
""" + FOOTER
    (case/"constant"/"transportProperties").write_text(s)
    

def write_control(case, prm, n_seg):
    s = header("controlDict") + f"""application     pimpleFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {prm['T']};
deltaT          {prm['dt']};
adjustTimeStep  no;

writeControl    adjustable;
writeInterval   0.1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;

functions
{{
    forcesTurek
    {{
        type            forces;
        libs            ("libforces.so");
        writeControl    timeStep;
        writeInterval   1;
        patches         (cylinder "segment.*Wall");
        rho             rhoInf;
        rhoInf          {prm['rho_f']:.6g};
        CofR            ({XC} {YC} 0);
        log             false;
    }}

    tipA
    {{
        type            coded;
        libs            ("libutilityFunctionObjects.so");
        name            tipA;

        codeData
        #{{
            label tipPoint_ = -2;
            vector A0_ = vector({X_TIP}, {YC}, 0);
        #}};

        codeExecute
        #{{
            const pointField& pts = mesh().points();
            if (tipPoint_ == -2)
            {{
                scalar dmin = GREAT;
                label imin = -1;
                forAll(pts, i)
                {{
                    const scalar d = Foam::sqrt(sqr(pts[i].x() - A0_.x()) + sqr(pts[i].y() - A0_.y()));
                    if (d < dmin) {{ dmin = d; imin = i; }}
                }}
                const scalar gmin = returnReduce(dmin, minOp<scalar>());
                tipPoint_ = (dmin <= gmin + SMALL) ? imin : -1;
                if (gmin > 1e-6)
                {{
                    WarningInFunction << "point A is not a mesh point, distance " << gmin << endl;
                }}
            }}
            vector p(Zero);
            label cnt = 0;
            if (tipPoint_ >= 0) {{ p = pts[tipPoint_]; cnt = 1; }}
            reduce(p, sumOp<vector>());
            reduce(cnt, sumOp<label>());
            if (cnt > 0) {{ p /= scalar(cnt); }}
            Info<< "TipA " << mesh().time().value() << " "
                << p.x() - A0_.x() << " " << p.y() - A0_.y() << endl;
        #}};
    }}
}}
""" + FOOTER
    (case/"system"/"controlDict").write_text(s)


def write_solution(case, n_outer):
    s = header("fvSolution") + f"""solvers
{{
    "cellDisplacement.*"
    {{
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-9;
        relTol          0;
    }}

    "(p|pcorr)"
    {{
        solver          GAMG;
        smoother        DICGaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }}

    "(p|pcorr)Final"
    {{
        $p;
        relTol          0;
    }}

    "U.*"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }}
}}

PIMPLE
{{
    momentumPredictor       yes;
    nOuterCorrectors        {n_outer};
    nCorrectors             2;
    nNonOrthogonalCorrectors 1;
    correctPhi              yes;
    moveMeshOuterCorrectors yes;
}}

relaxationFactors
{{
    equations
    {{
        ".*"            1;
    }}
}}
""" + FOOTER
    (case/"system"/"fvSolution").write_text(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", choices=sorted(CASES), default="FSI2")
    ap.add_argument("--n-segments", type=int, default=20)
    ap.add_argument("--cells-per-segment", type=int, default=6)
    ap.add_argument("--box-segments", type=int, default=5)
    ap.add_argument("--eta", type=float, default=0.0)
    ap.add_argument("--n-outer", type=int, default=None)
    ap.add_argument("--relax", type=float, default=None)
    a = ap.parse_args()

    prm = CASES[a.case]
    n_outer = a.n_outer or {"FSI2": 6, "FSI3": 12}[a.case]
    relax = a.relax or {"FSI2": 0.5, "FSI3": 0.1}[a.case]

    case = CASE_DIR
    for d in ("system", "constant", "0.orig"):
        (case/d).mkdir(exist_ok=True)

    bm, errs, info = build_mesh(a.n_segments, a.cells_per_segment, a.box_segments)
    if errs:
        print("blockMesh topology errors:")
        for e in errs:
            print("  -", e)
        raise SystemExit(1)
    bm.write(case/"system"/"blockMeshDict")
    write_topo_and_patch(case, a.n_segments, info["dx"])
    EI, m = write_dynamic_mesh_dict(case, a.n_segments, prm, a.eta, relax)
    write_U(case, prm["Um"])
    write_transport(case)
    write_control(case, prm, a.n_segments)
    write_solution(case, n_outer)

    rhoA = prm["rho_s"]*H_FLAG
    bl = np.array([1.875104068711961, 4.694091132974175, 7.854757438237613])
    f_dry = bl**2/(2*np.pi*L_FLAG**2)*math.sqrt(EI/rhoA)
    summary = f"""Case: {a.case}
Flag: L = {L_FLAG:.6f} m, h = {H_FLAG} m, clamp x = {X_ROOT:.6f}, point A = ({X_TIP}, {YC})
Structure: rho_s = {prm['rho_s']}, E = {prm['E']:.4g}, nu_s = {NU_S} -> EI (plane strain) = {EI:.6g} N*m^2/m, rho*A = {rhoA:.6g} kg/m
Fluid: rho_f = {prm['rho_f']:.4g}, nu_f = {NU_F}, U_mean = {prm['Um']}  (Re = {prm['Um']*2*R/NU_F:.0f})
Segments: N = {a.n_segments}, dx = {info['dx']:.6f} m, segment mass = {m:.6g} kg, eta = {a.eta} s
Natural frequencies of the dry cantilever (continuum, Euler-Bernoulli): {', '.join(f'{f:.4f}' for f in f_dry)} Hz
Mesh: {info['cells']} cells; cell size along the flag {info['hf']:.5f} m; n_r={info['n_r']}, n_q={info['n_q']}, n_ur={info['n_ur']},
      n_tip={info['n_tip']}, n_flag={info['n_fx']}, n_wake={info['n_w']} (wake grading {info['g_w']:.2f}), grading towards the flag {info['g_f']:.2f}
Coupling: nOuterCorrectors = {n_outer}, accelerationRelaxation = {relax}, dt = {prm['dt']}, endTime = {prm['T']}
Reference (Turek & Hron 2006): {REFERENCE[a.case]}
"""
    (case/"caseInfo.txt").write_text(summary)
    print("blockMesh topology checked: block orientation, edge and face matching, patches")


if __name__ == "__main__":
    main()
