"""
GiD-free test of the $$GHOST_STABILISATION block of createCPS.tcl.

The block between the comment "# GHOST STABILISATION:" and the line that writes the value line is extracted from
the problem type and evaluated in a plain Tcl interpreter (tkinter.Tcl) with mocks of the GiD and tDOM calls it uses:
``$stageNode selectNodes <xpath>`` returns the value of the field named in the xpath, ``GiD_WriteCalculationFile puts``
appends to a list, and ``=`` is ``format``. Field values come from the defaults of xml/CalculationData.xml unless a
case overrides them; a missing field returns "" as for models created before the field existed.

Usage: python test_createcps_ghost.py   (exit code 0 if all cases pass)

>>> ghost_line({})[:6]
['0', '0.05', '0.05', '0.6', '1', '2']
"""
import os
import re
import sys
import tkinter
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
GID = os.path.join(HERE, "..", "..", "GiD_Problemtype", "Anura3D_2025.gid")
TCL = os.path.join(GID, "scripts", "createCPS.tcl")
XML = os.path.join(GID, "xml", "CalculationData.xml")

MOCKS = r"""
proc stageNodeMock {method xpath} {
    if {![regexp {@n="([^"]+)"\]/@v\)$} $xpath -> name]} {error "unexpected xpath $xpath"}
    if {[dict exists $::fields $name]} {return [dict get $::fields $name]}
    return ""
}
set stageNode stageNodeMock
proc GiD_WriteCalculationFile {cmd text} {lappend ::out $text}
proc = {fmt args} {format $fmt {*}$args}
"""


def block_source():
    """Return the Tcl source of the ghost block of createCPS.tcl."""
    with open(TCL, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    start = next(i for i, s in enumerate(lines) if s.strip().startswith("# GHOST STABILISATION:"))
    end = next(i for i in range(start, len(lines)) if "$Ghost_flag $Ghost_gK" in lines[i])
    return "\n".join(lines[start:end + 1])


def xml_defaults():
    """Return {field name: default value} of the STABILISATION container of CalculationData.xml."""
    root = ET.parse(XML).getroot()
    cont = next(c for c in root.iter("container") if c.get("n") == "STABILISATION")
    return {v.get("n"): v.get("v") for v in cont.iter("value")}


def ghost_line(fields):
    """Evaluate the block with the given field values and return the tokens of the written value line."""
    tcl = tkinter.Tcl()
    tcl.eval(MOCKS)
    tcl.eval("set ::out {}")
    tcl.eval("set ::fields [dict create]")
    for k, v in fields.items():
        tcl.call("dict", "set", "::fields", k, v)
    tcl.eval(block_source())
    out = tcl.splitlist(tcl.eval("set ::out"))
    assert out[0] == "$$GHOST_STABILISATION", out
    return out[1].split()


def main():
    d = xml_defaults()
    on = dict(d, GHOST_STABILISATION="apply ghost stabilisation")
    cases = [
        ("model without STABILISATION container (old model)", {}, "0 0.05 0.05 0.6 1 2 0 0 0.25 0 0 0"),
        ("XML defaults, ghost off", d, "0 0.05 0.05 0.6 1 2 0 0 0.25 0 0 0"),
        ("XML defaults, ghost on", on, "1 0.05 0.05 0.6 1 2 0 0 0.25 0 0 0"),
        ("old model with ghost on, new fields missing",
         {k: v for k, v in on.items() if k in ("GHOST_STABILISATION", "ghost_preset", "ghost_fill_threshold",
                                               "ghost_mass_solver", "ghost_jacobi_iterations")},
         "1 0.05 0.05 0.6 1 2 0 0 0.25 0 0 0"),
        ("node-based, S4a on 0.2, border file, diagonal solver, 3 sweeps",
         dict(on, ghost_border_detection="node-based topological", ghost_nodal_mass_net="on",
              ghost_nodal_mass_threshold="0.2", ghost_diagnostics="extended with border file",
              ghost_mass_solver="diagonal", ghost_jacobi_iterations="3"),
         "1 0.05 0.05 0.6 0 3 1 1 0.2 0 2 0"),
        ("extended diagnostics, user penalties", dict(on, ghost_preset="user", ghost_gamma_K="0.02",
                                                      ghost_gamma_M="0.03", ghost_diagnostics="extended"),
         "1 0.02 0.03 0.6 1 2 0 0 0.25 0 1 0"),
        ("symmetric F1 scheme, coded bound", dict(on, ghost_mass_scheme="symmetric (F1)"),
         "1 0.05 0.05 0.6 1 2 0 0 0.25 1 0 0"),
        ("symmetric F1 scheme with the mass-independent bound, node-based",
         dict(on, ghost_mass_scheme="symmetric (F1)", ghost_time_step_bound="mass-independent (F1 only)",
              ghost_border_detection="node-based topological"),
         "1 0.05 0.05 0.6 1 2 1 0 0.25 1 0 1"),
    ]
    fails = 0
    for name, fields, expected in cases:
        got = " ".join(ghost_line(fields))
        ok = got == expected
        fails += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}: '{got}'" + ("" if ok else f" expected '{expected}'"))
    for key in ("ghost_border_detection", "ghost_nodal_mass_net", "ghost_nodal_mass_threshold", "ghost_mass_scheme",
                "ghost_diagnostics", "ghost_time_step_bound"):
        ok = key in d
        fails += not ok
        print(f"{'PASS' if ok else 'FAIL'}  XML field {key} present (default '{d.get(key)}')")
    print(f"SUMMARY failed {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
