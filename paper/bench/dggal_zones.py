"""Every zone of a DGGAL grid at one level: its centroid and its vertices.

DGGAL (https://github.com/ecere/dggal) reads coordinates as WGS84 geodetic
and maps them to the authalic sphere before projecting, so latitudes are
written back on the authalic sphere, where hexify's grids live.

Usage: python dggal_zones.py <DGGRS> <level> <out.csv>
       e.g. python dggal_zones.py IVEA3H 4 ivea3h_4.csv
Needs the dggal Python package (pip install dggal).
"""

import csv
import math
import sys

import dggal

# WGS84
A_E = 6378137.0
F = 1.0 / 298.257223563
E2 = F * (2.0 - F)
E = math.sqrt(E2)


def q(phi):
    s = math.sin(phi)
    return (1.0 - E2) * (s / (1.0 - E2 * s * s)
                         - math.log((1.0 - E * s) / (1.0 + E * s)) / (2.0 * E))


Q_P = q(math.pi / 2.0)


def authalic_deg(lat_deg):
    """Authalic latitude of a geodetic one, in degrees."""
    return math.degrees(math.asin(max(-1.0, min(1.0, q(math.radians(lat_deg)) / Q_P))))


def main():
    name, level, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    app = dggal.Application(appGlobals=globals())
    dggal.pydggal_setup(app)
    grid = getattr(dggal, name)()
    zones = grid.listZones(level, dggal.wholeWorld)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["zone", "kind", "k", "lon", "lat"])
        for z in zones:
            zid = grid.getZoneTextID(z)
            c = grid.getZoneWGS84Centroid(z)
            w.writerow([zid, "centroid", 0, repr(float(c.lon)), repr(authalic_deg(float(c.lat)))])
            for k, v in enumerate(grid.getZoneWGS84Vertices(z)):
                w.writerow([zid, "vertex", k + 1, repr(float(v.lon)),
                            repr(authalic_deg(float(v.lat)))])


if __name__ == "__main__":
    main()
