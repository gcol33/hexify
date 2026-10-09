"""Every zone of a DGGAL grid at one level: its textZIRS and uint64ZIRS
identifiers and its centroid in WGS84 geodetic longitude and latitude.

Usage: python dggal_zone_ids.py <DGGRS> <level> <out.csv>
       e.g. python dggal_zone_ids.py ISEA3H 6 isea3h_6.csv
Needs the dggal Python package (pip install dggal).
"""

import csv
import sys

import dggal


def main():
    name, level, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    app = dggal.Application(appGlobals=globals())
    dggal.pydggal_setup(app)
    grid = getattr(dggal, name)()
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["text", "uint64", "lon", "lat"])
        for z in grid.listZones(level, dggal.wholeWorld):
            c = grid.getZoneWGS84Centroid(z)
            w.writerow([grid.getZoneTextID(z), str(int(z)), repr(float(c.lon)),
                        repr(float(c.lat))])


if __name__ == "__main__":
    main()
