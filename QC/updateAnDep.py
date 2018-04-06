#!/usr/bin/env python
import sys
import argparse
from obsProperties import observation,obsQuality
from db import sqlite
from firstGuessCheck import calculateDeparture
import sqlite3

def main():
    parser = argparse.ArgumentParser(description='Update analysis departures')

    parser.add_argument('-d','--db', dest="sqlitefile", default=None, type=str,help='Name of database to write',required=True)
    parser.add_argument('-a','--an', dest="anfile", default=None, type=str,help='Name of analysis file',required=True)
    parser.add_argument('-v', '--var', dest="var", default=None, type=str, help='Name of variable',required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    conn = sqlite3.connect(args.sqlitefile)
    c = conn.cursor()
    c.execute("SELECT lon, lat,name, value, fgdep, elevation, status FROM obs WHERE name == '"+args.var+"'")
    data = c.fetchall()
    observations=list()
    for row in data:
        lon = row[0]
        lat = row[1]
        name= row[2]
        value= row[3]
        fgdep= row[4]
        elevation=row[5]
        obsQ=obsQuality(15)
        observations+=[observation(name,lon,lat,value,obsQ,elevation=elevation,fgdep=fgdep)]

    conn.close()

    db=sqlite(args.sqlitefile)
    calculateDeparture(args.anfile,args.var,observations,"an")
    db.update(observations)


if __name__ == "__main__":
    main()
