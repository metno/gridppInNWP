#!/usr/bin/env python
import os
import sys
import numpy as np
from math import isnan
import argparse
import netCDF4
from obsProperties import observation,obsQuality
from firstGuessCheck import calculateDeparture,FGCheck
from db import sqlite

def writeParamFile(filename,observations,method):

    lons = list()
    lats = list()
    values = list()
    elevs = list()
    cis = list()
    for obs in observations:
        if obs.status <= 0:
            if not isnan(obs.elevation):
                if obs.elevation > 7000:
                    print "MT Everest",obs.lon,obs.lat,obs.value,obs.elevation
                else:
                    lons += [obs.lon]
                    lats += [obs.lat]
                    values += [obs.value]
                    elevs += [obs.elevation]
                    cis += [obs.ci]

    S = len(lons)
    file = netCDF4.Dataset(filename, 'w')
    file.createDimension("time")
    file.createDimension("y", S)
    file.createDimension("x", 1)
    Nc = 1
    if method == "oi":
        Nc = 2
    file.createDimension("coefficient", Nc)

    vTime = file.createVariable("time", "f8", ("time"))
    vLat = file.createVariable("latitude", "f8", ("y", "x"))
    vLon = file.createVariable("longitude", "f8", ("y", "x"))
    vElev = file.createVariable("altitude", "f8", ("y", "x"))
    vCoeff = file.createVariable("coefficients", "f8", ("time", "y", "x", "coefficient"))

    if S > 0:

        vLat[:] = lats
        vLon[:] = lons
        vElev[:] = elevs
        print("Writing %d stations" % S)
        if method == "oi":
            vCoeff[0, :, 0, 0] = np.array(values)
            vCoeff[0, :, 0, 1] = np.array(cis)
        elif method == "override":
            vCoeff[0, :, 0, 0] = np.array(values)
        else:
            raise NotImplementedError
    else:
        print("No stations written to file")
    file.close()


def readInput(var,file,lonrangeInput,latrangeInput,keepInput,providersInput,delim,default_ci,override_ci,add,multiply):

    latrange = [-180, 180]
    lonrange = [-180, 180]
    if latrangeInput is not None:
        latrange = [float(x) for x in latrangeInput.split(',')]
    if lonrangeInput is not None:
        lonrange = [float(x) for x in lonrangeInput.split(',')]

    #print file
    observations=list()
    ifile = open(file, 'r')
    header = ifile.readline().strip().split(delim)
    Ilat = header.index("lat")
    Ilon = header.index("lon")
    Ielev = header.index("elev")
    Ici = None
    if "rep" in header:
        Ici = header.index("rep")
    Ivalue = header.index("value")
    keep = None
    if keepInput is not None:
        keep = [int(q) for q in keepInput.split(',')]
        if "dqc" not in header:
            print("File '%s' missing 'dqc' column. Cannot select based on dqc." % file)
            return observations
        Idqc = header.index("dqc")

    providers = None
    if providersInput is not None:
        providers = [int(q) for q in providersInput.split(',')]
        if "prid" not in header:
            print("File '%s' missing 'prid' column. Cannot select based on provider." % file)
            return observations
        Iprovider = header.index("prid")

    for line in ifile:
        words = line.strip().split(delim)
        lat = float(words[Ilat])
        lon = float(words[Ilon])
        oQ = obsQuality(15)

        status=-1
        if keep is not None:
            dqc = int(words[Idqc])
            status=dqc
            #print status
            #if dqc not in keep:
            #    continue
        if providers is not None:
            provider = int(words[Iprovider])
            if provider not in providers:
                continue
        elev=float(words[Ielev])
        value=(float(words[Ivalue])+add)*multiply
        ci_value = default_ci
        if override_ci is not None:
            ci_value = override_ci
        elif Ici is not None:
            try:
                ci_value = float(words[Ici])
            except Exception as e:
                ci_value = default_ci

        #print lon,lat,value,add,multiply
        obs=observation(var, lon, lat, value, oQ, elevation=elev, ci=ci_value)
        obs.status=status
        if lat > latrange[0] and lat < latrange[1] and lon > lonrange[0] and lon < lonrange[1]:
            observations += [obs]
        else:
            obs.status=11
            observations += [obs]


    return observations

def main():
    parser = argparse.ArgumentParser(description='Create parameter file for calibration methods (currently OI and override) in gridpp.i The program reads TITAN output and puts observed and CI values into the parameter file. If CI (rep column) is missing in the input, then a default value is used.')
    parser.add_argument('files', help='Input file', nargs="*")
    parser.add_argument('-o', help='Output file', dest="ofilename", required=True)
    parser.add_argument('-k', metavar="FLAGS", help='Only keep rows with these DQC flags. Comma-separated list accepted.', dest="keep")
    parser.add_argument('-m', default="oi", help='Write parameters for this method (default oi).', dest="method", choices=["oi", "override"])
    parser.add_argument('-p', metavar="PROVIDERS", help='Only keep rows with these providers.  Comma-separated list accepted.', dest="providers")
    parser.add_argument('-latrange', help='Only keep stations within these latitudes (min,max)', dest="latrange")
    parser.add_argument('-lonrange', help='Only keep stations within these longitudes (min,max)', dest="lonrange")
    parser.add_argument('--add', type=float, default=0, help='Add this value to the values (e.g. to change units)', dest="add")
    parser.add_argument('--delim', default=';', help='Delimiter separating columns', dest="delim")
    parser.add_argument('--multiply', type=float, default=1, help='Multiply this value to the values (e.g. to change units)', dest="multiply")
    parser.add_argument('--override_ci', type=float, help='Override CI values with this value (if -m oi)', dest="override_ci")
    parser.add_argument('--default_ci', default=1, type=float, help='Use this value if CI is missing (if -m oi; default 1)', dest="default_ci")
    parser.add_argument('--debug', help='Show debug information', action="store_true")
    parser.add_argument('-d','--db', dest="sqlitefile", default=None, type=str,help='Name of database to write')
    parser.add_argument('-f','--fg', dest="fgfile", default=None, type=str,help='Name of first guess file')
    parser.add_argument('-v', '--var', dest="var", default=None, type=str, help='Name of variable')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Setup database
    db=None
    if args.sqlitefile is not None:
        db = sqlite(args.sqlitefile)

    print args.files


    delim=[",",";"]
    keep=[None,"0"]
    files=args.files
    add=[0,273.15]
    multiply=[1,1]
    fgcheck=[True,False]

    if any(fgcheck) and args.fgfile is None:
        print "You must specify a first guess file"
        sys.exit(1)

    observations = list()
    for file in range(0,len(files)):

        print files[file]
        # Read observations
        observations_read=readInput(args.var,files[file], args.lonrange, args.latrange, keep[file],
                           args.providers, delim[file], args.default_ci, args.override_ci,
                           add[file], multiply[file])


        # First guess check
        if args.fgfile is not None:
            calculateDeparture(args.fgfile, args.var, observations_read, "fg")

            if fgcheck[file]:
                print "Performing First Guess Check"
                FGCheck(observations_read,db)

        observations=observations+observations_read


    if args.sqlitefile is not None:
        db.update(observations)

    # Write parameter file
    writeParamFile(args.ofilename, observations, method=args.method)



if __name__ == "__main__":
    main()
