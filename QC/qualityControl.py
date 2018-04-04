#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse
import netCDF4
from obsProperties import observation,obsQuality
from firstGuessCheck import calculateFirstGuess,FGCheck,plotStatus
from db import sqlite

def writeParamFile(filename,observations,method):

    lons = list()
    lats = list()
    values = list()
    elevs = list()
    cis = list()
    for obs in observations:
        if obs.status <= 0:
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


def readInput(var,files,lonrangeInput,latrangeInput,keepInput,providersInput,delim,default_ci,override_ci,add,multiply):

    latrange = [-180, 180]
    lonrange = [-180, 180]
    if latrangeInput is not None:
        latrange = [float(x) for x in latrangeInput.split(',')]
    if lonrangeInput is not None:
        lonrange = [float(x) for x in lonrangeInput.split(',')]

    observations=list()
    for file in files:
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
                continue
            Idqc = header.index("dqc")

        providers = None
        if providersInput is not None:
            providers = [int(q) for q in providersInput.split(',')]
            if "prid" not in header:
                print("File '%s' missing 'prid' column. Cannot select based on provider." % file)
                continue
            Iprovider = header.index("prid")

        for line in ifile:
            words = line.strip().split(delim)
            lat = float(words[Ilat])
            lon = float(words[Ilon])
            oQ = obsQuality(15)

            if keep is not None:
                dqc = int(words[Idqc])
                if dqc not in keep:
                    continue
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

            if lat > latrange[0] and lat < latrange[1] and lon > lonrange[0] and lon < lonrange[1]:
                observations += [observation(var, lon, lat,value,oQ,elevation=elev,ci=ci_value)]
            else:
                obs=observation(var, lon, lat,value,oQ,elevation=elev, ci=ci_value)
                obs.status=1
                observations += [obs]


    return observations

def main():
    parser = argparse.ArgumentParser(description='Create parameter file for calibration methods (currently OI and override) in gridpp.i The program reads TITAN output and puts observed and CI values into the parameter file. If CI (rep column) is missing in the input, then a default value is used.')
    parser.add_argument('files', help='Netatmo input file', nargs="*")
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

    # Need information on:
    # Variable with quality information
    # first guess file
    # sqlite database name
    #

    fgfile = "/home/trygveasp/PycharmProjects/gridppInNWP/firstGuess/raw.nc"
    #fgfile=None

    var = "air_temperature_2m"
    #var = "surface_snow_thickness"
    sqlitefile="test.db"
    #sqlitefile=None

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Setup database
    db=None
    if sqlitefile is not None:
        db = sqlite(sqlitefile)

    # Read observations
    observations=readInput(var, args.files, args.lonrange, args.latrange, args.keep,
                           args.providers, args.delim, args.default_ci, args.override_ci,
                           args.add, args.multiply)



    # First guess check
    if fgfile is not None:
        print "Performing First Guess Check"
        calculateFirstGuess(fgfile,var,observations,db=db)
        FGCheck(observations,db)

    #for obs in observations:
    #    print obs.getStatusText()

    netatmo=["/home/trygveasp/PycharmProjects/gridppInNWP/test/data/20180312T06Z_air_temperature_2m_netatmo.txt"]
    observations_netatmo = readInput(var, netatmo, None, None,None,None,";",1,None,0,1)
    observations=observations+observations_netatmo

    print len(observations)
    if sqlitefile is not None:
        db.update(observations)

    # Write parameter file
    writeParamFile(args.ofilename, observations, method=args.method)



if __name__ == "__main__":
    main()
