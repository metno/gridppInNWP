#!/usr/bin/env python
import sys
import numpy as np
import argparse
from firstGuessCheck import calculateDeparture
from db import sqlite


class observation():

    def __init__(self,var,words,Ilon,Ilat,Ielev,Ivalue,Ici,Idqc,Iprovider,Isct):
        if var is None:
            self.name="unknown"
        else:
            self.name=var
        self.lon=float(words[Ilon])
        self.lat=float(words[Ilat])
        self.elev=wordValue(words,Ielev)
        self.value = wordValue(words, Ivalue)
        self.ci = wordValue(words, Ici)
        self.dqc = wordValue(words, Idqc)
        self.provider = wordValue(words, Iprovider)
        self.sct= wordValue(words,Isct)
        self.fgdep=np.nan
        self.andep=np.nan


def wordValue(words,ind):
    val=np.nan
    if ind >= 0:
        if words[ind] == 'NA':
            val=np.nan
        else:
            val=float(words[ind])
    return val

def hdrIndex(header,col):
    ind=-1
    if col in header:
        ind=header.index(col)
    else:
        print("WARNING: File is missing '"+col+"' column.")
    return ind

def readInput(var,file,delim):

    #print file
    observations=list()
    ifile = open(file, 'r')
    header = ifile.readline().strip().split(delim)
    Ilat = hdrIndex(header,"lat")
    Ilon = hdrIndex(header,"lon")
    Ielev = hdrIndex(header,"elev")
    Ivalue = hdrIndex(header, "value")
    Ici = hdrIndex(header,"rep")
    Idqc = hdrIndex(header,"dqc")
    Isct = hdrIndex(header, "sct")
    Iprovider = hdrIndex(header,"prid")

    l=1
    for line in ifile:
        words = line.strip().split(delim)
        if Ilon >= 0 and Ilat >= 0:
            observations += [observation(var,words,Ilon,Ilat,Ielev,Ivalue,Ici,Idqc,Iprovider,Isct)]
        else:
            print "Longitude at latitude not found in line "+str(l)
        l=l+1


    return observations

def main():
    parser = argparse.ArgumentParser(description='Convert TITAN output to SQLite')
    parser.add_argument('-i', dest="input", help='Input file',required=True )
    parser.add_argument('-o', dest="output", help='Output SQLite file', required=True)
    parser.add_argument('--debug', dest="debug", help='Show debug information', action="store_true")
    parser.add_argument('--delim', dest="delim", type=str, help='Delimiter', default=";")
    parser.add_argument('-f','--fg', dest="fgfile", default=None, type=str,help='Name of first guess file')
    parser.add_argument('-v', '--var', dest="var", default=None, type=str, help='Name of variable')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    db = sqlite(args.output)
    observations=readInput(args.var,args.input,args.delim)


    # First guess check
    if args.fgfile is not None:
        calculateDeparture(args.fgfile, args.var, observations, "fg")

    db.update(observations)


if __name__ == "__main__":
    main()
