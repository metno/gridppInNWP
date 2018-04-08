#!/usr/bin/env python2.7

from __future__ import print_function
import traceback
import argparse
import numpy as np
CAN_PLOT=True
try:
    from cartopy import crs as ccrs
    import matplotlib.pyplot as plt
except:
    CAN_PLOT=False

import sys
from eccodes import codes_set,codes_bufr_new_from_file,CodesInternalError,codes_release,codes_get,CODES_MISSING_DOUBLE,CODES_MISSING_LONG


class observation():
    def __init__(self, lon, lat, value, elevation):
        self.lon = lon
        self.lat = lat
        self.value = value
        self.elevation = elevation

def readBufrFile(bufrFile,var,lonrange,latrange):
    # open bufr file
    f = open(bufrFile)
 
    # define the keys to be printed
    keys = [
        #'blockNumber',
        #'stationNumber',
        'latitude',
        'longitude',
        'heightOfStation'
    ]
    #'airTemperatureAt2M',
    #    'relativeHumidity',
    #    'totalSnowDepth'
    keys.append(var)
 
    # The cloud information is stored in several blocks in the
    # SYNOP message and the same key means a different thing in different
    # parts of the message. In this example we will read the first
    # cloud block introduced by the key
    # verticalSignificanceSurfaceObservations=1.
    # We know that this is the first occurrence of the keys we want to
    # read so in the list above we used the # (occurrence) operator
    # accordingly.


    print("Reading "+bufrFile)
    print("Looking for keys: "+str(keys))
    cnt = 0
    observations=list()

    # loop for the messages in the file
    not_found = 0
    while 1:
        # get handle for message
        bufr = codes_bufr_new_from_file(f)
        if bufr is None:
            break
 
        #print("message: %s" % cnt)
 
        # we need to instruct ecCodes to expand all the descriptors
        # i.e. unpack the data values
        codes_set(bufr, 'unpack', 1)
 
        # print the values for the selected keys from the message
        
        lat=np.nan
        lon=np.nan
        value=np.nan
        elev=np.nan
        for key in keys:
            try:
                val=codes_get(bufr, key)
                #if val != CODES_MISSING_DOUBLE:
                #    print('  %s: %s' % (key,val))
                if val == CODES_MISSING_DOUBLE or val == CODES_MISSING_LONG: val=np.nan
                if key == "latitude": lat=val
                if key == "longitude": lon=val
                if key == "heightOfStation": elev=val
                if key == var:
                    value = val
                    if var == "relativeHumidity":
                        if value > 100: values = 100.
                    elif var == "airTemperatureAt2M":
                        value=value-273.15

            except CodesInternalError as err:
                if key == var:
                    not_found=not_found+1
                #print('Error with key="%s" : %s' % (key, err.msg))

        if lat > latrange[0] and lat < latrange[1] and lon > lonrange[0] and lon < lonrange[1]:
            if not np.isnan(value):
                observations.append(observation(lon,lat,value,elev))

        cnt += 1
 
        # delete handle
        codes_release(bufr)

    print("Found "+str(len(observations))+"/"+str(cnt))
    print("Not encoded for "+str(var)+": "+str(not_found))
    # close the file
    f.close()
    return observations

def writeObs(fname,observations):

    fh=open(fname,"w")
    fh.write("lon;lat;elev;value\n")
    for obs in observations:
        fh.write(str(obs.lat)+";"+str(obs.lon)+";"+str(obs.elevation)+";"+str(obs.value)+"\n")
    fh.close()


def two_floats(value):
    values = value.split(",")
    #print(values)
    if len(values) != 2:
        raise argparse.ArgumentError
    values = map(float, values)
    return values

def main(argv):

    parser = argparse.ArgumentParser(description='Parse airTemperatureAt2M,relativeHumidity or totalSnowDepth from bufr file')
    parser.add_argument('-f',dest="input",help='Inputfile file', required=True)
    parser.add_argument('-o', dest="output", help='Output file', required=True)
    parser.add_argument('-v',dest="var",help='variabel',required=True)
    parser.add_argument('--plot',dest="plot",default=False,action='store_true')
    parser.add_argument('-latrange', help='Only keep stations within these latitudes (min,max)', dest="latrange", action = 'store', type=two_floats,default=[-180,180])
    parser.add_argument('-lonrange', help='Only keep stations within these longitudes (min,max)',dest="lonrange", action = 'store', type=two_floats,default=[-180,180])
    parser.add_argument('--debug', dest="debug", help='Debug', action='store')

    if len(argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    try:
        observations=readBufrFile(args.input,args.var,args.lonrange,args.latrange)
    except CodesInternalError as err:
        if args.debug:
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
 
        return 1


    writeObs(args.output,observations)

    if args.plot and CAN_PLOT:
        lons=list()
        lats=list()
        values=list()
        elevs=list()
        for obs in observations:
            lats.append(obs.lat)
            lons.append(obs.lon)
            values.append(obs.value)
            elevs.append(obs.elevation)

        proj = ccrs.LambertConformal(central_longitude=15., central_latitude=63., standard_parallels=[63.])
        # Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r,
        # GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r,
        # Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r,
        # Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r,
        # Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r,
        # YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg,
        # brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r,
        # gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow,
        # gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray,
        # gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r,
        # ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral,
        # spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
        plt.close()
        ax = plt.axes(projection=proj)
        ax.coastlines(resolution="10m")
        plt.scatter(lons, lats, c=values, transform=ccrs.PlateCarree(),cmap="coolwarm")
        plt.title(args.var)
        plt.colorbar()
        plt.show()
        plt.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))


