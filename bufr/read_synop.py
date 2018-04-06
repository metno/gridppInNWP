#!/usr/bin/env python2.7

from __future__ import print_function
import traceback
import sys
import numpy as np
from cartopy import crs as ccrs
import matplotlib.pyplot as plt
 
from eccodes import *

VERBOSE=False

def example(INPUT):
    # open bufr file
    f = open(INPUT)
 
    # define the keys to be printed
    keys = [
        #'blockNumber',
        #'stationNumber',
        'latitude',
        'longitude',
        'airTemperatureAt2M',
        'relativeHumidity',
        'totalSnowDepth'
    ]
 
    # The cloud information is stored in several blocks in the
    # SYNOP message and the same key means a different thing in different
    # parts of the message. In this example we will read the first
    # cloud block introduced by the key
    # verticalSignificanceSurfaceObservations=1.
    # We know that this is the first occurrence of the keys we want to
    # read so in the list above we used the # (occurrence) operator
    # accordingly.
 
    cnt = 0

    lons=list()
    lats=list()
    t2m=list()
    sd=list()
    rh2m=list() 
    # loop for the messages in the file
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
        t=np.nan
        r=np.nan
        s=np.nan
        for key in keys:
            try:
                val=codes_get(bufr, key)
                #if val != CODES_MISSING_DOUBLE:
                #    print('  %s: %s' % (key,val))
                if val == CODES_MISSING_DOUBLE: val=np.nan
                if key == "latitude": lat=val
                if key == "longitude": lon=val
                if key == "airTemperatureAt2M": t=val
                if key == "relativeHumidity": r=val
                if key == "totalSnowDepth": s=val
            except CodesInternalError as err:
                dummy=1
                print('Error with key="%s" : %s' % (key, err.msg))
 
        lats.append(lat)
        lons.append(lon)
        t2m.append(t)
        rh2m.append(r)
        sd.append(s)

        cnt += 1
 
        # delete handle
        codes_release(bufr)
 
    # close the file
    f.close()
 
    print(len(lons))
    print(len(lats))
    print(len(t2m))
    print(len(sd))
    print(len(rh2m))

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
    print(t2m)
    print(rh2m)
    for t in range(0,len(rh2m)):
        if rh2m[t] > 100: rh2m[t]=100.
    print(rh2m)
    print(sd)
    plt.scatter(lons, lats, c=t2m, transform=ccrs.PlateCarree(),cmap="Paired")
    plt.colorbar()
    plt.show()
    plt.close()
    ax = plt.axes(projection=proj)
    ax.coastlines(resolution="10m")
    plt.scatter(lons, lats, c=rh2m, transform=ccrs.PlateCarree(),cmap="Paired")
    plt.colorbar()
    plt.show()
    plt.close()
    ax = plt.axes(projection=proj)
    ax.coastlines(resolution="10m")
    plt.scatter(lons, lats, c=sd, transform=ccrs.PlateCarree(),cmap="Paired")
    plt.colorbar()
    plt.show()

def main(INPUT):
    try:
        example(INPUT)
    except CodesInternalError as err:
        if VERBOSE:
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
 
        return 1
 
 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))


