import netCDF4 as nc
import eccodes as ec
import yaml
from pyproj import Proj
import numpy as np
import sys
import os
from argparse import ArgumentParser
from datetime import datetime,timedelta
import calendar

class converter():
    def __init__(self,name,var,config,ua_name,sfx_name=None):
        self.name=name
        self.var=var
        self.config=config
        self.ua_name=ua_name
        self.sfx_name=sfx_name
        #print("Constructed converter")

    def get_grib_definitions(self,var=None):
        fname = self.ua_name
        if var == None:
            dict=self.config[self.var]["converter"][self.name]
        else:
            dict=self.config[self.var]["converter"][self.name][var]

        if "ftype" in dict:
            if dict["ftype"] == "sfx": fname = self.sfx_name
            if dict["ftype"] == "ua": fname = self.ua_name

        if fname == None:
            print("The file is not set for "+self.var+" converter: "+self.name)
            sys.exit(1)
        if not os.path.isfile(fname):
            print("The file "+fname+" does not exist!")
            sys.exit(1)

        par = dict["par"]
        lev = dict["lev"]
        typ = dict["typ"]
        tri = dict["tri"]
        return fname,par,lev,typ,tri

    def read_field(self):
        if self.name == "none":
            fname,par,lev,typ,tri=self.get_grib_definitions()
            var=variable(fname,par,lev,typ,tri)
            lons, lats, X, Y,dt,field=var.read_field()
            return lons,lats,X,Y,dt,field
        elif self.name == "swe2sd":
            fname, par, lev, typ, tri = self.get_grib_definitions("swe")
            print fname,par,lev,type,tri
            var = variable(fname, par, lev, typ, tri)
            lons, lats, X, Y,dt,swe = var.read_field()
            fname, par, lev, typ, tri = self.get_grib_definitions("rho")
            var = variable(fname, par, lev, typ, tri)
            lons, lats, X, Y,dt,rho = var.read_field()
            field=np.divide(swe,rho)
            return lons,lats,X,Y,dt,field
        elif self.name == "sweclim":
            fname, par, lev, typ, tri = self.get_grib_definitions()
            var = variable(fname, par, lev, typ, tri)
            lons, lats, X, Y, dt,field = var.read_field()
            rhoclim={"01":222.,"02":233.,"03":240.,"04":278.,"05":212.,"06":312.,"07":312.,"08":143.,"09":143.,"10":161.,"11":182.,"12":213.}
            month=dt.strftime("%m")
            if month in rhoclim:
                field=np.divide(field,rhoclim[month])
                return lons, lats, X, Y, dt,field
            else:
                print("Could not found climatological mean for month "+str(month))
                sys.exit(1)

        elif self.name == "phi2m":
            fname, par, lev, typ, tri = self.get_grib_definitions()
            var = variable(fname, par, lev, typ, tri)
            lons, lats, X, Y, dt,field = var.read_field()
            field=np.divide(field,9.81)
            return lons, lats, X, Y, dt,field
        else:
            print("Converter not "+self.name+" not defined")
            sys.exit(1)


class variable():

    def __init__(self,fname,par,lev,typ,tri):
        self.fname=fname
        self.par=par
        self.lev=lev
        self.typ=typ
        self.tri=tri
        #print("Constructed variable")


    def read_field(self):

        geography = ["bitmapPresent",
                     "Nx",
                     "Ny",
                     "latitudeOfFirstGridPointInDegrees",
                     "longitudeOfFirstGridPointInDegrees",
                     "LoVInDegrees",
                     "DxInMetres",
                     "DyInMetres",
                     "iScansNegatively",
                     "jScansPositively",
                     "jPointsAreConsecutive",
                     "Latin1InDegrees",
                     "LaDInDegrees",
                     "Latin2InDegrees",
                     "latitudeOfSouthernPoleInDegrees",
                     "longitudeOfSouthernPoleInDegrees",
                     "gridType"
        ]

        if self.fname == None or not os.path.isfile(self.fname):
            print("The file "+str(self.fname)+" does not exist!")
            sys.exit(1)

        print("Reading file: "+self.fname)
        f = open(self.fname, "r")
        while 1:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break


            par = ec.codes_get(gid, "indicatorOfParameter")
            lev = ec.codes_get(gid, "level")
            typ = ec.codes_get(gid, "indicatorOfTypeOfLevel")
            tri = ec.codes_get(gid, "timeRangeIndicator")


            #print("Search::", w_par, w_lev, w_typ, w_tri)
            if self.par == par and self.lev == lev and self.typ == typ and self.tri == tri:
                print("Found:", self.par, self.lev, self.typ, self.tri)
                geo = {}
                for key in geography:
                    try:
                        geo.update({key: ec.codes_get(gid, key)})
                    except ec.CodesInternalError as err:
                        print('Error with key="%s" : %s' % (key, err.msg))

                print('There are %d values, average is %f, min is %f, max is %f' % (
                    ec.codes_get_size(gid, 'values'),
                    ec.codes_get(gid, 'average'),
                    ec.codes_get(gid, 'min'),
                    ec.codes_get(gid, 'max')
                ))

                # Date/time
                d=ec.codes_get(gid,"validityDate")
                t=ec.codes_get(gid,"validityTime")
                h=int(t)/100
                m=t%h
                s=(h*3600)+(m*60)
                date=datetime.strptime(str(d),"%Y%m%d")
                time=timedelta(seconds=s)
                dt=date+time

                # Missing values
                mv=None
                try:
                    mv=ec.codes_get(gid,"missingValue")
                except:
                    print("Field does not contanin missing values")

                if geo["gridType"].lower() == "lambert":
                    values = ec.codes_get_values(gid)
                    nx = geo["Nx"]
                    ny = geo["Ny"]

                    lonCenter = geo["LoVInDegrees"]
                    latCenter = geo["LaDInDegrees"]
                    latRef = geo["Latin2InDegrees"]
                    lon0 = geo["longitudeOfFirstGridPointInDegrees"]
                    lat0 = geo["latitudeOfFirstGridPointInDegrees"]
                    dx = geo["DxInMetres"]
                    dy = geo["DyInMetres"]

                    proj4_string = "+proj=lcc +lat_0=" + str(latCenter) + " +lon_0=" + str(lonCenter) + " +lat_1=" + str(
                    latRef) + " +lat_2=" + str(latRef) + " +no_defs +units=m +R=6.371e+06"
                    proj4 = Proj(proj4_string)

                    x0, y0 = proj4(lon0, lat0)
                    x0 = int(round(x0))
                    y0 = int(round(y0))
                    field = np.empty([nx, ny])
                    lons = np.empty([nx, ny])
                    lats = np.empty([nx, ny])
                    X = np.arange(x0, x0 + (nx * dx), dx)
                    Y = np.arange(y0, y0 + (ny * dy), dy)
                    ii = 0
                    for i in range(0, nx):
                        for j in range(0, ny):
                            field[i, j] = values[ii]
                            lons[i, j], lats[i, j] = proj4(X[i], Y[j], inverse=True)
                            # print i,j,lons[i, j], lats[i, j]
                            ii = ii + 1

                    if mv is not None:
                        field[field==mv]=np.nan
                    ec.codes_release(gid)
                    f.close()
                    return(lons, lats, X,Y,dt,field)
                else:
                    print(geo["gridType"] + " not implemented yet!")

            ec.codes_release(gid)
        f.close()

def create_template(vars,nx,ny,fname="raw.nc"):

    fg = nc.Dataset(fname, "w")
    fg.createDimension("y", ny)
    fg.createDimension("x", nx)
    fg.createDimension("time", 1)
    fg.createVariable("time", "f8", ("time",))
    fg["time"].long_name = "time"
    fg["time"].standard_name = "time"
    fg["time"].units = "seconds since 1970-01-01 00:00:00 +00:00"
    fg.createVariable("longitude", "f8", ("y", "x"))
    fg["longitude"].units = "degree_east"
    fg["longitude"].long_name = "longitude"
    fg["longitude"].standard_name = "longitude"
    fg.createVariable("latitude", "f8", ("y", "x"))
    fg["latitude"].units = "degree_north"
    fg["latitude"].long_name = "latitude"
    fg["latitude"].standard_name = "latitude"
    fg.createVariable("x", "f4", ("x"))
    fg["x"].long_name = "x-coordinate in Cartesian system"
    fg["x"].standard_name = "projection_x_coordinate"
    fg["x"].units = "m"
    fg.createVariable("y", "f4", ("y"))
    fg["y"].long_name = "y-coordinate in Cartesian system"
    fg["y"].standard_name = "projection_y_coordinate"
    fg["y"].units = "m"

    standard_name = {"air_temperature_2m": "air_temperature", "relative_humidity_2m": "relative_humidity",
                     "altitude": "altitude", "surface_snow_thickness": "surface_snow_thickness"}
    long_name = {"air_temperature_2m": "Screen level temperature (T2M)",
                 "relative_humidity_2m": "Screen level relative humidity (RH2M)",
                 "altitude": "Altitude",
                 "surface_snow_thickness": "Surface snow thickness"}
    units = {"air_temperature_2m": "K", "relative_humidity_2m": "1", "altitude": "m",
             "surface_snow_thickness": "m"}
    fillvalue = {"air_temperature_2m": "9.96921e+36", "relative_humidity_2m": "9.96921e+36",
               "altitude": "9.96921e+36", "surface_snow_thickness": "9.96921e+36"}

    for var in vars:
        fg.createVariable(var, "f4", ("time", "y", "x"), fill_value=fillvalue[var])
        fg[var].long_name = long_name[var]
        fg[var].standard_name = standard_name[var]
        fg[var].units = units[var]

    return fg


parser = ArgumentParser(description="Create first guess file for gridpp from grib")
parser.add_argument('ua_gribfile', type=str, help="Upper air grib file")
parser.add_argument('-sfx', type=str, default=None, help="SURFEX grib file",nargs="?")
parser.add_argument('--sd_converter', type=str, default="none", help="",nargs="?",choices=["sweclim","swe2sd"])
parser.add_argument('--altitude_converter', type=str, default="phi2m", help="",nargs="?",choices=["none","phi2m"])
parser.add_argument('--config','-c', type=str, help="YAML config file", default="grib_codes.yaml", nargs="?")
args = parser.parse_args(sys.argv[1:])


grib_codes = yaml.load(open(args.config))
ftype="ua"
if "ftype" in grib_codes: ftype=grib_codes[ftype]
vars=["air_temperature_2m","relative_humidity_2m","surface_snow_thickness","altitude"]
#vars=["air_temperature_2m"]

first=True
for var in vars:
    ua_name = args.ua_gribfile
    sfx_name=args.sfx
    convertName="none"
    if var == "surface_snow_thickness": convertName=args.sd_converter
    if var == "altitude": convertName = args.altitude_converter

    lons,lats,X,Y,dt,field=converter(convertName,var,grib_codes,ua_name,sfx_name).read_field()

    # Create file
    if first: fg=fg=create_template(vars,X.shape[0],Y.shape[0])
    if var == "altitude":
        field[field<0]=0

    fg[var][:]=field
    if first:
        epoch=calendar.timegm(dt.timetuple())
        fg["time"][:]=epoch
        fg["longitude"][:]=np.transpose(lons)
        fg["latitude"][:]=np.transpose(lats)
        fg["x"][:]=X
        fg["y"][:]=Y

    first=False
fg.close()
