import netCDF4
import numpy as np
from scipy.interpolate import NearestNDInterpolator
import abc
import scipy.spatial.qhull as qhull
import sys

def FGCheck(observations,db=None):

    for obs in observations:
        if obs.status == 0 and not obs.isOK(): obs.status=2

    if db is not None: db.update(observations)


def readFile(fname,var):
    nc = netCDF4.Dataset(fname, 'r')
    var_lons = nc["longitude"]
    var_lats = nc["latitude"]
    field = nc[var]
    return nc,var_lons,var_lats,field

def closeFile(nc):
    nc.close()

def interpolate(method,observations,var_lons,var_lats,field):

    # Find first guess for point (nearest neighbour)

    lons = list()
    lats = list()
    for obs in observations:
        lons += [obs.lon]
        lats += [obs.lat]

    print("Interpolating "+str(len(observations))+" observations." )
    print("Input file dimensions: var_lons="+str(var_lons.shape)+" var_lats="+str(var_lons.shape))

    if method == "nearest":
        nearest = NearestNeighbour(lons, lats, var_lons, var_lats)

        i=0
        for obs in observations:
            ind_y = nearest.index[i][0]
            ind_x = nearest.index[i][1]
            obs.ind_x = ind_x
            obs.ind_y = ind_y
            i=i+1

def calculateDeparture(fname,var,observations,dep,interp=True):

    nc,var_lons, var_lats, field=readFile(fname,var)
    method="nearest"
    if interp: interpolate(method, observations, var_lons, var_lats, field)

    if method == "nearest":

        i=0
        for obs in observations:

            #print field.shape
            if obs.ind_x > 0 and obs.ind_y > 0:
                interpolated_value=field[0][obs.ind_y][obs.ind_x]
                depval= interpolated_value - (obs.value+273.15)

            else:
                depval = np.nan

            #print depval
            if dep == "fg":
                obs.fgdep = depval
            elif dep == "an":
                if obs.dqc == 0:
                    obs.andep = depval
                else:
                    obs.andep=np.nan
            else:
                print "Departure type " + str(dep) + " not defined!"
                sys.exit(1)
            i=i+1

    closeFile(nc)

class Interpolation(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self,type,nx,ny,var_lons,var_lats):
        self.type=type
        self.nx=nx
        self.ny=ny
        self.var_lons=var_lons
        self.var_lats=var_lats

    def distance(self,lon1, lat1, lon2, lat2):
        """
        Computes the great circle distance between two points using the
        haversine formula. Values can be vectors.
        """
        # Convert from degrees to radians
        pi = 3.14159265
        lon1 = lon1 * 2 * pi / 360
        lat1 = lat1 * 2 * pi / 360
        lon2 = lon2 * 2 * pi / 360
        lat2 = lat2 * 2 * pi / 360
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
        c = 2 * np.arcsin(np.sqrt(a))
        distance = 6.367e6 * c
        return distance

class NearestNeighbour(Interpolation):

    def __init__(self,interpolated_lons,interpolated_lats,var_lons,var_lats):
        nx = var_lons.shape[1]
        ny = var_lats.shape[0]
        self.index=self.create_index(interpolated_lons, interpolated_lats, var_lons,var_lats)
        super(NearestNeighbour, self).__init__("nearest",nx,ny,var_lons,var_lats)


    def create_index(self,interpolated_lons, interpolated_lats, var_lons,var_lats):

        dim_x = var_lons.shape[1]
        dim_y = var_lats.shape[0]
        npoints = len(interpolated_lons)

        lons_vec = np.reshape(var_lons,dim_x*dim_y)
        lats_vec = np.reshape(var_lats,dim_x*dim_y)

        points=np.empty([dim_x*dim_y,2])
        points[:,0]=lats_vec
        points[:,1]=lons_vec


        values_vec=np.empty([dim_x*dim_y])
        x = []
        y = []
        ii = 0
        for j in range(0, dim_y):
            for i in range(0, dim_x):
                values_vec[ii] = ii
                x.append(i)
                y.append(j)
                ii = ii + 1

        #print("Interpolating..." + str(len(interpolated_lons)) + " points")
        nn = NearestNDInterpolator(points, values_vec)

        # Set max distance as sanity
        if len(lons_vec) > 1 and len(lats_vec) > 1:
            max_distance=1.5*self.distance(lons_vec[0],lats_vec[0],lons_vec[1],lats_vec[1])
        else:
            print("You only have one point is your input field!")
            sys.exit(1)

        grid_points = []
        for n in range(0, npoints):
            ii = nn(interpolated_lats[n], interpolated_lons[n])
            #print ii
            ii = int(round(float(ii),0))
            i = x[ii]
            j = y[ii]
            #print ii,i,j,dim_x,dim_y,interpolated_lons[n],interpolated_lats[n],lons_vec[ii],lats_vec[ii]
            dist=self.distance(interpolated_lons[n],interpolated_lats[n],lons_vec[[ii]],lats_vec[ii])
            #print dist,max_distance
            if dist > max_distance:
                #print("Point is too far away from nearest point: "+str(dist)+" Max distance="+str(max_distance))
                i=-1
                j=-1
                #sys.exit(1)

            #sys.exit(1)
            grid_points.append([j, i])
        return grid_points

"""
Not tested yet
"""
class Linear(Interpolation):

    def __init__(self,int_lons,int_lats,var_lons,var_lats):

        nx,ny=self.setup_weights(int_lons, int_lats, var_lons,var_lats)
        super(Linear, self).__init__("linear",nx,ny,var_lons,var_lats)


    def setup_weights(self, int_lons, int_lats, var_lons,var_lats):
        print("Setup weights for linear interpolation")

        # Posistions in input file
        lons=var_lons
        lats=var_lats
        xy = np.zeros([lons.shape[0] * lons.shape[1], 2])
        xy[:, 0] = lons.flatten()
        xy[:, 1] = lats.flatten()

        # Target positions
        [Xi, Yi] = [np.asarray(int_lons), np.asarray(int_lats)]
        uv = np.zeros([len(Xi), 2])
        uv[:, 0] = Xi.flatten()
        uv[:, 1] = Yi.flatten()

        # Setup the weights
        self.interp_weights(xy, uv)
        return lons.shape[1],lons.shape[0]

    def interp_weights(self,xy, uv, d=2):
        tri = qhull.Delaunay(xy)
        simplex = tri.find_simplex(uv)
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uv - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
        self.vtx=vertices
        self.wts=np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    def interpolate(self,values):
        values=values.flatten()
        return np.einsum('nj,nj->n', np.take(values, self.vtx), self.wts)
