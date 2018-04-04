class obsQuality(object):

    def __init__(self,threshold):
        self.threshold=threshold


class observation():
        def __init__(self,name,lon,lat,value,obsQ,stid=None,elevation=None,fgdep=None,ci=1.):
            self.name=name
            self.lon=lon
            self.lat=lat
            self.value=value
            self.stid=stid
            self.elevation=elevation
            self.ci=ci
            self.fgdep=fgdep
            self.interpolation=None
            self.ind_x=-1
            self.ind_y=-1
            self.obsQ=obsQ
            self.status=-1
            #print("Constructed observation")

        """ First guess departure """
        @property
        def fgdep(self):
            return self.fgdep

        @fgdep.setter
        def fgdep(self, fgdep):
            self.fgdep = fgdep

        @property
        def ind_x(self):
            return self.ind_x

        @ind_x.setter
        def ind_x(self, ind_x):
            self.ind_x = ind_x

        @property
        def ind_y(self):
            return self.ind_y

        @ind_y.setter
        def ind_y(self, ind_y):
            self.ind_y = ind_y

        @property
        def fg(self):
            if self.fgdep == None:
                return np.nan
            else:
                return self.value+self.fgdep

        def getStatusText(self):
            statusDesc={
               -1:"Not checked",
                0:"OK",
                1:"Domain check",
                2:"First Guess check"
            }
            stat="Not defined"
            if self.status in statusDesc: stat=statusDesc[self.status]
            return stat

        def isOK(self):
            #print self.obsQ.threshold,self.fgdep,self.value
            if self.status != 0: return False
            if abs(self.obsQ.threshold) > abs(self.fgdep):
                return True
            else:
                return False
