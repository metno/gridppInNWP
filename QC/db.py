import sqlite3

class sqlite(object):

    def __init__(self,fname):
        self.open(fname)
        self.create_layout()
        #self.close()

    def open(self,fname):
        self.conn = sqlite3.connect(fname)

    def close(self):
        self.conn.close()


    def create_layout(self):
        c = self.conn.cursor()

        # Create table
        #self.ci = wordValue(words, Ici)
        #self.dqc = wordValue(words, Idqc)
        #self.provider = wordValue(words, Iprovider)
        #self.sct = wordValue(words, Isct)
        c.execute("CREATE TABLE IF NOT EXISTS obs (lon, lat, name, value, elev, ci,dqc,provider,sct, fgdep, andep)")

        # Save (commit) the changes
        self.conn.commit()

    def update(self,observations):

        c = self.conn.cursor()
        for obs in observations:
            #print obs.lon,obs.lat,obs.value,obs.fgdep,obs.andep,obs.status
            c.execute("SELECT lon,lat FROM obs WHERE lon == "+str(obs.lon)+" AND lat == "+str(obs.lat)+"")
            if c.fetchone() == None:
                #print "new"
                c.execute("INSERT INTO obs VALUES (?, ?, ?, ?, ?, ?, ?, ? ,?, ?, ?) ",(obs.lon,obs.lat,obs.name,obs.value,obs.elev,obs.ci,obs.dqc,obs.provider,obs.sct,obs.fgdep,obs.andep))
            else:
                #print "update"
                c.execute("UPDATE obs SET lon=?, lat=?, name=?, value=?, fgdep=?, andep=?, elev=?, ci=?, dqc=?, provider=?, sct=?  WHERE lon == "+str(obs.lon)+" AND lat == "+str(obs.lat),(obs.lon, obs.lat, obs.name, obs.value, obs.fgdep, obs.andep,obs.elev,obs.ci,obs.dqc,obs.provider,obs.sct))

        # Save (commit) the changes
        self.conn.commit()