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
        c.execute("CREATE TABLE IF NOT EXISTS obs (lon, lat,name, value, fgdep, andep, elevation, status)")

        # Save (commit) the changes
        self.conn.commit()

    def update(self,observations):

        c = self.conn.cursor()
        for obs in observations:
            #print obs.lon,obs.lat,obs.value,obs.fgdep,obs.andep,obs.status
            c.execute("SELECT lon,lat FROM obs WHERE lon == "+str(obs.lon)+" AND lat == "+str(obs.lat)+"")
            if c.fetchone() == None:
                #print "new"
                c.execute("INSERT INTO obs VALUES (?, ?, ?, ?, ?, ?,?, ?) ",(obs.lon,obs.lat,obs.name,obs.value,obs.fgdep,obs.andep,obs.elevation,obs.status))
            else:
                #print "update"
                c.execute("UPDATE obs SET lon=?, lat=?, name=?, value=?, fgdep=?, andep=?, elevation=?, status=? WHERE lon == "+str(obs.lon)+" AND lat == "+str(obs.lat),(obs.lon, obs.lat, obs.name, obs.value, obs.fgdep, obs.andep,obs.elevation, obs.status))

        # Save (commit) the changes
        self.conn.commit()