#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A universal mathematical model of non-photochemical quenching

Copyright (C) 2015-2016  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""

import sqlite3
from datetime import datetime
import numpy as np


class DB:

    """ Connect to the data base and retrieve single experiment results for further analysis """

    def __init__(self, dblocation='DataBase/paperdata.db'):
        self.conn = sqlite3.connect(dblocation, detect_types=sqlite3.PARSE_DECLTYPES)
        self.cursor = self.conn.cursor()

    def __enter__(self):
        return self

    def __exit__(self):
        self.conn.close()

    def retrieve_file_names(self, searchfields):
        """
        :param searchfields: a dictionary of search criteria, e.g. {'lightIntensity':900,'darkDuration':60}
        :return: list of file names from the data base matching the searchFields requirements
        """

        clauses = []
        for k, v in searchfields.items():
            if k == 'specie':
                clauses.append('('+k+' like \''+v+'\')')
            else:
                clauses.append('('+k+' == '+str(v)+')')

        whereclause = ' AND '.join(clauses)

        c = self.cursor
        fnames = []
        for row in c.execute('SELECT DISTINCT(filename) FROM lightmemory WHERE '+whereclause):
            fnames.append(row)

        return fnames

    def retrieve_data_sets(self, searchfields,
                           retrievefields=['expid', 'time', 'Ft', 'Fm', 'PAR', 'qN', 'Yield', 'ETR'],
                           infofields=['specie', 'darkduration', 'lightintensity', 'replicate']):
        """
        :param searchfields: dictionary with the search criteria. Cannot be empty
        :param retrievefields: a list of database fields stored in the resulting dataSet object
        :param infofields: database fields with useful information characterising the whole data set,
        which are therefore expected to be identical in all entries in the data set and thus only
        stored as additional information once
        :return: data sets according to the search criteria
        """

        fnames = self.retrieve_file_names(searchfields)

        c = self.cursor
        datasets = []

        selectclause = ', '.join(retrievefields + infofields)

        for fn in fnames:
            d = []
            for row in c.execute('SELECT '+selectclause+' FROM lightmemory WHERE filename like ?', fn):
                d.append(row)
            
            optinfo = {'filename': fn}
            for i in range(len(infofields)):
                optinfo[infofields[i]] = d[0][len(retrievefields) + i]

            datasets.append(DataSet(d, retrievefields, optinfo))

        return datasets


class DataSet:

    """ Methods performed on single datasets

        Attributes:
            expid: number of the measurement during the experiment
            filename: name of the .dat file from which the data were collected into the data base
            lightintensity: intensity of the light during both phases [in micro Einsteins per m2 per s]
            darkduration: duration of the in-between relaxation phase [in minutes]
            numericFields: returns the list of fields with float values
            replicate: number of the replicate [1,2,3]
            specie: name of the specie ['Arabidopsis', 'Pothos']

            ETR: electron transport rate defined by: Yield ∗ PAR ∗ 0.5 ∗ ETR-factor
            Fm: maximal fluorescence
            Ft: steady fluorescence
            PAR: Photosynthetically active radiation at the time shortly preceding the last saturation pulse
            time: exact time of the measurement [hh:mm:ss]
            T: list of time points as the seconds passed from the forst measurements
            Yield: effective quantum yield of energy conversion (first value for the experiments stands for the Fv/Fm)
            qN: coefficient of non-photochemical quenching defined by: (Fm' - F)/(Fm' - F0)
    """

    numericFields = ['expid', 'Ft', 'Fm', 'PAR', 'qN', 'Yield', 'ETR']
 
    def __init__(self, data, fields, optinfo={}):
        """
        a DataSet object contains all the data in an easily accessible way.
        When time is imported, automatically a field 'T' will be
        generated containing the delta-T with respect to the first time point
        """
        
        for i in range(len(fields)):
            dbdata = [data[j][i] for j in range(len(data))]
            if fields[i] in DataSet.numericFields:
                setattr(self, fields[i], np.array(dbdata))
            else:
                setattr(self, fields[i], dbdata)

        for k, v in optinfo.items():
            setattr(self, k, v)

        if hasattr(self, 'time'):
            t0 = datetime.strptime(self.time[0], '%H:%M:%S')
            self.T = np.array([(datetime.strptime(self.time[j], '%H:%M:%S')-t0).seconds
                               for j in range(len(self.time))])

if __name__ == '__main__':
    db = DB()

    # load the data for Arabidopsis into the ds list
    ds = db.retrieve_data_sets({'specie': 'ARA%'})
    print('Loaded all experiments for Arabidopsis thaliana. Type ds[i] for i in range(26) to access specific set')
