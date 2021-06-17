import time,sqlite3
import pandas as pd
import sqlalchemy as sqla
from . import config as cfg
__version__ = '0.1.2'

class CypTyper:
    '''This is really just a db interface.  Typing is done in the db'''
    def __init__(self,alleles,variants,dbFile,dbStore=False):
        self.alleles  = alleles
        self.variants = variants
        self.dbFile   = dbFile
        self.db       = self._openDB(dbFile)
        self.dbStore  = dbStore
        self._import2DB()

    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        if not self.dbStore:
            uuids = tuple(self.alleles.index)
            for table in ['variantTable','alleleTable']:
                sql   = f'''DELETE
                          FROM {cfg.database[table]}
                          WHERE uuid in {uuids}'''
                try:
                    con = self.db.connect()
                    con.execute(sql)
                except (sqla.exc.OperationalError,sqlite3.OperationalError) as e:
                    raise Typer_Error(f'Unable to delete records in {self.dbFile}. DB error as follows:\n\n{e}')            
    
    def getSummary(self,table,failed=False):
        uuids    = tuple(self.alleles.index)
        sql = f'SELECT * FROM {table}'
        try:
            result = pd.read_sql(sql,con=self.db)
            if 'uuid' in result.columns: #full table
                query = 'uuid in @uuids'
                if not failed:
                    query += ' and status == "passed"'
                return result.query(query) \
                             .iloc[:,1:] \
                             .sort_values(cfg.typer['sortColumns'])
            else:
                return result  #need some better filtering here
        except (sqla.exc.OperationalError,sqlite3.OperationalError) as e:
            raise Typer_Error(f'Unable to execute query in {self.dbFile}. DB error as follows:\n\n{e}')

    def _openDB(self,dbFile):
        return sqla.create_engine(f'sqlite:///{dbFile}', echo=False)

    def _import2DB(self):
        tables   = [cfg.database[t] for t in ['alleleTable','variantTable']]
        support  = cfg.database['supportField'] 
        maxTries = cfg.database['maxTries']
        for pydf,sqldf in zip([self.alleles,self.variants],tables):
            pbaa = pydf.reset_index()\
                       .reindex(columns=list(cfg.tableMap[sqldf].values()))
            pbaa.columns = list(cfg.tableMap[sqldf].keys())
            #make sure support dtyep is "string"
            if support in pbaa.columns:
                pbaa[support] = pbaa[support].astype(str)
            tries = 0
            while tries < maxTries:
                try:
                    pbaa.set_index('uuid').to_sql(sqldf, con=self.db, if_exists='append')
                    break
                #except sqlite3.OperationalError as e:
                except sqla.exc.OperationalError as e:
                    tries += 1
                    print(f'WARNING: sqlite3 import error try #{tries} of {maxTries}')
                    if tries == maxTries:
                        raise Typer_Error(f'Unable to import {pydf.source.unique()} to {self.dbFile}')
                    else:
                        time.sleep(1)
        return None

class Typer_Error(Exception):
    pass
