from guppi_utils import *

# TODO parse command line args allow manual setting

g = guppi_status()
cur_scan = g["SCAN"]
cur_scan += 1
g.update("SCAN", cur_scan)
g.update("BASENAME", "guppi_test_%s_%04d"%(g['SRC_NAME'], g['SCAN']))
g.write()
