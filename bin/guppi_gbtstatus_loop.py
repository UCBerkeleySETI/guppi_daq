#!/usr/bin/env python
import sys
import time
import atexit
import signal
from guppi_daq.guppi_utils import *

# Attach to status shared mem
g = guppi_status()

def close_and_exit(signum, frame):
    print 'caught signal %d, cleaning and exiting' % signum
    g.close_gbtstatus()
    sys.exit(0)

atexit.register(g.close_gbtstatus)
signal.signal(signal.SIGINT, close_and_exit)
signal.signal(signal.SIGTERM, close_and_exit)

while (1):
    try:
        g.lock()
        g.read(lock=False)
        g.update_with_gbtstatus()
        g.write(lock=False)
        g.unlock()
    except:
        g.unlock()
    time.sleep(1)
