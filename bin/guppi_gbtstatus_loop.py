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

for i in range(5):
    try:
        print 'connecting to gbtstatus database...',
        g.connect_to_gbtstatus()
        print 'ok'
        break
    except Exception as e:
        print 'got exception:'
        print e
        if i < 4:
            print 'will retry in 1 second'
            time.sleep(1)

if not g.is_connected_to_gbtstatus():
    print "unable to connect to gbtstatus database"
    sys.exit(1)

while (1):
    try:
        g.lock()
        g.read(lock=False)
        g.update_with_gbtstatus()
        g.write(lock=False)
    except:
        pass
    finally:
        g.unlock()

    time.sleep(1)
