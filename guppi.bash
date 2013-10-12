#!/bin/bash
# Set environment variables for GUPPI, bash version
echo "This script is specific to the Shanghai Observatory"
echo "Setting GUPPI_DIR, PATH, PYTHONPATH, LD_LIBRARY_PATH, TEMPO, PRESTO and PGPLOT_DIR for GUPPI..."

export HEADAS=/opt/dibas/pulsar64/src/heasoft-6.6.2/x86_64-unknown-linux-gnu-libc2.3.4
alias ftools=". $HEADAS/headas-init.sh"

PSR64=/opt/dibas/pulsar
OPT64=/opt/dibas/dibaslibs

export CUDA=/usr/local/cuda
export GUPPI=$OPT64/dibas_repos
export GUPPI_DIR=$GUPPI/guppi_daq
export PRESTO=$PSR64/src/presto
export PATH=$PSR64/bin:$PRESTO/bin:$GUPPI_DIR/bin:$GUPPI/bin:$OPT64/bin:$PATH:$CUDA/bin
export PYTHONPATH=$PSR64/lib/python:$PSR64/lib/python/site-packages:$PRESTO/lib/python:$GUPPI/lib/python/site-packages:$GUPPI/lib/python:$GUPPI_DIR/python
export PGPLOT_DIR=$PSR64/pgplot
#export LD_LIBRARY_PATH=$PSR64/lib:$OPT64/lib:$PGPLOT_DIR:$PRESTO/lib:$CUDA/lib64
export LD_LIBRARY_PATH=$PSR64/lib:$PGPLOT_DIR:$PRESTO/lib:$CUDA/lib64
export TEMPO=$PSR64/tempo
