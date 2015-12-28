CFLAGS = -g -O3 -Wall -I$(OPT64)/include -I$(CUDA)/include
CFLAGS+= -DUSE_VEGAS_STATUSID -msse2 # -DUSE_SSE_TRANSPOSE -msse2
#CFLAGS+= -DTRANSPOSE_NETBLKS_ON_GPU
NVCCFLAGS = -O3 -g -arch compute_35 -code=sm_35 -I$(CUDA)/include -I$(OPT64)/include
PROGS = check_guppi_databuf check_guppi_status clean_guppi_shmem \
	test_udp_recv test_psrfits test_psrfits_read \
	fix_psrfits_polyco \
	raw_quant raw_read raw_head
OBJS  = guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
	guppi_params.o guppi_time.o guppi_thread_args.o \
	write_psrfits.o read_psrfits.o misc_utils.o \
	fold.o polyco.o hget.o hput.o sla.o downsample.o \
	dedisperse_utils.o quantization.o privilege_management.o
THREAD_PROGS = test_net_thread guppi_daq guppi_daq_fold
THREAD_OBJS  = guppi_net_thread.o guppi_net_thread_codd.o guppi_pktsock_thread_codd.o \
	       guppi_rawdisk_thread.o \
	       guppi_psrfits_thread.o guppi_fold_thread.o \
	       guppi_null_thread.o
CUDA_OBJS = guppi_dedisp_thread.o guppi_dedisp_ds_thread.o \
	    fold_gpu.o downsample_gpu.o dedisperse_gpu.o
LIBS = -L$(OPT64)/lib -lcfitsio -L$(PRESTO)/lib -lsla -lm -lpthread -lcap -lhashpipe
CUDA_LIBS = -L$(CUDA)/lib64 -lcufft -lcuda -lcudart -lrt -lm
all: $(PROGS) $(THREAD_PROGS) guppi_daq psrfits_subband \
	test_dedisp_thread guppi_daq_dedisp guppi_daq_server \
	libguppi_status.so
clean:
	rm -f $(PROGS) $(THREAD_PROGS) psrfits.tgz *~ *.o test_psrfits_0*.fits *.ptx
tags:
	ctags -R .
.PHONY: all clean tags
INSTALL_DIR = /usr/local/guppi_daq
BIN_DIR = $(INSTALL_DIR)/bin
INC_DIR = $(INSTALL_DIR)/include
LIB_DIR = $(INSTALL_DIR)/lib
install: $(PROGS) $(THREAD_PROGS) guppi_daq psrfits_subband test_dedisp_thread guppi_daq_dedisp guppi_daq_server
	mkdir -p $(BIN_DIR) && \
	cp -f $(PROGS) $(THREAD_PROGS) psrfits_subband test_dedisp_thread guppi_daq_dedisp guppi_daq_server $(BIN_DIR) && \
	mkdir -p $(INC_DIR) && \
	cp -f guppi_error.h guppi_status.h fitshead.h $(INC_DIR) && \
	mkdir -p $(LIB_DIR) && \
	cp -f libguppi_status.so $(LIB_DIR)
psrfits.tgz: psrfits.h read_psrfits.c write_psrfits.c polyco.c polyco.h \
	guppi_PSRFITS_v3.4_fold_template.txt \
	guppi_PSRFITS_v3.4_search_template.txt
	tar cvzf $@ $^
libguppi_status.so: guppi_status.h guppi_status.c guppi_error.h guppi_error.c hput.c hget.c
	gcc -shared -fPIC -pthread -DUSE_VEGAS_STATUSID -o $@ $^ -lm
find_dropped_blocks: find_dropped_blocks.o
	$(CC) $(CFLAGS) $< -o $@ -L$(OPT64)/lib -lcfitsio -lm
psrfits_subband: psrfits_subband.c psrfits_subband_cmd.o $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ psrfits_subband_cmd.o $(OBJS) $(LIBS)
%.o : %.cu
	nvcc -c $(NVCCFLAGS) $< -o $@

test_dedisp_thread: test_dedisp_thread.c $(THREAD_OBJS) $(OBJS) $(CUDA_OBJS)
	$(CC) $(CFLAGS) $(CUDA_CFLAGS) $< -o $@ $(THREAD_OBJS) \
		$(CUDA_OBJS) $(OBJS) $(LIBS) $(CUDA_LIBS)
guppi_daq_dedisp: guppi_daq_dedisp.c $(THREAD_OBJS) $(OBJS) $(CUDA_OBJS)
	$(CC) $(CFLAGS) $(CUDA_CFLAGS) $< -o $@ $(THREAD_OBJS) \
		$(CUDA_OBJS) $(OBJS) $(LIBS) $(CUDA_LIBS)
guppi_daq_server: guppi_daq_server.c $(THREAD_OBJS) $(OBJS) $(CUDA_OBJS)
	$(CC) $(CFLAGS) $(CUDA_CFLAGS) $< -o $@ $(THREAD_OBJS) \
		$(CUDA_OBJS) $(OBJS) $(LIBS) $(CUDA_LIBS)
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS) $(THREAD_LIBS)
$(THREAD_PROGS): $$@.c $(THREAD_OBJS) $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(THREAD_OBJS) $(OBJS) $(LIBS)
