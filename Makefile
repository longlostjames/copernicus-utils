CC=gcc
LIBS= -lnetcdf -lm
CFLAGS= -Wall -O3

# Uncomment this for Darwin/PPC (Mac OS X)
#NETCDF_LIB_PATH=/sw/lib
#NETCDF_INCLUDE_PATH=/sw/include

proc-copernicus : proc-copernicus.c
	$(CC) $(CFLAGS) -I$(NETCDF_INCLUDE_PATH) -L$(NETCDF_LIB_PATH) $(LIBS) -o proc-copernicus proc-copernicus.c

clean :
	rm proc-copernicus
