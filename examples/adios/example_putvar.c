#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define NUM_NETCDF_FLAVORS 2
#define NDIM 1
#define DIM_LEN 16

int main(int argc, char* argv[])
{
    int my_rank;
    int ntasks;
    int format[NUM_NETCDF_FLAVORS];
    int niotasks = 2;
    int ioproc_stride = 2;
    int ioproc_start = 0;
    int dimid;
    int iosysid;
    int ncid;
    int varid;
    PIO_Offset start[NDIM];
    PIO_Offset count[NDIM];
    int buffer[4] = {0};
    char filename[NC_MAX_NAME + 1];
    int num_flavors = 0;

#ifdef TIMING
    GPTLinitialize();
#endif

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if (ntasks != 4)
    {
        if (my_rank == 0)
            fprintf(stderr, "Number of processors must be 4!\n");

        MPI_Finalize();
        return 1;
    }

    PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                        PIO_REARR_SUBSET, &iosysid);

#ifdef _PNETCDF
    format[num_flavors++] = PIO_IOTYPE_PNETCDF;
#endif

#ifdef _ADIOS2
    format[num_flavors++] = PIO_IOTYPE_ADIOS;
#endif

	for (int ii=0;ii<4;ii++) {
		buffer[ii] = ii;
	}

    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        sprintf(filename, "test_put_vara_%d.nc", fmt);

        PIOc_createfile(iosysid, &ncid, &(format[fmt]), filename, PIO_CLOBBER);

        PIOc_def_dim(ncid, "row", DIM_LEN, &dimid);
        PIOc_def_var(ncid, "foo", PIO_INT, NDIM, &dimid, &varid);
        PIOc_enddef(ncid);

        count[0] = 4;

        start[0] = 0;
        PIOc_put_vara_int(ncid, varid, start, count, buffer);

        start[0] = 4;
        PIOc_put_vara_int(ncid, varid, start, count, buffer);

        start[0] = 8;
        PIOc_put_vara_int(ncid, varid, start, count, buffer);

        start[0] = 12;
        PIOc_put_vara_int(ncid, varid, start, count, buffer);

        PIOc_closefile(ncid);
    }

    PIOc_finalize(iosysid);

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}

