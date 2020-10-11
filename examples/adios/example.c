#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define NUM_NETCDF_FLAVORS 5
#define NDIM 2
#define DIM_LEN 2

int main(int argc, char* argv[])
{
    int my_rank;
    int ntasks;
    int format[NUM_NETCDF_FLAVORS];
    int niotasks = 4;
    int ioproc_stride = 4;
    int ioproc_start = 0;
    int dimid[NDIM];
    int iosysid;
    int ncid;
    int varid;
    PIO_Offset start[NDIM];
    PIO_Offset count[NDIM];
    double buffer[NDIM];
    char filename[NC_MAX_NAME + 1];
    int num_flavors = 0;

#ifdef TIMING
    GPTLinitialize();
#endif

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if (ntasks != 16)
    {
        if (my_rank == 0)
            fprintf(stderr, "Number of processors must be 16!\n");

        MPI_Finalize();
        return 1;
    }

    PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                        PIO_REARR_SUBSET, &iosysid);

#ifdef _PNETCDF
    format[num_flavors++] = PIO_IOTYPE_PNETCDF;
#endif
#ifdef _NETCDF
    format[num_flavors++] = PIO_IOTYPE_NETCDF;
#endif
#ifdef _NETCDF4
    format[num_flavors++] = PIO_IOTYPE_NETCDF4C;
    format[num_flavors++] = PIO_IOTYPE_NETCDF4P;
#endif
#ifdef _ADIOS2
    format[num_flavors++] = PIO_IOTYPE_ADIOS;
#endif

    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        sprintf(filename, "test_put_vars_%d.nc", fmt);

        PIOc_createfile(iosysid, &ncid, &(format[fmt]), filename, PIO_CLOBBER);

        PIOc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]);

        PIOc_def_dim(ncid, "nbnd", DIM_LEN, &dimid[1]);
        PIOc_def_var(ncid, "time_bnds", PIO_DOUBLE, NDIM, dimid, &varid);
        PIOc_enddef(ncid);

        start[1] = 0;
        count[1] = 2;

        start[0] = 0;
        count[0] = 1;
        buffer[0] = 0.0;
        buffer[1] = 0.0;
        PIOc_put_vars_double(ncid, varid, start, count, NULL, buffer);

        start[0] = 1;
        count[0] = 1;
        buffer[0] = 0.0;
        buffer[1] = 0.0833333333333333;
        PIOc_put_vars_double(ncid, varid, start, count, NULL, buffer);

        start[0] = 2;
        count[0] = 1;
        buffer[0] = 0.0833333333333333;
        buffer[1] = 0.25;
        PIOc_put_vars_double(ncid, varid, start, count, NULL, buffer);

        start[0] = 3;
        count[0] = 1;
        buffer[0] = 0.25;
        buffer[1] = 0.333333333333333;
        PIOc_put_vars_double(ncid, varid, start, count, NULL, buffer);

        PIOc_closefile(ncid);
    }

    PIOc_finalize(iosysid);

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}

