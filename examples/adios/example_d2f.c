#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define NUM_NETCDF_FLAVORS 2
#define NDIM 1
#define DIM_LEN 15

int main(int argc, char* argv[])
{
    int my_rank;
    int ntasks;
    int format[NUM_NETCDF_FLAVORS];
    int niotasks = 4;
    int ioproc_stride = 4;
    int ioproc_start = 1;
    int dimid[NDIM];
    int iosysid;
    int ncid;
    int varid;
    PIO_Offset start[NDIM];
    PIO_Offset count[NDIM];
    double buffer[DIM_LEN] = {0.007100635, 0.027925, 0.06225858, 0.1188651, 0.2121934, 0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607, 4.739157, 7.829766, 12.92532, 21.32647, 35.17762};
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
                        PIO_REARR_BOX, &iosysid);

#ifdef _PNETCDF
    format[num_flavors++] = PIO_IOTYPE_PNETCDF;
#endif
#ifdef _ADIOS2
    format[num_flavors++] = PIO_IOTYPE_ADIOS;
#endif

    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        sprintf(filename, "test_put_vars_%d.nc", fmt);

        PIOc_createfile(iosysid, &ncid, &(format[fmt]), filename, PIO_CLOBBER);

        PIOc_def_dim(ncid, "levdcmp", DIM_LEN, &dimid[0]);
        PIOc_def_var(ncid, "levdcmp", PIO_FLOAT, NDIM, dimid, &varid);
        PIOc_enddef(ncid);

        start[0] = 0;
        count[0] = 15;
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
