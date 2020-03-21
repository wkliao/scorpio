#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define NUM_NETCDF_FLAVORS 2
#define NDIM 2
#define DIM_LEN 1
#define VAR_NUM 10

int main(int argc, char* argv[])
{
    int my_rank;
    int ntasks;
    int format[NUM_NETCDF_FLAVORS];
    int niotasks;
    int ioproc_stride = 1;
    int ioproc_start = 0;
    int dimid[NDIM];
    PIO_Offset elements_per_pe;
    int dim_len[1] = {DIM_LEN};
    int iosysid;
    int ncid;
    int varid;
    int wr_iodesc;
    PIO_Offset *compdof = NULL;
    int *buffer;
    char filename[NC_MAX_NAME + 1];
    char varname[NC_MAX_NAME + 1];
    int num_flavors = 0;

#ifdef TIMING
    GPTLinitialize();
#endif

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if (ntasks != 1)
    {
        if (my_rank == 0)
            fprintf(stderr, "Number of processors must be 1!\n");

        MPI_Finalize();
        return 1;
    }

    niotasks = ntasks;
    PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                        PIO_REARR_SUBSET, &iosysid);

    elements_per_pe = DIM_LEN / ntasks;
    compdof = malloc(elements_per_pe * sizeof(PIO_Offset));

    /* Decomp : [1] */
    compdof[0] = 1;

    PIOc_InitDecomp(iosysid, PIO_INT, 1, dim_len, (PIO_Offset)elements_per_pe,
                    compdof, &wr_iodesc, NULL, NULL, NULL);

    buffer = malloc(elements_per_pe * sizeof(int));
    for (int i = 0; i < elements_per_pe; i++)
        buffer[i] = compdof[i];

    free(compdof);

#ifdef _PNETCDF
    format[num_flavors++] = PIO_IOTYPE_PNETCDF;
#endif

#ifdef _ADIOS2
    format[num_flavors++] = PIO_IOTYPE_ADIOS;
#endif

    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        sprintf(filename, "test_ioid_%d.nc", fmt);

        PIOc_createfile(iosysid, &ncid, &(format[fmt]), filename, PIO_CLOBBER);

        PIOc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]);
        PIOc_def_dim(ncid, "lndgrid", DIM_LEN, &dimid[1]);

        for (int i = 0; i < VAR_NUM; i++)
        {
            sprintf(varname, "var_id_%d", i);
            PIOc_def_var(ncid, varname, PIO_INT, NDIM, dimid, &varid);
        }

        PIOc_enddef(ncid);

        for (varid = 0; varid < VAR_NUM; varid++)
        {
            PIOc_setframe(ncid, varid, 0);
            PIOc_write_darray(ncid, varid, wr_iodesc, elements_per_pe, buffer, NULL);
        }

        PIOc_closefile(ncid);
    }

    free(buffer);

    PIOc_freedecomp(iosysid, wr_iodesc);

    PIOc_finalize(iosysid);

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}
