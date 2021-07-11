#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

/* Number of elements of the local darray data that
   will be handled by each processor. */
#define ELEMENTS_PER_PE 8

/* Length of variables to be put/get */
#define PUT_GET_VAR_LEN 10

int main(int argc, char* argv[])
{
  /* Zero-based rank of processor. */
  int my_rank;

  /* Number of processors involved in current execution. */
  int ntasks;

  /* This example only tests two formats (IO types): PnetCDF and ADIOS */
  int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};

  /* Number of processors that will do IO. In this example we
     will do IO from all processors. */
  int niotasks;

  /* Stride in the MPI rank between IO tasks. Always 1 in this example. */
  const int ioproc_stride = 1;

  /* Zero based rank of first processor to be used for I/O. */
  const int ioproc_start = 0;

  /* The dimension IDs. */
  int dimid_text_var_len;
  int dimid_darray_var_len;
  int dimid_put_get_var_len;

  /* This simple example uses one-dimensional data. */
  const int NDIMS = 1;

  /* Lengths of the global dimensions */
  int gdimlen[NDIMS];

  /* The ID for the parallel I/O system. It is set by
     PIOc_Init_Intracomm(). It references an internal structure
     containing the general IO subsystem data and MPI
     structure. It is passed to PIOc_finalize() to free
     associated resources, after all I/O, but before
     MPI_Finalize is called. */
  int iosysid;

  /* The ncid of the netCDF file created and read in this example. */
  int ncid_write;
  int ncid_read;

  /* The IDs of the netCDF variables in the example file. */
  int varid_dummy_scalar_var_int;
  int varid_dummy_scalar_var_float;

  int varid_dummy_text_var;

  int varid_dummy_darray_var_int;
  int varid_dummy_darray_var_float;
  int varid_dummy_put_get_var_int;
  int varid_dummy_put_get_var_float;

  /* start/count arrays for get/put var */
  PIO_Offset start[NDIMS];
  PIO_Offset count[NDIMS];

  /* The I/O description IDs as passed back by PIOc_InitDecomp()
     and freed in PIOc_freedecomp(). */
  int ioid_int;
  int ioid_double;

  /* Sample data for global attributes. */
  int put_global_att_int_data = -10;
  int get_global_att_int_data = 0;
  float put_global_att_float_data = -9.9;
  float get_global_att_float_data = 0.0;
  char put_global_att_text_data[PIO_MAX_NAME] = "Dummy global attribute string";
  char get_global_att_text_data[PIO_MAX_NAME] = "\0";

  /* Sample data for variable attributes. */
  int put_att_int_data = -100;
  int get_att_int_data = 0;
  float put_att_float_data = -99.99;
  float get_att_float_data = 0.0;
  char put_att_text_data[PIO_MAX_NAME] = "Dummy variable attribute string";
  char get_att_text_data[PIO_MAX_NAME] = "\0";

  /* Sample data for scalar variables. */
  int put_scalar_int_data = -1000;
  int get_scalar_int_data = 0;
  float put_scalar_float_data = -999.999;
  float get_scalar_float_data = 0.0;

  /* Sample string for text variable. */
  char put_text_data[PIO_MAX_NAME] = "Dummy text variable string";
  char get_text_data[PIO_MAX_NAME] = "\0";

  /* Buffers for sample write/read darray data. */
  int write_darray_buffer_int[ELEMENTS_PER_PE];
  int read_darray_buffer_int[ELEMENTS_PER_PE];
  double write_darray_buffer_double[ELEMENTS_PER_PE];
  double read_darray_buffer_double[ELEMENTS_PER_PE];

  /* Buffers for sample put/get var data. */
  int put_var_buffer_int[PUT_GET_VAR_LEN];
  int get_var_buffer_int[PUT_GET_VAR_LEN];
  double put_var_buffer_double[PUT_GET_VAR_LEN];
  double get_var_buffer_double[PUT_GET_VAR_LEN];

  /* A 1-D array which holds the decomposition mapping for this example. */
  PIO_Offset *compmap;

  /* Test filename. */
  char filename[PIO_MAX_NAME];

  float diff_float;
  double diff_double;

  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /* Set lengths of the global dimensions (1D in this example). */
  gdimlen[0] = ELEMENTS_PER_PE * ntasks;

  /* Keep things simple - 1 IO task per MPI process */
  niotasks = ntasks;

  /* Initialize the PIO IO system. This specifies how
     many and which processors are involved in I/O. */
  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                      PIO_REARR_SUBSET, &iosysid); ERR

  /* Describe the decomposition. This is a 1-based array. */
  compmap = malloc(ELEMENTS_PER_PE * sizeof(PIO_Offset));
  for (int i = 0; i < ELEMENTS_PER_PE; i++)
    compmap[i] = my_rank * ELEMENTS_PER_PE + i + 1;

  /* Create the decomposition for this example. */
  ret = PIOc_InitDecomp(iosysid, PIO_INT, NDIMS, gdimlen, ELEMENTS_PER_PE, compmap, &ioid_int, NULL, NULL, NULL); ERR
  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, NDIMS, gdimlen, ELEMENTS_PER_PE, compmap, &ioid_double, NULL, NULL, NULL); ERR
  free(compmap);

  /* Prepare sample data for write buffers and initialize read buffers. */
  for (int i = 0; i < ELEMENTS_PER_PE; i++) {
    write_darray_buffer_int[i] = my_rank;
    write_darray_buffer_double[i] = my_rank * 0.1;
    read_darray_buffer_int[i] = 0;
    read_darray_buffer_double[i] = 0.0;
  }

  for (int i = 0; i < PUT_GET_VAR_LEN; i++) {
    put_var_buffer_int[i] = i + 1;
    put_var_buffer_double[i] = (i + 1) * 0.1;
    get_var_buffer_int[i] = 0;
    get_var_buffer_double[i] = 0.0;
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    /* Create a filename to write. */
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    /* Put some global attributes. */
    ret = PIOc_put_att(ncid_write, PIO_GLOBAL, "dummy_global_att_int", PIO_INT, 1, &put_global_att_int_data); ERR
    ret = PIOc_put_att(ncid_write, PIO_GLOBAL, "dummy_global_att_float", PIO_FLOAT, 1, &put_global_att_float_data); ERR
    ret = PIOc_put_att_text(ncid_write, PIO_GLOBAL, "dummy_global_att_text", strlen(put_global_att_text_data), put_global_att_text_data); ERR

    /* Define some scalar variables. */
    ret = PIOc_def_var(ncid_write, "dummy_scalar_var_int", PIO_INT, 0, NULL, &varid_dummy_scalar_var_int); ERR
    ret = PIOc_def_var(ncid_write, "dummy_scalar_var_float", PIO_FLOAT, 0, NULL, &varid_dummy_scalar_var_float); ERR

    /* Define a text variable. */
    ret = PIOc_def_dim(ncid_write, "text_var_len", PIO_MAX_NAME, &dimid_text_var_len); ERR
    ret = PIOc_def_var(ncid_write, "dummy_text_var", PIO_CHAR, 1, &dimid_text_var_len, &varid_dummy_text_var); ERR

    /* Define some variables for PIOc_write_darray. */
    ret = PIOc_def_dim(ncid_write, "darray_var_len", (PIO_Offset)gdimlen[0], &dimid_darray_var_len); ERR
    ret = PIOc_def_var(ncid_write, "dummy_darray_var_int", PIO_INT, NDIMS, &dimid_darray_var_len, &varid_dummy_darray_var_int); ERR
    ret = PIOc_def_var(ncid_write, "dummy_darray_var_float", PIO_FLOAT, NDIMS, &dimid_darray_var_len, &varid_dummy_darray_var_float); ERR

    /* Put some local attributes for variable dummy_darray_var_int. */
    ret = PIOc_put_att(ncid_write, varid_dummy_darray_var_int, "dummy_att_float", PIO_FLOAT, 1, &put_att_float_data); ERR
    ret = PIOc_put_att(ncid_write, varid_dummy_darray_var_int, "dummy_att_int", PIO_INT, 1, &put_att_int_data); ERR
    ret = PIOc_put_att_text(ncid_write, varid_dummy_darray_var_int, "dummy_att_text", strlen(put_att_text_data), put_att_text_data); ERR

    /* Define some variables for PIOc_put_vars. */
    ret = PIOc_def_dim(ncid_write, "put_get_var_len", PUT_GET_VAR_LEN, &dimid_put_get_var_len); ERR
    ret = PIOc_def_var(ncid_write, "dummy_put_get_var_int", PIO_INT, NDIMS, &dimid_put_get_var_len, &varid_dummy_put_get_var_int); ERR
    ret = PIOc_def_var(ncid_write, "dummy_put_get_var_float", PIO_FLOAT, NDIMS, &dimid_put_get_var_len, &varid_dummy_put_get_var_float); ERR

    ret = PIOc_enddef(ncid_write); ERR

    /* Put some scalar variables. */
    ret = PIOc_put_var_int(ncid_write, varid_dummy_scalar_var_int, &put_scalar_int_data); ERR
    ret = PIOc_put_var_float(ncid_write, varid_dummy_scalar_var_float, &put_scalar_float_data); ERR

    /* Put a text variable. */
    ret = PIOc_put_var_text(ncid_write, varid_dummy_text_var, put_text_data); ERR

    /* Put int type data to an int type variable, type conversions will not be performed. */
    /* Put 1st half data first */
    start[0] = 0;
    count[0] = PUT_GET_VAR_LEN / 2;
    ret = PIOc_put_vars_int(ncid_write, varid_dummy_put_get_var_int, start, count, NULL, put_var_buffer_int); ERR
    /* Put 2nd half data */
    start[0] = PUT_GET_VAR_LEN / 2;
    count[0] = PUT_GET_VAR_LEN / 2;
    ret = PIOc_put_vars_int(ncid_write, varid_dummy_put_get_var_int, start, count, NULL, put_var_buffer_int + (PUT_GET_VAR_LEN / 2)); ERR

    /* Put double type data to a float type variable, type conversions will be performed. */
    /* Put 2nd half data first */
    start[0] = PUT_GET_VAR_LEN / 2;
    count[0] = PUT_GET_VAR_LEN / 2;
    ret = PIOc_put_vars_double(ncid_write, varid_dummy_put_get_var_float, start, count, NULL, put_var_buffer_double + (PUT_GET_VAR_LEN / 2)); ERR
    /* Put 1st half data */
    start[0] = 0;
    count[0] = PUT_GET_VAR_LEN / 2;
    ret = PIOc_put_vars_double(ncid_write, varid_dummy_put_get_var_float, start, count, NULL, put_var_buffer_double); ERR

    /* Write to int type variable with int type decomposition, type conversions will not be performed. */
    ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_int, ioid_int, ELEMENTS_PER_PE, write_darray_buffer_int, NULL); ERR

    /* Write to float type variable with double type decomposition, type conversions will be performed. */
    ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float, ioid_double, ELEMENTS_PER_PE, write_darray_buffer_double, NULL); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  /* Read support is not implemented for ADIOS type yet: change to "fmt < 2" for testing this new feature later.
     Currently, ADIOS type in SCORPIO simply changes actual type to PnetCDF or NetCDF for reading. */
  for (int fmt = 0; fmt < 1; fmt++) {
    /* Note: for ADIOS type, the actual file name on disk is a BP directory named example1_1.nc.bp.dir,
       client code like E3SM still assumes example1_1.nc as the file name to be read.  */
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    /* Get int type global attribute. */
    ret = PIOc_get_att(ncid_read, PIO_GLOBAL, "dummy_global_att_int", &get_global_att_int_data); ERR
    if (get_global_att_int_data != put_global_att_int_data)
      printf("rank = %d, read wrong data for dummy_global_att_int\n", my_rank);

    /* Get float type global attribute. */
    ret = PIOc_get_att(ncid_read, PIO_GLOBAL, "dummy_global_att_float", &get_global_att_float_data); ERR
    diff_float = get_global_att_float_data - put_global_att_float_data;
    if (fabs(diff_float) > 1E-5)
      printf("rank = %d, read wrong data for dummy_global_att_float\n", my_rank);

    /* Get text type global attribute. */
    ret = PIOc_get_att_text(ncid_read, PIO_GLOBAL, "dummy_global_att_text", get_global_att_text_data); ERR
    if (strncmp(get_global_att_text_data, put_global_att_text_data, strlen(put_global_att_text_data)))
      printf("rank = %d, read wrong data for dummy_global_att_text\n", my_rank);

    /* Get int type scalar variable. */
    ret = PIOc_inq_varid(ncid_read, "dummy_scalar_var_int", &varid_dummy_scalar_var_int); ERR
    ret = PIOc_get_var_int(ncid_read, varid_dummy_scalar_var_int, &get_scalar_int_data); ERR
    if (get_scalar_int_data != put_scalar_int_data)
      printf("rank = %d, read wrong data for dummy_scalar_var_int\n", my_rank);

    /* Get float type scalar variable. */
    ret = PIOc_inq_varid(ncid_read, "dummy_scalar_var_float", &varid_dummy_scalar_var_float); ERR
    ret = PIOc_get_var_float(ncid_read, varid_dummy_scalar_var_float, &get_scalar_float_data); ERR
    diff_float = get_scalar_float_data - put_scalar_float_data;
    if (fabs(diff_float) > 1E-5)
      printf("rank = %d, read wrong data for dummy_scalar_var_float\n", my_rank);

    /* Get text type variable. */
    ret = PIOc_inq_varid(ncid_read, "dummy_text_var", &varid_dummy_text_var); ERR
    ret = PIOc_get_var_text(ncid_read, varid_dummy_text_var, get_text_data); ERR
    if (strncmp(get_text_data, put_text_data, strlen(put_text_data)))
      printf("rank = %d, read wrong data for dummy_text_var\n", my_rank);

    /* Read int type variable with int type decomposition, type conversions will not be performed. */
    ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_int", &varid_dummy_darray_var_int); ERR
    ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_int, ioid_int, ELEMENTS_PER_PE, read_darray_buffer_int); ERR
    for (int i = 0; i < ELEMENTS_PER_PE; i++) {
      if (read_darray_buffer_int[i] != write_darray_buffer_int[i]) {
          printf("rank = %d, read wrong data for dummy_darray_var_int at index %d\n", my_rank, i);
          break;
      }
    }

    /* Get int type attribute of variable dummy_darray_var_int. */
    ret = PIOc_get_att(ncid_read, varid_dummy_darray_var_int, "dummy_att_int", &get_att_int_data); ERR
    if (get_att_int_data != put_att_int_data)
      printf("rank = %d, read wrong data for dummy_att_int of dummy_darray_var_int\n", my_rank);

    /* Get float type attribute of variable dummy_darray_var_int. */
    ret = PIOc_get_att(ncid_read, varid_dummy_darray_var_int, "dummy_att_float", &get_att_float_data); ERR
    diff_float = get_att_float_data - put_att_float_data;
    if (fabs(diff_float) > 1E-5)
      printf("rank = %d, read wrong data for dummy_att_float of dummy_darray_var_int\n", my_rank);

    /* Get text type attribute of variable dummy_darray_var_int. */
    ret = PIOc_get_att_text(ncid_read, varid_dummy_darray_var_int, "dummy_att_text", get_att_text_data); ERR
    if (strncmp(get_att_text_data, put_att_text_data, strlen(put_att_text_data)))
      printf("rank = %d, read wrong data for dummy_att_text of dummy_darray_var_int\n", my_rank);

    /* Read float type variable with double type decomposition, type conversions will be performed. */
    ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float", &varid_dummy_darray_var_float); ERR
    ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float, ioid_double, ELEMENTS_PER_PE, read_darray_buffer_double); ERR
    for (int i = 0; i < ELEMENTS_PER_PE; i++) {
      diff_double = read_darray_buffer_double[i] - write_darray_buffer_double[i];
      if (fabs(diff_double) > 1E-5) {
          printf("rank = %d, read wrong data for dummy_darray_var_double at index %d\n", my_rank, i);
          break;
      }
    }

    /* Get int type variable with int type decomposition, type conversions will not be performed. */
    ret = PIOc_inq_varid(ncid_read, "dummy_put_get_var_int", &varid_dummy_put_get_var_int); ERR
    /* Partial get: excluding the first and the last elements */
    start[0] = 1;
    count[0] = PUT_GET_VAR_LEN - 2;
    ret = PIOc_get_vars_int(ncid_read, varid_dummy_put_get_var_int, start, count, NULL, get_var_buffer_int + 1); ERR
    for (int i = 1; i < PUT_GET_VAR_LEN - 1; i++) {
      if (get_var_buffer_int[i] != put_var_buffer_int[i]) {
          printf("rank = %d, get wrong data for dummy_put_get_var_int at index %d\n", my_rank, i);
          break;
      }
    }

    /* Get float type variable with double type decomposition, type conversions will be performed. */
    ret = PIOc_inq_varid(ncid_read, "dummy_put_get_var_float", &varid_dummy_put_get_var_float); ERR
    /* Partial get: excluding the first and the last elements */
    start[0] = 1;
    count[0] = PUT_GET_VAR_LEN - 2;
    ret = PIOc_get_vars_double(ncid_read, varid_dummy_put_get_var_float, start, count, NULL, get_var_buffer_double + 1); ERR
    for (int i = 1; i < PUT_GET_VAR_LEN - 1; i++) {
      diff_double = get_var_buffer_double[i] - put_var_buffer_double[i];
      if (fabs(diff_double) > 1E-5) {
          printf("rank = %d, get wrong data for dummy_put_get_var_float at index %d\n", my_rank, i);
          break;
      }
    }

    ret = PIOc_closefile(ncid_read); ERR
  }

  ret = PIOc_freedecomp(iosysid, ioid_int); ERR
  ret = PIOc_freedecomp(iosysid, ioid_double); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
