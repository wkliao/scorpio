#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

int main(int argc, char* argv[])
{
  int world_rank;
  int ntasks;
  MPI_Comm comp_comm_all = MPI_COMM_NULL;
  MPI_Comm comp_comm_1st_half = MPI_COMM_NULL;
  MPI_Group group = MPI_GROUP_NULL;
  MPI_Group compgroup = MPI_GROUP_NULL;
  int *compranks_1st_half = NULL;
  int num_comptasks_1st_half = 8;
  int format = PIO_IOTYPE_ADIOS;
  int iosysid1, iosysid2;
  int ncid;
  int i;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  if (ntasks != 16) {
    if (world_rank == 0)
      fprintf(stderr, "Number of processors must be 16!\n");

    MPI_Finalize();
    return 1;
  }

  compranks_1st_half = malloc(num_comptasks_1st_half * sizeof(int));
  for (i = 0; i < num_comptasks_1st_half; i++)
    compranks_1st_half[i] = i;

  /* Create a group for all MPI tasks */
  MPI_Comm_group(MPI_COMM_WORLD, &group);

  /* Create an MPI communicator for all comp tasks: world rank 0 to 15 */
  MPI_Comm_create(MPI_COMM_WORLD, group, &comp_comm_all);

  /* Create a sub-group for 1st half comp tasks: world rank 0 to 7 */
  MPI_Group_incl(group, num_comptasks_1st_half, compranks_1st_half, &compgroup);

  /* Create an MPI communicator for 1st half comp tasks: world rank 0 to 7 */
  MPI_Comm_create(MPI_COMM_WORLD, compgroup, &comp_comm_1st_half);

printf("I am here...\n");

  /* comptasks = 16, iotasks = 4, stride = 4, base = 1 */
  PIOc_Init_Intracomm(comp_comm_all, 4, 4, 1, PIO_REARR_BOX, &iosysid1);


printf("I am here now...\n");

  if (world_rank < 8) {
    /* comptasks = 8, iotasks = 4, stride = 2, base = 1 */
    PIOc_Init_Intracomm(comp_comm_1st_half, 4, 2, 1, PIO_REARR_BOX, &iosysid2);

    /* adios2_open() hangs in this call */
    PIOc_createfile(iosysid2, &ncid, &format, "test_adios.nc", PIO_CLOBBER);

    PIOc_closefile(ncid);
  }


printf("I am here out...\n");

  PIOc_finalize(iosysid1);


printf("I am here out2...\n");

  if (world_rank < 8)
    PIOc_finalize(iosysid2);

printf("I am here out3...\n");

  if (group != MPI_GROUP_NULL)
    MPI_Group_free(&group);

  if (compgroup != MPI_GROUP_NULL)
    MPI_Group_free(&compgroup);

  if (comp_comm_all != MPI_COMM_NULL)
    MPI_Comm_free(&comp_comm_all);

  if (comp_comm_1st_half != MPI_COMM_NULL)
    MPI_Comm_free(&comp_comm_1st_half);

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
