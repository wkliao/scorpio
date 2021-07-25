#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define NUM_NETCDF_FLAVORS 2

int main(int argc, char* argv[])
{
    int my_rank;
    int ntasks;
    int format[NUM_NETCDF_FLAVORS];
    int niotasks = 4;
    int ioproc_stride = 4;
    int ioproc_start = 0;
    int iosysid;
    int ncid;
    int varids[268];
    PIO_Offset start[2];
    PIO_Offset count[2];
    char filename[NC_MAX_NAME + 1];
    int num_flavors = 0;
    int mode = 512;

    int cal_strlen_id;
    int lev_id;
    int ilev_id ;
    int ncol_d_id;
    int timelevels_id;
    int ncol_id;
    int pbuf_00072_id;
    int pbuf_00073_id;
    int pbuf_00288_id;
    int max_chars_id;
    int prescraero_randn_seed_dim_id;
    int pcnst_id;
    int ptapes_id;
    int max_string_len_id;
    int fieldname_lenp2_id;
    int pflds_id;
    int max_fieldname_len_id;
    int maxnflds_id;
    int maxvarmdims_id;
    int registeredmdims_id;
    int max_hcoordname_len_id;

    int dimids[3];

    char dummy_chars[256];
    int dummy_int;
    double dummy_double;
    double val1 = 1e+36;
    double val2 = 1e+30;

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
        sprintf(filename, "test_%d.nc", fmt);

        PIOc_createfile(iosysid, &ncid, &(format[fmt]), filename, mode);

        /* Dimensions */
        PIOc_def_dim(ncid, "cal_strlen", 32, &cal_strlen_id);
        PIOc_def_dim(ncid, "lev", 72, &lev_id);
        PIOc_def_dim(ncid, "ilev", 73, &ilev_id);
        PIOc_def_dim(ncid, "ncol_d", 1536, &ncol_d_id);
        PIOc_def_dim(ncid, "timelevels", PIO_UNLIMITED, &timelevels_id);
        PIOc_def_dim(ncid, "ncol", 866, &ncol_id);
        PIOc_def_dim(ncid, "pbuf_00072", 72, &pbuf_00072_id);
        PIOc_def_dim(ncid, "pbuf_00073", 73, &pbuf_00073_id);
        PIOc_def_dim(ncid, "pbuf_00288", 288, &pbuf_00288_id);
        PIOc_def_dim(ncid, "max_chars", 256, &max_chars_id);
        PIOc_def_dim(ncid, "prescraero_randn_seed_dim", 4, &prescraero_randn_seed_dim_id);
        PIOc_def_dim(ncid, "pcnst", 40, &pcnst_id);
        PIOc_def_dim(ncid, "ptapes", 12, &ptapes_id);
        PIOc_def_dim(ncid, "max_string_len", 256, &max_string_len_id);
        PIOc_def_dim(ncid, "fieldname_lenp2", 26, &fieldname_lenp2_id);
        PIOc_def_dim(ncid, "pflds", 1000, &pflds_id);
        PIOc_def_dim(ncid, "max_fieldname_len", 27, &max_fieldname_len_id);
        PIOc_def_dim(ncid, "maxnflds", 384, &maxnflds_id);
        PIOc_def_dim(ncid, "maxvarmdims", 1, &maxvarmdims_id);
        PIOc_def_dim(ncid, "registeredmdims", 2, &registeredmdims_id);
        PIOc_def_dim(ncid, "max_hcoordname_len", 16, &max_hcoordname_len_id);

        /* 268 variables */
        dimids[0] = cal_strlen_id;
        PIOc_def_var(ncid, "rst_calendar", PIO_CHAR, 1, dimids, &varids[0]);

        PIOc_def_var(ncid, "rst_nstep", PIO_INT, 0, NULL, &varids[1]);
        PIOc_def_var(ncid, "rst_step_days", PIO_INT, 0, NULL, &varids[2]);
        PIOc_def_var(ncid, "rst_step_sec", PIO_INT, 0, NULL, &varids[3]);
        PIOc_def_var(ncid, "rst_start_ymd", PIO_INT, 0, NULL, &varids[4]);
        PIOc_def_var(ncid, "rst_start_tod", PIO_INT, 0, NULL, &varids[5]);
        PIOc_def_var(ncid, "rst_stop_ymd", PIO_INT, 0, NULL, &varids[6]);
        PIOc_def_var(ncid, "rst_stop_tod", PIO_INT, 0, NULL, &varids[7]);
        PIOc_def_var(ncid, "rst_ref_ymd", PIO_INT, 0, NULL, &varids[8]);
        PIOc_def_var(ncid, "rst_ref_tod", PIO_INT, 0, NULL, &varids[9]);
        PIOc_def_var(ncid, "rst_curr_ymd", PIO_INT, 0, NULL, &varids[10]);
        PIOc_def_var(ncid, "rst_curr_tod", PIO_INT, 0, NULL, &varids[11]);
        PIOc_def_var(ncid, "rst_perp_ymd", PIO_INT, 0, NULL, &varids[12]);
        PIOc_def_var(ncid, "rst_perp_cal_int", PIO_INT, 0, NULL, &varids[13]);

        dimids[0] = ilev_id;
        PIOc_def_var(ncid, "hyai", PIO_DOUBLE, 1, dimids, &varids[14]);

        dimids[0] = lev_id;
        PIOc_def_var(ncid, "hyam", PIO_DOUBLE, 1, dimids, &varids[15]);

        dimids[0] = ilev_id;
        PIOc_def_var(ncid, "hybi", PIO_DOUBLE, 1, dimids, &varids[16]);

        dimids[0] = lev_id;
        PIOc_def_var(ncid, "hybm", PIO_DOUBLE, 1, dimids, &varids[17]);

        dimids[0] = timelevels_id;
        PIOc_def_var(ncid, "time", PIO_DOUBLE, 1, dimids, &varids[18]);

        dimids[0] = timelevels_id;
        dimids[1] = lev_id;
        dimids[2] = ncol_d_id;
        PIOc_def_var(ncid, "U", PIO_DOUBLE, 3, dimids, &varids[19]);
        PIOc_def_var(ncid, "V", PIO_DOUBLE, 3, dimids, &varids[20]);
        PIOc_def_var(ncid, "T", PIO_DOUBLE, 3, dimids, &varids[21]);

        dimids[0] = lev_id;
        dimids[1] = ncol_d_id;
        PIOc_def_var(ncid, "OMEGA", PIO_DOUBLE, 2, dimids, &varids[22]);

        dimids[0] = timelevels_id;
        dimids[1] = ncol_d_id;
        PIOc_def_var(ncid, "PS", PIO_DOUBLE, 2, dimids, &varids[23]);

        dimids[0] = ncol_d_id;
        PIOc_def_var(ncid, "PHIS", PIO_DOUBLE, 1, dimids, &varids[24]);

        dimids[0] = timelevels_id;
        dimids[1] = lev_id;
        dimids[2] = ncol_d_id;
        PIOc_def_var(ncid, "dpQ", PIO_DOUBLE, 3, dimids, &varids[25]);
        PIOc_def_var(ncid, "dpCLDLIQ", PIO_DOUBLE, 3, dimids, &varids[26]);
        PIOc_def_var(ncid, "dpCLDICE", PIO_DOUBLE, 3, dimids, &varids[27]);
        PIOc_def_var(ncid, "dpNUMLIQ", PIO_DOUBLE, 3, dimids, &varids[28]);
        PIOc_def_var(ncid, "dpNUMICE", PIO_DOUBLE, 3, dimids, &varids[29]);
        PIOc_def_var(ncid, "dpRAINQM", PIO_DOUBLE, 3, dimids, &varids[30]);
        PIOc_def_var(ncid, "dpSNOWQM", PIO_DOUBLE, 3, dimids, &varids[31]);
        PIOc_def_var(ncid, "dpNUMRAI", PIO_DOUBLE, 3, dimids, &varids[32]);
        PIOc_def_var(ncid, "dpNUMSNO", PIO_DOUBLE, 3, dimids, &varids[33]);
        PIOc_def_var(ncid, "dpO3", PIO_DOUBLE, 3, dimids, &varids[34]);
        PIOc_def_var(ncid, "dpH2O2", PIO_DOUBLE, 3, dimids, &varids[35]);
        PIOc_def_var(ncid, "dpH2SO4", PIO_DOUBLE, 3, dimids, &varids[36]);
        PIOc_def_var(ncid, "dpSO2", PIO_DOUBLE, 3, dimids, &varids[37]);
        PIOc_def_var(ncid, "dpDMS", PIO_DOUBLE, 3, dimids, &varids[38]);
        PIOc_def_var(ncid, "dpSOAG", PIO_DOUBLE, 3, dimids, &varids[39]);
        PIOc_def_var(ncid, "dpso4_a1", PIO_DOUBLE, 3, dimids, &varids[40]);
        PIOc_def_var(ncid, "dppom_a1", PIO_DOUBLE, 3, dimids, &varids[41]);
        PIOc_def_var(ncid, "dpsoa_a1", PIO_DOUBLE, 3, dimids, &varids[42]);
        PIOc_def_var(ncid, "dpbc_a1", PIO_DOUBLE, 3, dimids, &varids[43]);
        PIOc_def_var(ncid, "dpdst_a1", PIO_DOUBLE, 3, dimids, &varids[44]);
        PIOc_def_var(ncid, "dpncl_a1", PIO_DOUBLE, 3, dimids, &varids[45]);
        PIOc_def_var(ncid, "dpmom_a1", PIO_DOUBLE, 3, dimids, &varids[46]);
        PIOc_def_var(ncid, "dpnum_a1", PIO_DOUBLE, 3, dimids, &varids[47]);
        PIOc_def_var(ncid, "dpso4_a2", PIO_DOUBLE, 3, dimids, &varids[48]);
        PIOc_def_var(ncid, "dpsoa_a2", PIO_DOUBLE, 3, dimids, &varids[49]);
        PIOc_def_var(ncid, "dpncl_a2", PIO_DOUBLE, 3, dimids, &varids[50]);
        PIOc_def_var(ncid, "dpmom_a2", PIO_DOUBLE, 3, dimids, &varids[51]);
        PIOc_def_var(ncid, "dpnum_a2", PIO_DOUBLE, 3, dimids, &varids[52]);
        PIOc_def_var(ncid, "dpdst_a3", PIO_DOUBLE, 3, dimids, &varids[53]);
        PIOc_def_var(ncid, "dpncl_a3", PIO_DOUBLE, 3, dimids, &varids[54]);
        PIOc_def_var(ncid, "dpso4_a3", PIO_DOUBLE, 3, dimids, &varids[55]);
        PIOc_def_var(ncid, "dpbc_a3", PIO_DOUBLE, 3, dimids, &varids[56]);
        PIOc_def_var(ncid, "dppom_a3", PIO_DOUBLE, 3, dimids, &varids[57]);
        PIOc_def_var(ncid, "dpsoa_a3", PIO_DOUBLE, 3, dimids, &varids[58]);
        PIOc_def_var(ncid, "dpmom_a3", PIO_DOUBLE, 3, dimids, &varids[59]);
        PIOc_def_var(ncid, "dpnum_a3", PIO_DOUBLE, 3, dimids, &varids[60]);
        PIOc_def_var(ncid, "dppom_a4", PIO_DOUBLE, 3, dimids, &varids[61]);
        PIOc_def_var(ncid, "dpbc_a4", PIO_DOUBLE, 3, dimids, &varids[62]);
        PIOc_def_var(ncid, "dpmom_a4", PIO_DOUBLE, 3, dimids, &varids[63]);
        PIOc_def_var(ncid, "dpnum_a4", PIO_DOUBLE, 3, dimids, &varids[64]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "lat", PIO_DOUBLE, 1, dimids, &varids[65]);
        PIOc_def_var(ncid, "lon", PIO_DOUBLE, 1, dimids, &varids[66]);
        PIOc_def_var(ncid, "area", PIO_DOUBLE, 1, dimids, &varids[67]);
        PIOc_def_var(ncid, "static_ener_ac", PIO_DOUBLE, 1, dimids, &varids[68]);
        PIOc_def_var(ncid, "water_vap_ac", PIO_DOUBLE, 1, dimids, &varids[69]);
        PIOc_def_var(ncid, "TEOUT", PIO_DOUBLE, 1, dimids, &varids[70]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "DTCORE", PIO_DOUBLE, 2, dimids, &varids[71]);
        PIOc_def_var(ncid, "CLDO", PIO_DOUBLE, 2, dimids, &varids[72]);
        PIOc_def_var(ncid, "PRER_EVAP", PIO_DOUBLE, 2, dimids, &varids[73]);
        PIOc_def_var(ncid, "CC_T", PIO_DOUBLE, 2, dimids, &varids[74]);
        PIOc_def_var(ncid, "CC_qv", PIO_DOUBLE, 2, dimids, &varids[75]);
        PIOc_def_var(ncid, "CC_ql", PIO_DOUBLE, 2, dimids, &varids[76]);
        PIOc_def_var(ncid, "CC_qi", PIO_DOUBLE, 2, dimids, &varids[77]);
        PIOc_def_var(ncid, "CC_nl", PIO_DOUBLE, 2, dimids, &varids[78]);
        PIOc_def_var(ncid, "CC_ni", PIO_DOUBLE, 2, dimids, &varids[79]);
        PIOc_def_var(ncid, "CC_qlst", PIO_DOUBLE, 2, dimids, &varids[80]);
        PIOc_def_var(ncid, "am_evp_st", PIO_DOUBLE, 2, dimids, &varids[81]);
        PIOc_def_var(ncid, "evprain_st", PIO_DOUBLE, 2, dimids, &varids[82]);
        PIOc_def_var(ncid, "evpsnow_st", PIO_DOUBLE, 2, dimids, &varids[83]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "ACPRECL", PIO_DOUBLE, 1, dimids, &varids[84]);
        PIOc_def_var(ncid, "ACGCME", PIO_DOUBLE, 1, dimids, &varids[85]);
        PIOc_def_var(ncid, "ACNUM", PIO_INT, 1, dimids, &varids[86]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "RELVAR", PIO_DOUBLE, 2, dimids, &varids[87]);
        PIOc_def_var(ncid, "ACCRE_ENHAN", PIO_DOUBLE, 2, dimids, &varids[88]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "pblh", PIO_DOUBLE, 1, dimids, &varids[89]);

        dimids[0] = pbuf_00073_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "tke", PIO_DOUBLE, 2, dimids, &varids[90]);
        PIOc_def_var(ncid, "kvh", PIO_DOUBLE, 2, dimids, &varids[91]);
        PIOc_def_var(ncid, "kvm", PIO_DOUBLE, 2, dimids, &varids[92]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "tpert", PIO_DOUBLE, 1, dimids, &varids[93]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "AST", PIO_DOUBLE, 2, dimids, &varids[94]);
        PIOc_def_var(ncid, "AIST", PIO_DOUBLE, 2, dimids, &varids[95]);
        PIOc_def_var(ncid, "ALST", PIO_DOUBLE, 2, dimids, &varids[96]);
        PIOc_def_var(ncid, "QIST", PIO_DOUBLE, 2, dimids, &varids[97]);
        PIOc_def_var(ncid, "QLST", PIO_DOUBLE, 2, dimids, &varids[98]);
        PIOc_def_var(ncid, "CONCLD", PIO_DOUBLE, 2, dimids, &varids[99]);
        PIOc_def_var(ncid, "CLD", PIO_DOUBLE, 2, dimids, &varids[100]);
        PIOc_def_var(ncid, "RAD_CLUBB", PIO_DOUBLE, 2, dimids, &varids[101]);

        dimids[0] = pbuf_00073_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "WP2_nadv", PIO_DOUBLE, 2, dimids, &varids[102]);
        PIOc_def_var(ncid, "WP3_nadv", PIO_DOUBLE, 2, dimids, &varids[103]);
        PIOc_def_var(ncid, "WPTHLP_nadv", PIO_DOUBLE, 2, dimids, &varids[104]);
        PIOc_def_var(ncid, "WPRTP_nadv", PIO_DOUBLE, 2, dimids, &varids[105]);
        PIOc_def_var(ncid, "RTPTHLP_nadv", PIO_DOUBLE, 2, dimids, &varids[106]);
        PIOc_def_var(ncid, "RTP2_nadv", PIO_DOUBLE, 2, dimids, &varids[107]);
        PIOc_def_var(ncid, "THLP2_nadv", PIO_DOUBLE, 2, dimids, &varids[108]);
        PIOc_def_var(ncid, "UP2_nadv", PIO_DOUBLE, 2, dimids, &varids[109]);
        PIOc_def_var(ncid, "VP2_nadv", PIO_DOUBLE, 2, dimids, &varids[110]);
        PIOc_def_var(ncid, "UPWP", PIO_DOUBLE, 2, dimids, &varids[111]);
        PIOc_def_var(ncid, "VPWP", PIO_DOUBLE, 2, dimids, &varids[112]);
        PIOc_def_var(ncid, "THLM", PIO_DOUBLE, 2, dimids, &varids[113]);
        PIOc_def_var(ncid, "RTM", PIO_DOUBLE, 2, dimids, &varids[114]);
        PIOc_def_var(ncid, "UM", PIO_DOUBLE, 2, dimids, &varids[115]);
        PIOc_def_var(ncid, "VM", PIO_DOUBLE, 2, dimids, &varids[116]);
        PIOc_def_var(ncid, "WPTHVP", PIO_DOUBLE, 2, dimids, &varids[117]);
        PIOc_def_var(ncid, "WP2THVP", PIO_DOUBLE, 2, dimids, &varids[118]);
        PIOc_def_var(ncid, "RTPTHVP", PIO_DOUBLE, 2, dimids, &varids[119]);
        PIOc_def_var(ncid, "THLPTHVP", PIO_DOUBLE, 2, dimids, &varids[120]);
        PIOc_def_var(ncid, "RCM", PIO_DOUBLE, 2, dimids, &varids[121]);
        PIOc_def_var(ncid, "CLOUD_FRAC", PIO_DOUBLE, 2, dimids, &varids[122]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "vmag_gust", PIO_DOUBLE, 1, dimids, &varids[123]);

        dimids[0] = pbuf_00288_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "DGNUM", PIO_DOUBLE, 2, dimids, &varids[124]);
        PIOc_def_var(ncid, "DGNUMWET", PIO_DOUBLE, 2, dimids, &varids[125]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "num_c1", PIO_DOUBLE, 2, dimids, &varids[126]);
        PIOc_def_var(ncid, "so4_c1", PIO_DOUBLE, 2, dimids, &varids[127]);
        PIOc_def_var(ncid, "pom_c1", PIO_DOUBLE, 2, dimids, &varids[128]);
        PIOc_def_var(ncid, "soa_c1", PIO_DOUBLE, 2, dimids, &varids[129]);
        PIOc_def_var(ncid, "bc_c1", PIO_DOUBLE, 2, dimids, &varids[130]);
        PIOc_def_var(ncid, "dst_c1", PIO_DOUBLE, 2, dimids, &varids[131]);
        PIOc_def_var(ncid, "ncl_c1", PIO_DOUBLE, 2, dimids, &varids[132]);
        PIOc_def_var(ncid, "mom_c1", PIO_DOUBLE, 2, dimids, &varids[133]);
        PIOc_def_var(ncid, "num_c2", PIO_DOUBLE, 2, dimids, &varids[134]);
        PIOc_def_var(ncid, "so4_c2", PIO_DOUBLE, 2, dimids, &varids[135]);
        PIOc_def_var(ncid, "soa_c2", PIO_DOUBLE, 2, dimids, &varids[136]);
        PIOc_def_var(ncid, "ncl_c2", PIO_DOUBLE, 2, dimids, &varids[137]);
        PIOc_def_var(ncid, "mom_c2", PIO_DOUBLE, 2, dimids, &varids[138]);
        PIOc_def_var(ncid, "num_c3", PIO_DOUBLE, 2, dimids, &varids[139]);
        PIOc_def_var(ncid, "dst_c3", PIO_DOUBLE, 2, dimids, &varids[140]);
        PIOc_def_var(ncid, "ncl_c3", PIO_DOUBLE, 2, dimids, &varids[141]);
        PIOc_def_var(ncid, "so4_c3", PIO_DOUBLE, 2, dimids, &varids[142]);
        PIOc_def_var(ncid, "bc_c3", PIO_DOUBLE, 2, dimids, &varids[143]);
        PIOc_def_var(ncid, "pom_c3", PIO_DOUBLE, 2, dimids, &varids[144]);
        PIOc_def_var(ncid, "soa_c3", PIO_DOUBLE, 2, dimids, &varids[145]);
        PIOc_def_var(ncid, "mom_c3", PIO_DOUBLE, 2, dimids, &varids[146]);
        PIOc_def_var(ncid, "num_c4", PIO_DOUBLE, 2, dimids, &varids[147]);
        PIOc_def_var(ncid, "pom_c4", PIO_DOUBLE, 2, dimids, &varids[148]);
        PIOc_def_var(ncid, "bc_c4", PIO_DOUBLE, 2, dimids, &varids[149]);
        PIOc_def_var(ncid, "mom_c4", PIO_DOUBLE, 2, dimids, &varids[150]);

        dimids[0] = pbuf_00073_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "DP_FLXPRC", PIO_DOUBLE, 2, dimids, &varids[151]);
        PIOc_def_var(ncid, "DP_FLXSNW", PIO_DOUBLE, 2, dimids, &varids[152]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "DP_CLDLIQ", PIO_DOUBLE, 2, dimids, &varids[153]);
        PIOc_def_var(ncid, "DP_CLDICE", PIO_DOUBLE, 2, dimids, &varids[154]);
        PIOc_def_var(ncid, "HU_NM1", PIO_DOUBLE, 2, dimids, &varids[155]);
        PIOc_def_var(ncid, "CNV_NM1", PIO_DOUBLE, 2, dimids, &varids[156]);
        PIOc_def_var(ncid, "TM1", PIO_DOUBLE, 2, dimids, &varids[157]);
        PIOc_def_var(ncid, "QM1", PIO_DOUBLE, 2, dimids, &varids[158]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "cush", PIO_DOUBLE, 1, dimids, &varids[159]);

        dimids[0] = pbuf_00072_id;
        dimids[1] = ncol_id;
        PIOc_def_var(ncid, "QRS", PIO_DOUBLE, 2, dimids, &varids[160]);
        PIOc_def_var(ncid, "QRL", PIO_DOUBLE, 2, dimids, &varids[161]);
        PIOc_def_var(ncid, "ICIWP", PIO_DOUBLE, 2, dimids, &varids[162]);
        PIOc_def_var(ncid, "ICLWP", PIO_DOUBLE, 2, dimids, &varids[163]);
        PIOc_def_var(ncid, "T_TTEND", PIO_DOUBLE, 2, dimids, &varids[164]);
        PIOc_def_var(ncid, "FRONTGF", PIO_DOUBLE, 2, dimids, &varids[165]);
        PIOc_def_var(ncid, "FRONTGA", PIO_DOUBLE, 2, dimids, &varids[166]);

        PIOc_def_var(ncid, "pbuf_time_idx", PIO_INT, 0, NULL, &varids[167]);

        dimids[0] = max_chars_id;
        PIOc_def_var(ncid, "tracer_cnst_curr_fname", PIO_CHAR, 1, dimids, &varids[168]);
        PIOc_def_var(ncid, "linoz_data_curr_fname", PIO_CHAR, 1, dimids, &varids[169]);

        dimids[0] = prescraero_randn_seed_dim_id;
        PIOc_def_var(ncid, "prescraero_randn_seed", PIO_INT, 1, dimids, &varids[170]);

        dimids[0] = ncol_id;
        PIOc_def_var(ncid, "LANDM", PIO_DOUBLE, 1, dimids, &varids[171]);
        PIOc_def_var(ncid, "SGH", PIO_DOUBLE, 1, dimids, &varids[172]);
        PIOc_def_var(ncid, "SGH30", PIO_DOUBLE, 1, dimids, &varids[173]);
        PIOc_def_var(ncid, "TREFMXAV", PIO_DOUBLE, 1, dimids, &varids[174]);
        PIOc_def_var(ncid, "TREFMNAV", PIO_DOUBLE, 1, dimids, &varids[175]);
        PIOc_def_var(ncid, "FLWDS", PIO_DOUBLE, 1, dimids, &varids[176]);
        PIOc_def_var(ncid, "SOLS", PIO_DOUBLE, 1, dimids, &varids[177]);
        PIOc_def_var(ncid, "SOLL", PIO_DOUBLE, 1, dimids, &varids[178]);
        PIOc_def_var(ncid, "SOLSD", PIO_DOUBLE, 1, dimids, &varids[179]);
        PIOc_def_var(ncid, "SOLLD", PIO_DOUBLE, 1, dimids, &varids[180]);
        PIOc_def_var(ncid, "BCPHIDRY", PIO_DOUBLE, 1, dimids, &varids[181]);
        PIOc_def_var(ncid, "BCPHODRY", PIO_DOUBLE, 1, dimids, &varids[182]);
        PIOc_def_var(ncid, "OCPHIDRY", PIO_DOUBLE, 1, dimids, &varids[183]);
        PIOc_def_var(ncid, "OCPHODRY", PIO_DOUBLE, 1, dimids, &varids[184]);
        PIOc_def_var(ncid, "DSTDRY1", PIO_DOUBLE, 1, dimids, &varids[185]);
        PIOc_def_var(ncid, "DSTDRY2", PIO_DOUBLE, 1, dimids, &varids[186]);
        PIOc_def_var(ncid, "DSTDRY3", PIO_DOUBLE, 1, dimids, &varids[187]);
        PIOc_def_var(ncid, "DSTDRY4", PIO_DOUBLE, 1, dimids, &varids[188]);
        PIOc_def_var(ncid, "CFLX0001", PIO_DOUBLE, 1, dimids, &varids[189]);
        PIOc_def_var(ncid, "CFLX0002", PIO_DOUBLE, 1, dimids, &varids[190]);
        PIOc_def_var(ncid, "CFLX0003", PIO_DOUBLE, 1, dimids, &varids[191]);
        PIOc_def_var(ncid, "CFLX0004", PIO_DOUBLE, 1, dimids, &varids[192]);
        PIOc_def_var(ncid, "CFLX0005", PIO_DOUBLE, 1, dimids, &varids[193]);
        PIOc_def_var(ncid, "CFLX0006", PIO_DOUBLE, 1, dimids, &varids[194]);
        PIOc_def_var(ncid, "CFLX0007", PIO_DOUBLE, 1, dimids, &varids[195]);
        PIOc_def_var(ncid, "CFLX0008", PIO_DOUBLE, 1, dimids, &varids[196]);
        PIOc_def_var(ncid, "CFLX0009", PIO_DOUBLE, 1, dimids, &varids[197]);
        PIOc_def_var(ncid, "CFLX0010", PIO_DOUBLE, 1, dimids, &varids[198]);
        PIOc_def_var(ncid, "CFLX0011", PIO_DOUBLE, 1, dimids, &varids[199]);
        PIOc_def_var(ncid, "CFLX0012", PIO_DOUBLE, 1, dimids, &varids[200]);
        PIOc_def_var(ncid, "CFLX0013", PIO_DOUBLE, 1, dimids, &varids[201]);
        PIOc_def_var(ncid, "CFLX0014", PIO_DOUBLE, 1, dimids, &varids[202]);
        PIOc_def_var(ncid, "CFLX0015", PIO_DOUBLE, 1, dimids, &varids[203]);
        PIOc_def_var(ncid, "CFLX0016", PIO_DOUBLE, 1, dimids, &varids[204]);
        PIOc_def_var(ncid, "CFLX0017", PIO_DOUBLE, 1, dimids, &varids[205]);
        PIOc_def_var(ncid, "CFLX0018", PIO_DOUBLE, 1, dimids, &varids[206]);
        PIOc_def_var(ncid, "CFLX0019", PIO_DOUBLE, 1, dimids, &varids[207]);
        PIOc_def_var(ncid, "CFLX0020", PIO_DOUBLE, 1, dimids, &varids[208]);
        PIOc_def_var(ncid, "CFLX0021", PIO_DOUBLE, 1, dimids, &varids[209]);
        PIOc_def_var(ncid, "CFLX0022", PIO_DOUBLE, 1, dimids, &varids[210]);
        PIOc_def_var(ncid, "CFLX0023", PIO_DOUBLE, 1, dimids, &varids[211]);
        PIOc_def_var(ncid, "CFLX0024", PIO_DOUBLE, 1, dimids, &varids[212]);
        PIOc_def_var(ncid, "CFLX0025", PIO_DOUBLE, 1, dimids, &varids[213]);
        PIOc_def_var(ncid, "CFLX0026", PIO_DOUBLE, 1, dimids, &varids[214]);
        PIOc_def_var(ncid, "CFLX0027", PIO_DOUBLE, 1, dimids, &varids[215]);
        PIOc_def_var(ncid, "CFLX0028", PIO_DOUBLE, 1, dimids, &varids[216]);
        PIOc_def_var(ncid, "CFLX0029", PIO_DOUBLE, 1, dimids, &varids[217]);
        PIOc_def_var(ncid, "CFLX0030", PIO_DOUBLE, 1, dimids, &varids[218]);
        PIOc_def_var(ncid, "CFLX0031", PIO_DOUBLE, 1, dimids, &varids[219]);
        PIOc_def_var(ncid, "CFLX0032", PIO_DOUBLE, 1, dimids, &varids[220]);
        PIOc_def_var(ncid, "CFLX0033", PIO_DOUBLE, 1, dimids, &varids[221]);
        PIOc_def_var(ncid, "CFLX0034", PIO_DOUBLE, 1, dimids, &varids[222]);
        PIOc_def_var(ncid, "CFLX0035", PIO_DOUBLE, 1, dimids, &varids[223]);
        PIOc_def_var(ncid, "CFLX0036", PIO_DOUBLE, 1, dimids, &varids[224]);
        PIOc_def_var(ncid, "CFLX0037", PIO_DOUBLE, 1, dimids, &varids[225]);
        PIOc_def_var(ncid, "CFLX0038", PIO_DOUBLE, 1, dimids, &varids[226]);
        PIOc_def_var(ncid, "CFLX0039", PIO_DOUBLE, 1, dimids, &varids[227]);
        PIOc_def_var(ncid, "CFLX0040", PIO_DOUBLE, 1, dimids, &varids[228]);
        PIOc_def_var(ncid, "LHF", PIO_DOUBLE, 1, dimids, &varids[229]);
        PIOc_def_var(ncid, "SHF", PIO_DOUBLE, 1, dimids, &varids[230]);

        dimids[0] = ptapes_id;
        PIOc_def_var(ncid, "rgnht", PIO_INT, 1, dimids, &varids[231]);
        PIOc_def_var(ncid, "nhtfrq", PIO_INT, 1, dimids, &varids[232]);
        PIOc_def_var(ncid, "nflds", PIO_INT, 1, dimids, &varids[233]);
        PIOc_def_var(ncid, "nfils", PIO_INT, 1, dimids, &varids[234]);
        PIOc_def_var(ncid, "mfilt", PIO_INT, 1, dimids, &varids[235]);

        dimids[0] = ptapes_id;
        dimids[1] = max_string_len_id;
        PIOc_def_var(ncid, "nfpath", PIO_CHAR, 2, dimids, &varids[236]);
        PIOc_def_var(ncid, "cpath", PIO_CHAR, 2, dimids, &varids[237]);
        PIOc_def_var(ncid, "nhfil", PIO_CHAR, 2, dimids, &varids[238]);

        dimids[0] = ptapes_id;
        PIOc_def_var(ncid, "ndens", PIO_INT, 1, dimids, &varids[239]);

        dimids[0] = ptapes_id;
        dimids[1] = pflds_id;
        dimids[2] = max_chars_id;
        PIOc_def_var(ncid, "fincllonlat", PIO_CHAR, 3, dimids, &varids[240]);

        dimids[0] = ptapes_id;
        PIOc_def_var(ncid, "ncprec", PIO_INT, 1, dimids, &varids[241]);
        PIOc_def_var(ncid, "beg_time", PIO_DOUBLE, 1, dimids, &varids[242]);

        dimids[0] = ptapes_id;
        dimids[1] = pflds_id;
        dimids[2] = fieldname_lenp2_id;
        PIOc_def_var(ncid, "fincl", PIO_CHAR, 3, dimids, &varids[243]);
        PIOc_def_var(ncid, "fexcl", PIO_CHAR, 3, dimids, &varids[244]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        dimids[2] = max_fieldname_len_id;
        PIOc_def_var(ncid, "field_name", PIO_CHAR, 3, dimids, &varids[245]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "decomp_type", PIO_INT, 2, dimids, &varids[246]);
        PIOc_def_var(ncid, "numlev", PIO_INT, 2, dimids, &varids[247]);

        dimids[0] = ptapes_id;
        dimids[1] = max_string_len_id;
        PIOc_def_var(ncid, "hrestpath", PIO_CHAR, 2, dimids, &varids[248]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "hwrt_prec", PIO_INT, 2, dimids, &varids[249]);
        PIOc_def_var(ncid, "avgflag", PIO_CHAR, 2, dimids, &varids[250]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        dimids[2] = max_chars_id;
        PIOc_def_var(ncid, "sampling_seq", PIO_CHAR, 3, dimids, &varids[251]);
        PIOc_def_var(ncid, "long_name", PIO_CHAR, 3, dimids, &varids[252]);
        PIOc_def_var(ncid, "units", PIO_CHAR, 3, dimids, &varids[253]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "xyfill", PIO_INT, 2, dimids, &varids[254]);

        dimids[0] = ptapes_id;
        PIOc_def_var(ncid, "lcltod_start", PIO_INT, 1, dimids, &varids[255]);
        PIOc_def_var(ncid, "lcltod_stop", PIO_INT, 1, dimids, &varids[256]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "fillvalue", PIO_DOUBLE, 2, dimids, &varids[257]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        dimids[2] = maxvarmdims_id;
        PIOc_def_var(ncid, "mdims", PIO_INT, 3, dimids, &varids[258]);

        dimids[0] = registeredmdims_id;
        dimids[1] = max_hcoordname_len_id;
        PIOc_def_var(ncid, "mdimnames", PIO_CHAR, 2, dimids, &varids[259]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "is_subcol", PIO_INT, 2, dimids, &varids[260]);

        dimids[0] = ptapes_id;
        PIOc_def_var(ncid, "interpolate_output", PIO_INT, 1, dimids, &varids[261]);
        PIOc_def_var(ncid, "interpolate_type", PIO_INT, 1, dimids, &varids[262]);
        PIOc_def_var(ncid, "interpolate_gridtype", PIO_INT, 1, dimids, &varids[263]);
        PIOc_def_var(ncid, "interpolate_nlat", PIO_INT, 1, dimids, &varids[264]);
        PIOc_def_var(ncid, "interpolate_nlon", PIO_INT, 1, dimids, &varids[265]);

        dimids[0] = ptapes_id;
        dimids[1] = maxnflds_id;
        PIOc_def_var(ncid, "meridional_complement", PIO_INT, 2, dimids, &varids[266]);
        PIOc_def_var(ncid, "zonal_complement", PIO_INT, 2, dimids, &varids[267]);

        PIOc_enddef(ncid);

        strcpy(dummy_chars, "NO_LEAP");
        start[0] = 0;
        count[0] = strlen(dummy_chars);
        PIOc_put_vars_text(ncid, varids[0], start, count, NULL, dummy_chars);

        dummy_int = -999999999;
        PIOc_put_var_int(ncid, varids[1], &dummy_int);

        dummy_int = -999999999;
        PIOc_put_var_int(ncid, varids[2], &dummy_int);

        dummy_int = 7200;
        PIOc_put_var_int(ncid, varids[3], &dummy_int);

        dummy_int = 10101;
        PIOc_put_var_int(ncid, varids[4], &dummy_int);

        dummy_int = 0;
        PIOc_put_var_int(ncid, varids[5], &dummy_int);

        dummy_int = 10102;
        PIOc_put_var_int(ncid, varids[6], &dummy_int);

        dummy_int = 0;
        PIOc_put_var_int(ncid, varids[7], &dummy_int);

        dummy_int = 10101;
        PIOc_put_var_int(ncid, varids[8], &dummy_int);

        dummy_int = 0;
        PIOc_put_var_int(ncid, varids[9], &dummy_int);

        dummy_int = 10102;
        PIOc_put_var_int(ncid, varids[10], &dummy_int);

        dummy_int = 0;
        PIOc_put_var_int(ncid, varids[11], &dummy_int);

        dummy_int = -999999999;
        PIOc_put_var_int(ncid, varids[12], &dummy_int);

        dummy_int = 0;
        PIOc_put_var_int(ncid, varids[13], &dummy_int);

        start[0] = 0;
        count[0] = 1;
        dummy_double = 1.0;
        PIOc_put_vars_double(ncid, varids[18], start, count, NULL, &dummy_double);

        count[0] = 1;
        count[1] = 1;
        start[0] = 0;
        for (int i = 0; i < 1; i++)
        {
            start[1] = i;
            if (i != 174)
                PIOc_put_vars_double(ncid, varids[257], start, count, NULL, &val1);
            else
                PIOc_put_vars_double(ncid, varids[257], start, count, NULL, &val2);
        }

        start[0] = 1;
        for (int i = 0; i < 21; i++)
        {
            start[1] = i;
            PIOc_put_vars_double(ncid, varids[257], start, count, NULL, &val1);
        }

        start[0] = 11;
        for (int i = 0; i < 52; i++)
        {
            start[1] = i;
            PIOc_put_vars_double(ncid, varids[257], start, count, NULL, &val1);
        }

        PIOc_closefile(ncid);
    }

    PIOc_finalize(iosysid);

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}

