#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include <regex>
#include <unistd.h> // usleep
#include <mpi.h>
#include <sys/types.h>
#include <dirent.h>

#include <adios2.h>

extern "C" {
#include "pio.h"
#include "pio_internal.h"
}

#include "adios2pio-nm-lib.h"
#include "adios2pio-nm-lib-c.h"

using namespace std;

using AttributeVector = std::vector<std::vector<char> >;
using IOVector = std::vector<adios2::IO>;
using EngineVector = std::vector<adios2::Engine>;

#define _MERGE_SCB

/* Debug output */
static int debug_out = 0;

#define BP2PIO_NOERR PIO_NOERR
#define BP2PIO_ERROR -600

#define ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm) \
    err_val = 0; \
    if (ierr != BP2PIO_NOERR) \
    { \
        err_val = 1; \
    } \
    MPI_Allreduce(&err_val, &err_cnt, 1, MPI_INT, MPI_SUM, comm); \
    if (err_cnt != 0) \
    { \
        return BP2PIO_ERROR; \
    }

#define DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm) \
    err_val = 0; \
    if (ierr != BP2PIO_NOERR) \
    { \
        err_val = 1; \
    } \
    MPI_Allreduce(&err_val, &err_cnt, 1, MPI_INT, MPI_SUM, comm); \
    if (err_cnt != 0) \
    { \
        return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR}; \
    }

#define ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, error_msg) \
    err_val = 0; \
    if (ierr != BP2PIO_NOERR) \
    { \
        err_val = 1; \
    } \
    MPI_Allreduce(&err_val, &err_cnt, 1, MPI_INT, MPI_SUM, comm); \
    if (err_cnt != 0) \
    { \
        throw std::runtime_error(error_msg); \
    }

#define ERROR_CHECK_SINGLE_THROW(ierr, error_msg) \
    if (ierr != BP2PIO_NOERR) \
    { \
        throw std::runtime_error(error_msg); \
    }

nc_type PIOc_get_nctype_from_adios_type(const std::string &atype)
{
#define adios2_GET_TYPE(a_type, T, n_type) \
    if (a_type == adios2::GetType<T>()) \
    { \
        return n_type; \
    }

    adios2_GET_TYPE(atype, int8_t, PIO_BYTE);
    adios2_GET_TYPE(atype, int16_t, PIO_SHORT);
    adios2_GET_TYPE(atype, int32_t, PIO_INT);
    adios2_GET_TYPE(atype, int64_t, PIO_INT64);
    adios2_GET_TYPE(atype, uint8_t, PIO_UBYTE);
    adios2_GET_TYPE(atype, uint16_t, PIO_USHORT);
    adios2_GET_TYPE(atype, uint32_t, PIO_UINT);
    adios2_GET_TYPE(atype, uint64_t, PIO_UINT64);

    adios2_GET_TYPE(atype, char, PIO_BYTE);
    adios2_GET_TYPE(atype, short, PIO_SHORT);
    adios2_GET_TYPE(atype, int, PIO_INT);
    adios2_GET_TYPE(atype, float, PIO_FLOAT);
    adios2_GET_TYPE(atype, double, PIO_DOUBLE);
    adios2_GET_TYPE(atype, unsigned char, PIO_UBYTE);
    adios2_GET_TYPE(atype, unsigned short, PIO_USHORT);
    adios2_GET_TYPE(atype, unsigned int, PIO_UINT);
    adios2_GET_TYPE(atype, long, PIO_INT64);
    adios2_GET_TYPE(atype, long long, PIO_INT64);
    adios2_GET_TYPE(atype, unsigned long, PIO_UINT64);
    adios2_GET_TYPE(atype, unsigned long long, PIO_UINT64);
    adios2_GET_TYPE(atype, std::string, PIO_CHAR);

    return PIO_BYTE;

#undef adios2_GET_TYPE
}

int adios2_type_size_a2(const std::string &atype)
{
#define adios2_GET_SIZE(a_type, T, a_size) \
    if (a_type == adios2::GetType<T>()) \
    { \
        return a_size; \
    }

    adios2_GET_SIZE(atype, int8_t, sizeof(int8_t));
    adios2_GET_SIZE(atype, int16_t, sizeof(int16_t));
    adios2_GET_SIZE(atype, int32_t, sizeof(int32_t));
    adios2_GET_SIZE(atype, int64_t, sizeof(int64_t));
    adios2_GET_SIZE(atype, uint8_t, sizeof(uint8_t));
    adios2_GET_SIZE(atype, uint16_t, sizeof(uint16_t));
    adios2_GET_SIZE(atype, uint32_t, sizeof(uint32_t));
    adios2_GET_SIZE(atype, uint64_t, sizeof(uint64_t));

    adios2_GET_SIZE(atype, char, sizeof(char));
    adios2_GET_SIZE(atype, unsigned char, sizeof(unsigned char));
    adios2_GET_SIZE(atype, std::string, 1);
    adios2_GET_SIZE(atype, short, sizeof(short));
    adios2_GET_SIZE(atype, unsigned short, sizeof(unsigned short));
    adios2_GET_SIZE(atype, int, sizeof(int));
    adios2_GET_SIZE(atype, unsigned int, sizeof(unsigned int));
    adios2_GET_SIZE(atype, long, sizeof(long));
    adios2_GET_SIZE(atype, long long, sizeof(long long));
    adios2_GET_SIZE(atype, unsigned long, sizeof(unsigned long));
    adios2_GET_SIZE(atype, unsigned long long, sizeof(unsigned long long));
    adios2_GET_SIZE(atype, float, sizeof(float));
    adios2_GET_SIZE(atype, double, sizeof(double));
    adios2_GET_SIZE(atype, std::complex<float>, sizeof(std::complex<float>));
    adios2_GET_SIZE(atype, std::complex<double>, sizeof(std::complex<double>));

    return -1;

#undef adios2_GET_SIZE
}

template <class T>
int adios2_adios_get_attr_a2(adios2::Attribute<T> &a_base, std::string &atype, AttributeVector &adata)
{
    atype = a_base.Type();
    const std::vector<T> a_data = a_base.Data();
    adata.resize(1);
    adata[0].resize(a_data.size() * sizeof(T));
    memcpy(adata[0].data(), a_data.data(), a_data.size() * sizeof(T));

    return BP2PIO_NOERR;
}

int adios_get_attr_a2(adios2::IO &bpIO, char *aname, std::string &atype, AttributeVector &adata)
{
    int ierr = BP2PIO_NOERR;

    std::string a_type = bpIO.AttributeType(aname);
    if (a_type.empty())
    {
        return BP2PIO_ERROR;
    }
    else if (a_type == adios2::GetType<std::string>())
    {
        try
        {
            adios2::Attribute<std::string> a_base = bpIO.InquireAttribute<std::string>(aname);
			atype = a_base.Type();
            const std::vector<std::string> a_data = a_base.Data();
            adata.resize(a_data.size());
            for (size_t ii = 0; ii < a_data.size(); ii++)
            {
                adata[ii].resize(a_data[ii].length() + 1);
                memcpy(adata[ii].data(), a_data[ii].c_str(), a_data[ii].length() + 1);
            }
        }
        catch (const std::exception &e)
        {
			std::cout << e.what() << std::endl;
            return BP2PIO_ERROR;
        }
        catch (...)
        {
            return BP2PIO_ERROR;
        }

        return BP2PIO_NOERR;
    }

#define my_declare_template_instantiation(T) \
    else if (a_type == adios2::GetType<T>()) \
    { \
        adios2::Attribute<T> a_base = bpIO.InquireAttribute<T>(aname); \
        ierr = adios2_adios_get_attr_a2(a_base, atype, adata); \
        return ierr; \
    }

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(my_declare_template_instantiation)

#undef my_declare_template_instantiation

    return BP2PIO_NOERR;
}

/* Timer functions to separate time for read (ADIOS) and write (PIO) */
#ifdef ADIOS_TIMING
static double time_read, time_write;
static double time_temp_read, time_temp_write;

void TimerInitialize_nm()
{
    time_read = 0.0;
    time_write = 0.0;
}

#define TimerStart(x) { time_temp_##x = MPI_Wtime(); }
#define TimerStop(x) { time_##x += (MPI_Wtime() - time_temp_##x); }

void TimerReport_nm(MPI_Comm comm)
{
    int nproc, rank;
    double tr_sum, tr_max;
    double tw_sum, tw_max;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    MPI_Reduce(&time_read, &tr_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&time_read, &tr_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&time_write, &tw_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&time_write, &tw_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    if (!rank && debug_out)
    {
        cout << "Timing information:     Max     Sum of all\n";
        cout.precision(2);
        cout << "ADIOS read time   = " << std::fixed << std::setw(8) << tr_max << "s "
             << std::setw(8) << tr_sum << "s\n";
        cout << "PIO  write time   = " << std::fixed << std::setw(8) << tw_max << "s "
             << std::setw(8) << tw_sum << "s\n";
    }
}

void TimerFinalize_nm() {}
#else
void TimerInitialize_nm() {}

#define TimerStart(x) {}
#define TimerStop(x) {}

void TimerReport_nm(MPI_Comm comm) {}

void TimerFinalize_nm() {}
#endif

void SetDebugOutput(int val) { debug_out = val; }

struct Dimension
{
    int dimid;
    PIO_Offset dimvalue;
};

using DimensionMap = std::map<std::string, Dimension>;

struct Variable
{
    int         nc_varid;
    bool        is_timed;
    nc_type     nctype;
	int 		ndims;
	std::string op;
	std::string decomp_name;  /* decomposition */
	int         start_time_step; /* a timed variable may be spread across multiple adios time steps */
};

using VariableMap = std::map<std::string, Variable>;

struct Decomposition
{
    int ioid;
    int piotype;
};

using DecompositionMap = std::map<std::string, Decomposition>;
using DecompositionStepMap = std::map<std::string, uint64_t>;
#define NO_DECOMP "no_decomp"

int InitPIO(MPI_Comm comm, int mpirank, int nproc)
{
    int ret = PIO_NOERR;
    int iosysid;

    ret = PIOc_Init_Intracomm(comm, nproc, 1, 0, PIO_REARR_SUBSET, &iosysid);
    if (ret != PIO_NOERR)
        return BP2PIO_ERROR;

    return iosysid;
}

void FlushStdout_nm(MPI_Comm comm)
{
    cout << std::flush;
    usleep((useconds_t)100);
    MPI_Barrier(comm);
}

/* Set the currently encountered max number of steps if argument is given.
 * Return the max step currently
 */
int GlobalMaxSteps_nm(int nsteps_in = 0)
{
    static int nsteps_current = 1;
    if (nsteps_in > nsteps_current)
        nsteps_current = nsteps_in;

    return nsteps_current;
}

std::vector<int> AssignWriteRanks(int n_bp_writers, MPI_Comm comm, int mpirank, int nproc)
{
    if (!mpirank && debug_out)
        cout << "The BP file was written by " << n_bp_writers << " processes\n";

    int nwb = n_bp_writers / nproc; // Number of blocks to process
    int start_wb; // Starting wb (in [0 .. nwb-1])

    if (mpirank < n_bp_writers % nproc)
        nwb++;

    if (mpirank < n_bp_writers % nproc)
        start_wb = mpirank * nwb;
    else
        start_wb = mpirank * nwb + n_bp_writers % nproc;

    if (debug_out)
        cout << "Process " << mpirank << " start block = " << start_wb <<
                " number of blocks = " << nwb << endl;

    std::vector<int> blocks(nwb);
    for (int i = 0; i < nwb; ++i)
        blocks[i] = start_wb + i;

    return blocks;
}

int ProcessGlobalFillmode(adios2::IO &bpIO, int ncid, MPI_Comm comm)
{
    std::string atype;
    AttributeVector fillmode;
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;

    ierr = adios_get_attr_a2(bpIO, (char*)"/__pio__/fillmode", atype, fillmode);
    ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

    ret = PIOc_set_fill(ncid, *((int*)fillmode[0].data()), NULL);
    if (ret != PIO_NOERR)
        return BP2PIO_ERROR;

    return BP2PIO_NOERR;
}

int ProcessVarAttributes(adios2::IO &bpIO, 
						adios2::Engine &bpReader, 
						const std::string &varname,
                        int ncid, 
						int nc_varid, 
						MPI_Comm comm)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;

    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableAttributes(varname);
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        if (debug_out)
            cout << "    Attribute: " << a2_iter->first << endl;

        if (a2_iter->first.find("__pio__/") == string::npos)
        {
            std::string aa_name = varname + "/" + a2_iter->first;
            std::string atype = bpIO.AttributeType(aa_name.c_str());
            nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
            char *attname = (char*)a2_iter->first.c_str();

            if (debug_out)
                cout << "        define PIO attribute: " << attname << ""
                     << "  type=" << piotype << endl;

            AttributeVector adata;
            ierr = adios_get_attr_a2(bpIO, (char*)aa_name.c_str(), atype, adata);
            ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

            if (atype == adios2::GetType<std::string>())
            {
                ret = PIOc_put_att(ncid, nc_varid, attname, piotype, adata[0].size() - 1, adata[0].data());
            }
            else
            {
                ret = PIOc_put_att(ncid, nc_varid, attname, piotype, 1, adata[0].data());
            }
            if (ret != PIO_NOERR)
                return BP2PIO_ERROR;
        }
    }

    return BP2PIO_NOERR;
}

int ProcessGlobalAttributes(adios2::IO &bpIO, adios2::Engine &bpReader, 
							int ncid, 
							DimensionMap& dimension_map,
                            VariableMap& vars_map, 
							MPI_Comm comm,
    						std::map<std::string, char> &processed_attrs,
    						std::map<std::string, int> &var_att_map)
{
    if (debug_out)
        cout << "Process Global Attributes: \n";

    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret  = PIO_NOERR;
    std::string delimiter = "/";

    std::map<std::string, adios2::Params> a2_attr = bpIO.AvailableAttributes();
    size_t nattrs = a2_attr.size();

    std::map<std::string, adios2::Params>::iterator a2_iter = a2_attr.begin();

    int total_cnt = (int) nattrs;
    while (total_cnt > 0)
    {
        std::string a = a2_iter->first;
        if (++a2_iter == a2_attr.end())
            a2_iter = a2_attr.begin();
        char *attr_namelist = (char*)a.c_str();
        std::string token = a.substr(0, a.find(delimiter));

        if (token == "") /* Not a variable with attributes */
        {
            processed_attrs[a] = 1;
            total_cnt--;
        }
        else if (processed_attrs.find(a) == processed_attrs.end())
        {
			/* check if it is a global attribute */
            if (a.find("pio_global/") != string::npos)
            {
                if (debug_out)
                    cout << " GLOBAL Attribute: " << attr_namelist << endl;

                std::string atype = bpIO.AttributeType(attr_namelist);
                nc_type piotype   = PIOc_get_nctype_from_adios_type(atype);
                char *attname     = attr_namelist + strlen("pio_global/");

                if (debug_out)
                    cout << "        define PIO attribute: " << attname << ""
                         << "  type=" << piotype << endl;

                AttributeVector adata;
                ierr = adios_get_attr_a2(bpIO, attr_namelist, atype, adata);
                ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

                if (atype == adios2::GetType<std::string>())
                {
                    ret = PIOc_put_att(ncid, PIO_GLOBAL, attname, piotype, adata[0].size() - 1, adata[0].data());
                }
                else
                {
                    ret = PIOc_put_att(ncid, PIO_GLOBAL, attname, piotype, 1, adata[0].data());
                }
                if (ret != PIO_NOERR)
                    return BP2PIO_ERROR;

                processed_attrs[a] = 1;
                total_cnt--;
            }
            else
            {
                if (debug_out)
                    cout << "    Attribute: " << attr_namelist << endl;

                if (vars_map.find(token) == vars_map.end())
                {
                    if (var_att_map.find(token) == var_att_map.end())
                    {
                        // First encounter
                        std::string atype;
                        AttributeVector adata;
                        string attname = token + "/__pio__/nctype";
                        processed_attrs[attname] = 1;
                        total_cnt--;
                        ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                        ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

                        int nctype = *((int*)adata[0].data());

                        attname = token + "/__pio__/ndims";
                        processed_attrs[attname] = 1;
                        total_cnt--;
                        ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                        ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

                        int ndims = *((int*)adata[0].data());

                        int dimids[PIO_MAX_DIMS];
                        if (ndims)
                        {
                            attname = token + "/__pio__/dims";
                            processed_attrs[attname] = 1;
                            total_cnt--;
                            ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                            ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

                            for (int d = 0; d < ndims; d++)
                                dimids[d] = dimension_map[adata[d].data()].dimid;
                        }
                        else
                        {
                            ierr = BP2PIO_NOERR;
                            ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)
                        }

                        int varid;
                        ret = PIOc_def_var(ncid, token.c_str(), nctype, ndims, dimids, &varid);
                        if (ret != PIO_NOERR)
                            return BP2PIO_ERROR;

                        var_att_map[token] = varid;
                    }
                    else
                    {
                        if (processed_attrs.find(a) == processed_attrs.end())
                        {
                            processed_attrs[a] = 1;
                            total_cnt--;
                            std::string atype = bpIO.AttributeType(a.c_str());
                            nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
                            char *attname = ((char*)a.c_str()) + token.length() + 1;
                            AttributeVector adata;
                            ierr = adios_get_attr_a2(bpIO, (char*)a.c_str(), atype, adata);
                            ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

                            if (atype == adios2::GetType<std::string>())
                            {
                                ret = PIOc_put_att(ncid, var_att_map[token], attname, piotype, adata[0].size() - 1, adata[0].data());
                            }
                            else
                            {
                                ret = PIOc_put_att(ncid, var_att_map[token], attname, piotype, 1, adata[0].data());
                            }
                            if (ret != PIO_NOERR)
                                return BP2PIO_ERROR;
                        }
                        else
                        {
                            ierr = BP2PIO_NOERR;
                            ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)
                        }
                    }
                }
            }
        }
    }

    return BP2PIO_NOERR;
}

template <class T>
Decomposition adios2_ProcessOneDecomposition(adios2::Variable<T> *v_base,
									 	adios2::IO &bpIO, adios2::Engine &bpReader, 
										int ncid,
									 	const char *varname,
									 	const std::vector<int>& wfiles,
									 	int iosysid, 
										int mpirank, int nproc, MPI_Comm comm, 
										uint64_t time_step,
										std::vector<int> &local_proc_blocks, 
										int forced_type = NC_NAT)
{
    /* Read all decomposition blocks assigned to this process,
     * create one big array from them and create a single big
     * decomposition with PIO
     */

    /* Sum the sizes of blocks assigned to this process */
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;

	std::string v_type = bpIO.VariableType(varname);
	if (v_type.empty())
	{
		ierr = BP2PIO_ERROR;
		DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm);
	}

	/* get block merge status of the decomposition */
	std::vector<int> merge_blocks;
    try
    {
		std::string var_name(varname);
		adios2::Variable<int> tmp_base = bpIO.InquireVariable<int>("num_decomp_block_writers/"+var_name);
		const auto v_blocks = bpReader.BlocksInfo(tmp_base, time_step);
		merge_blocks.resize(v_blocks.size());
		for (size_t j = 0; j < v_blocks.size(); j++) {
			std::vector<int> v_value;
		    tmp_base.SetBlockSelection(j);
			bpReader.Get(tmp_base, v_value, adios2::Mode::Sync);
			merge_blocks[j] = (int)v_value[0];
		}
    }
    catch (const std::exception &e)
    {
        ierr = BP2PIO_ERROR;
    }
    catch (...)
    {
        ierr = BP2PIO_ERROR;
    }
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

	printf("MERGE: size: %d\n",merge_blocks.size());
	fflush(stdout);

	*v_base = bpIO.InquireVariable<T>(varname);
    uint64_t nelems = 0;
    try
    {
		if (local_proc_blocks.size()>0) {
			const auto v_blocks = bpReader.BlocksInfo(*v_base, time_step);
			int start_block = 0;
			for (int i = 0; i<local_proc_blocks[0]; i++) {
				start_block += merge_blocks[i];
			}
			for (int i = 0; i<local_proc_blocks.size(); i++) {
				for (int j=0; j<merge_blocks[i]; j++) {
					nelems += v_blocks[start_block].Count[0];
					start_block++;
				}
			}
		}
    }
    catch (const std::exception &e)
    {
        ierr = BP2PIO_ERROR;
    }
    catch (...)
    {
        ierr = BP2PIO_ERROR;
    }
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

	printf("NUM_ELEMS: %d\n",nelems);

    /* Allocate +1 to prevent d_data.data() from returning NULL. Otherwise, read/write operations fail */
    /* nelems may be 0, when some processes do not have any data */
    std::vector<T> d_out;
    d_out.reserve(nelems + 1);
	std::vector<T> v_data;
	try
	{
		if (local_proc_blocks.size()>0) {
			const auto v_blocks = bpReader.BlocksInfo(*v_base, time_step);
			int block_id = 0;
			for (int i = 0; i<local_proc_blocks[0]; i++) {
				block_id += merge_blocks[i];
			}
			for (int i = 0; i<local_proc_blocks.size(); i++) {
				for (int j=0; j<merge_blocks[i]; j++) {
					v_base->SetBlockSelection(block_id);
					adios2::Dims start = {0}, count = {v_blocks[block_id].Count[0]};
					v_base->SetSelection({start, count});
					bpReader.Get(*v_base, v_data, adios2::Mode::Sync);
					d_out.insert(d_out.end(), v_data.begin(), v_data.end());
					block_id++;
				}
			}
		}
	}
	catch (const std::exception &e)
	{
		ierr = BP2PIO_ERROR;
	}
	catch (...)
	{
		ierr = BP2PIO_ERROR;
	}
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

    std::string attname;
    int piotype = forced_type;
    std::string atype;
    AttributeVector adata;
    if (forced_type == NC_NAT)
    {
        attname = string(varname) + "/piotype";
        ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
        if (ierr == BP2PIO_NOERR)
        {
            piotype = *((int*)adata[0].data());
        }
    }
    else
    {
        ierr = BP2PIO_NOERR;
    }
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

    attname = string(varname) + "/ndims";
    ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

    int decomp_ndims = *((int*)adata[0].data());

    attname = string(varname) + "/dimlen";
    ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
    DECOMPOSITION_ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

    int *decomp_dims = (int*)adata[0].data();

	printf("BEFORE: %d\n",nelems);
	fflush(stdout);

    int ioid;
    ret = PIOc_InitDecomp(iosysid, piotype, decomp_ndims, decomp_dims, (PIO_Offset)nelems,
                          (PIO_Offset*)d_out.data(), &ioid, NULL, NULL, NULL);

	printf("ERROR: %d %d %d\n",ret,PIO_NOERR,ioid); fflush(stdout);

    if (ret != PIO_NOERR)
        return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR};

    return Decomposition{ioid, piotype};
}

Decomposition ProcessOneDecomposition(adios2::IO &bpIO, adios2::Engine &bpReader, 
									int ncid,
							  		const char *varname,
							  		const std::vector<int>& wfiles,
							  		int iosysid, 
									int mpirank, int nproc, MPI_Comm comm, 
									uint64_t time_step,
									std::vector<int> &local_proc_blocks,
									int forced_type = NC_NAT)
{
    std::string v_type = bpIO.VariableType(varname);

    if (v_type.empty())
    {
        return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR};
    }

#define declare_template_instantiation(T) \
    else if (v_type == adios2::GetType<T>()) \
    { \
        adios2::Variable<T> v_base; \
        return adios2_ProcessOneDecomposition(&v_base, bpIO, bpReader, ncid, \
                                              varname, wfiles, iosysid, mpirank, nproc, comm, \
											  time_step, local_proc_blocks, forced_type); \
    }

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)

#undef declare_template_instantiation

    return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR};
}

#if 0
void FindDecompositionSteps(adios2::IO &bpIO, adios2::Engine &bpReader,
							uint64_t time_step,
							DecompositionStepMap &decomp_step_map)
{
    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        string v = a2_iter->first;
        if (v.find("/__pio__/decomp/") != string::npos)
        {
            string decompname = v.substr(16); // Skip by strlen("/__pio__/decomp/")
			decomp_step_map[decompname] = time_step;
        }
    }
}

void FindDecompositionInfo(adios2::IO &bpIO, adios2::Engine &bpReader,
							uint64_t time_step,
							DecompositionInfo &decomp_info)
{
    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        string v = a2_iter->first;
        if (v.find("/__pio__/decomp/") != string::npos)
        {
            string decompname = v.substr(16); // Skip by strlen("/__pio__/decomp/")
			decomp_step_map[decompname] = time_step;
        }
    }
}
#endif


void ProcessDecompositions(adios2::IO &bpIO, adios2::Engine &bpReader, 
							int ncid,
							const std::vector<int>& wfiles, 
							int iosysid, 
							MPI_Comm comm, int mpirank, int nproc, 
							uint64_t time_step,
							DecompositionMap &decomp_map,
							std::vector<int> &local_proc_blocks)
{
    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        string v = a2_iter->first;
        if (v.find("num_decomp_block_writers")==string::npos && v.find("/__pio__/decomp/") != string::npos)
        {
			const char *varname = v.c_str();
            string decompname = v.substr(16); // Skip by strlen("/__pio__/decomp/")
            if (!mpirank && debug_out)
                cout << "Process decomposition " << decompname << endl;

			printf("PROCESSING: %s %d\n",v.c_str(),time_step);
			fflush(stdout);

            Decomposition d = ProcessOneDecomposition(bpIO, bpReader, ncid, varname,
                                                      wfiles, iosysid, mpirank, nproc, comm, 
													  time_step, local_proc_blocks);


			printf("DONE_PROCESSING: %s %d\n",v.c_str(),time_step);
			fflush(stdout);

            if (d.ioid == BP2PIO_ERROR)
            {
                throw std::runtime_error("ProcessDecompositions failed.");
            }
            decomp_map[decompname] = d;

			printf("DECOMP: %d\n",d.ioid); fflush(stdout);
        }
        FlushStdout_nm(comm);
    }
}

Decomposition LoadDecomposition(DecompositionMap& decompmap,
							  const string &decompname,
							  adios2::IO &bpIO, adios2::Engine &bpReader, int ncid,
							  const std::vector<int>& wfiles,
							  int nctype, int iosysid,
							  int mpirank, int nproc, MPI_Comm comm,
							  std::string file0, adios2::ADIOS &adios,
							  std::vector<int> &local_proc_blocks)
{
    Decomposition d;

	/* Find the decomp in the file */
	uint64_t time_step = 0;
	while (bpReader.BeginStep()==adios2::StepStatus::OK) {
   		std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
		int found_it = 0;
   		for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
   		{
       		string v = a2_iter->first;
			if (0==v.compare(decompname)) { /* Found it! */
				found_it = 1;
       			d = ProcessOneDecomposition(bpIO, bpReader, ncid, (char*)decompname.c_str(), wfiles,
                  		                    iosysid, mpirank, nproc, comm, time_step, local_proc_blocks, nctype);
       			if (d.ioid == BP2PIO_ERROR)
       			{
           			return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR};
       			}
       			// decompmap[decompname] = d;
				break;
			}
		}
		bpReader.EndStep();
		time_step++;
		if (found_it) {
			bpReader.Close();
			std::string bpIO_name = bpIO.Name();
			adios.RemoveIO(bpIO_name);
            bpIO = adios.DeclareIO(file0 + std::to_string(rand()));
			bpIO.SetParameter("StreamReader","ON");
			bpIO.SetEngine("FileStream");
           	bpReader = bpIO.Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
			break;
		}
	}

    return d;
}

Decomposition GetNewDecomposition(DecompositionMap& decompmap,
							  const string &decompname,
							  adios2::IO &bpIO, adios2::Engine &bpReader, int ncid,
							  const std::vector<int>& wfiles,
							  int nctype, int iosysid,
							  int mpirank, int nproc, MPI_Comm comm,
							  std::string file0, adios2::ADIOS &adios,
							  std::vector<int> &local_proc_blocks)
{
    char ss[PIO_MAX_NAME];
    sprintf(ss, "%s_%d", decompname.c_str(), nctype);
    string key(ss);

    Decomposition d;
    auto it = decompmap.find(key);
    if (it == decompmap.end())
    {
        string varname = "/__pio__/decomp/" + decompname;
		/* Find the decomp in the file */
		uint64_t time_step = 0;
		while (bpReader.BeginStep()==adios2::StepStatus::OK) {
    		std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
			int found_it = 0;
    		for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    		{
        		string v = a2_iter->first;
				if (0==v.compare(varname)) { /* Found it! */
					found_it = 1;
        			d = ProcessOneDecomposition(bpIO, bpReader, ncid, (char*)varname.c_str(), wfiles,
                   				                 iosysid, mpirank, nproc, comm, time_step, local_proc_blocks, nctype);
        			if (d.ioid == BP2PIO_ERROR)
        			{
            			return Decomposition{BP2PIO_ERROR, BP2PIO_ERROR};
        			}
        			decompmap[key] = d;
					break;
				}
			}
			bpReader.EndStep();
			time_step++;
			if (found_it) {
				bpReader.Close();
				std::string bpIO_name = bpIO.Name();
				adios.RemoveIO(bpIO_name);
            	bpIO = adios.DeclareIO(file0 + std::to_string(rand()));
				bpIO.SetParameter("StreamReader","ON");
				bpIO.SetEngine("FileStream");
            	bpReader = bpIO.Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
				break;
			}
		}
    }
    else
    {
        d = it->second;
    }

    return d;
}

void ProcessDimensions(adios2::IO &bpIO, adios2::Engine &bpReader, 
						int ncid, 
						MPI_Comm comm, int mpirank, int nproc,
						DimensionMap &dimensions_map)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret  = PIO_NOERR;

    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        string v = a2_iter->first;
        const char *varname = v.c_str();
        if (v.find("/__pio__/dim/") != string::npos)
        {
            /* For each dimension stored, define a dimension variable with PIO */
            string dimname = v.substr(13); /* 13 = strlen("/__pio__/dim/") */
            if (!mpirank && debug_out)
                cout << "Process dimension " << dimname << endl;

            std::string v_type = bpIO.VariableType(a2_iter->first);
            if (v_type.empty())
            {
                ierr = BP2PIO_ERROR;
            }

#define declare_template_instantiation(T) \
            else if (v_type == adios2::GetType<T>()) \
            { \
                adios2::Variable<T> v_base = bpIO.InquireVariable<T>(varname); \
                std::vector<T> dimval; \
                try \
                { \
                    bpReader.Get(v_base, dimval, adios2::Mode::Sync); \
                } \
                catch (const std::exception &e) \
                { \
                    ierr = BP2PIO_ERROR; \
                } \
                catch (...) \
                { \
                    ierr = BP2PIO_ERROR; \
                } \
                if (ierr == BP2PIO_NOERR) { \
                    int dimid; \
                    PIO_Offset *d_val = (PIO_Offset*)dimval.data(); \
                    ret = PIOc_def_dim(ncid, dimname.c_str(), *d_val, &dimid); \
                    if (ret != PIO_NOERR) \
                    { \
                        throw std::runtime_error("ProcessDimensions failed."); \
                    } \
                    dimensions_map[dimname] = Dimension{dimid, (PIO_Offset)*d_val}; \
                } \
            }

            ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)

#undef declare_template_instantiation
        }
        else
        {
            ierr = BP2PIO_NOERR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessDimensions failed.")
    }
}

void ProcessVariableDefinitions(adios2::IO &bpIO, adios2::Engine &bpReader, 
								int ncid, 
								DimensionMap& dimension_map,
                                MPI_Comm comm, int mpirank, int nproc,
								VariableMap &vars_map,
								std::set<std::string> &var_processed_set)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;

    std::map<std::string, adios2::Params> a2_vi = bpIO.AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
    {
        string v = a2_iter->first;
        if (!mpirank && debug_out)
            cout << "BEFORE Process variable " << v << endl;

		/* check if this variable has been processed before */
		if (var_processed_set.count(v)!=0) {
			continue;
		}
		var_processed_set.insert(v);

        if (v.find("/__") == string::npos)
        {
            /* For each variable written define it with PIO */
            if (v.find("decomp_id/") == string::npos &&
                v.find("frame_id/") == string::npos &&
                v.find("fillval_id/") == string::npos && 
				v.find("num_data_block_writers/") == string::npos &&
				v.find("num_decomp_block_writers/") == string::npos)
            {
                std::string atype;
                AttributeVector adata;

                string attname = v + "/__pio__/nctype";
                ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessVariableDefinitions failed.")
                int nctype = *((int*)adata[0].data());

                attname = v + "/__pio__/ndims";
                ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessVariableDefinitions failed.")
                int ndims = *((int*)adata[0].data());

                int dimids[PIO_MAX_DIMS];
                bool timed = false;
                if (ndims>0)
                {
                    attname = v + "/__pio__/dims";
                    ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
                    ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessVariableDefinitions failed.")

                    for (int d = 0; d < ndims; d++)
                    {
                        dimids[d] = dimension_map[adata[d].data()].dimid;
                        if (dimension_map[adata[d].data()].dimvalue == PIO_UNLIMITED)
                        {
                            timed = true;
                        }
                    }
                }

 				attname = v + "/__pio__/ncop";
			    ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
			    ERROR_CHECK_SINGLE_THROW(ierr, "adios_get_attr_a2 failed.")
			    std::string op(adata[0].data());
				std::string decomp_name = NO_DECOMP;
				if (op == "darray") 
				{
 					attname = v + "/__pio__/decomp";
			    	ierr = adios_get_attr_a2(bpIO, (char*)attname.c_str(), atype, adata);
			    	ERROR_CHECK_SINGLE_THROW(ierr, "adios_get_attr_a2 failed.")
					decomp_name.assign(adata[0].data());
				}

                int varid;
                ret = PIOc_def_var(ncid, (char*)v.c_str(), nctype, ndims, dimids, &varid);
                if (ret != PIO_NOERR)
                    ierr = BP2PIO_ERROR;
                ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessVariableDefinitions failed.")

                vars_map[v] = Variable{varid, timed, nctype, ndims, op, decomp_name, 0};

                ierr = ProcessVarAttributes(bpIO, bpReader, v, ncid, varid, comm);
                ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "ProcessVariableDefinitions failed.")
            }
        }
    }
}

int put_var_nm(int ncid, int varid, int nctype, const std::string &memtype, const void* buf)
{
    int ret = PIO_NOERR;

    if (memtype == adios2::GetType<char>())
    {
        if (nctype == PIO_CHAR)
            ret = PIOc_put_var_text(ncid, varid, (const char*)buf);
        else
            ret = PIOc_put_var_schar(ncid, varid, (const signed char*)buf);
    }
    else if (memtype == adios2::GetType<int8_t>())
    {
        if (nctype == PIO_CHAR)
            ret = PIOc_put_var_text(ncid, varid, (const char*)buf);
        else
            ret = PIOc_put_var_schar(ncid, varid, (const signed char*)buf);
    }
    else if (memtype == adios2::GetType<short int>())
    {
        ret = PIOc_put_var_short(ncid, varid, (const signed short*)buf);
    }
    else if (memtype == adios2::GetType<int16_t>())
    {
        ret = PIOc_put_var_short(ncid, varid, (const signed short*)buf);
    }
    else if (memtype == adios2::GetType<int>())
    {
        ret = PIOc_put_var_int(ncid, varid, (const signed int*)buf);
    }
    else if (memtype == adios2::GetType<int32_t>())
    {
        ret = PIOc_put_var_int(ncid, varid, (const signed int*)buf);
    }
    else if (memtype == adios2::GetType<float>())
    {
        ret = PIOc_put_var_float(ncid, varid, (const float *)buf);
    }
    else if (memtype == adios2::GetType<double>())
    {
        ret = PIOc_put_var_double(ncid, varid, (const double *)buf);
    }
    else if (memtype == adios2::GetType<unsigned char>())
    {
        ret = PIOc_put_var_uchar(ncid, varid, (const unsigned char *)buf);
    }
    else if (memtype == adios2::GetType<uint8_t>())
    {
        ret = PIOc_put_var_uchar(ncid, varid, (const unsigned char*)buf);
    }
    else if (memtype == adios2::GetType<unsigned short>())
    {
        ret = PIOc_put_var_ushort(ncid, varid, (const unsigned short *)buf);
    }
    else if (memtype == adios2::GetType<uint16_t>())
    {
        ret = PIOc_put_var_ushort(ncid, varid, (const unsigned short *)buf);
    }
    else if (memtype == adios2::GetType<unsigned int>())
    {
        ret = PIOc_put_var_uint(ncid, varid, (const unsigned int *)buf);
    }
    else if (memtype == adios2::GetType<uint32_t>())
    {
        ret = PIOc_put_var_uint(ncid, varid, (const unsigned int *)buf);
    }
    else if (memtype == adios2::GetType<long long int>())
    {
        ret = PIOc_put_var_longlong(ncid, varid, (const signed long long *)buf);
    }
    else if (memtype == adios2::GetType<int64_t>())
    {
        ret = PIOc_put_var_longlong(ncid, varid, (const signed long long *)buf);
    }
    else if (memtype == adios2::GetType<unsigned long long int>())
    {
        ret = PIOc_put_var_ulonglong(ncid, varid, (const unsigned long long *)buf);
    }
    else if (memtype == adios2::GetType<uint64_t>())
    {
        ret = PIOc_put_var_ulonglong(ncid, varid, (const unsigned long long *)buf);
    }
    else if (memtype == adios2::GetType<std::string>())
    {
        ret = PIOc_put_var_text(ncid, varid, (const char *)buf);
    }
    else
    {
        /* We can't do anything here, hope for the best, i.e. memtype equals to nctype */
        ret = PIOc_put_var(ncid, varid, buf);
    }

    return ret;
}

int put_vara_nm(int ncid, int varid, int nctype, const std::string &memtype,
                const PIO_Offset *start, const PIO_Offset *count,
                const void* buf)
{
    int ret = PIO_NOERR;

#ifdef _MERGE_SCB
    if (nctype == PIO_BYTE)
    {
       ret = PIOc_put_vara_schar(ncid, varid, start, count, (const signed char*)buf);
    } 
	else if (nctype == PIO_STRING || nctype == PIO_CHAR) 
	{
       ret = PIOc_put_vara_text(ncid, varid, start, count, (const char*)buf);
	}
    else if (nctype == PIO_SHORT)
    {
        ret = PIOc_put_vara_short(ncid, varid, start, count, (const signed short*)buf);
    }
    else if (nctype == PIO_INT)
    {
        ret = PIOc_put_vara_int(ncid, varid, start, count, (const signed int*)buf);
    }
    else if (nctype == PIO_FLOAT || nctype == PIO_REAL)
    {
        ret = PIOc_put_vara_float(ncid, varid, start, count, (const float *)buf);
    }
    else if (nctype == PIO_DOUBLE)
    {
        ret = PIOc_put_vara_double(ncid, varid, start, count, (const double *)buf);
    }
    else if (nctype == PIO_UBYTE)
    {
        ret = PIOc_put_vara_uchar(ncid, varid, start, count, (const unsigned char *)buf);
    }
    else if (nctype == PIO_USHORT)
    {
        ret = PIOc_put_vara_ushort(ncid, varid, start, count, (const unsigned short *)buf);
    }
    else if (nctype == PIO_UINT)
    {
        ret = PIOc_put_vara_uint(ncid, varid, start, count, (const unsigned int *)buf);
    }
    else if (nctype == PIO_INT64)
    {
        ret = PIOc_put_vara_longlong(ncid, varid, start, count, (const signed long long *)buf);
    }
    else if (nctype == PIO_UINT64)
    {
        ret = PIOc_put_vara_ulonglong(ncid, varid, start, count, (const unsigned long long *)buf);
    }
    else
    {
        ret = PIOc_put_vara(ncid, varid, start, count, buf);
    }
#else
    if (memtype == adios2::GetType<char>())
    {
        if (nctype == PIO_BYTE)
            ret = PIOc_put_vara_schar(ncid, varid, start, count, (const signed char*)buf);
        else
            ret = PIOc_put_vara_text(ncid, varid, start, count, (const char*)buf);
    }
    else if (memtype == adios2::GetType<int8_t>())
    {
        if (nctype == PIO_BYTE)
            ret = PIOc_put_vara_schar(ncid, varid, start, count, (const signed char*)buf);
        else
            ret = PIOc_put_vara_text(ncid, varid, start, count, (const char*)buf);
    }
    else if (memtype == adios2::GetType<short>())
    {
        ret = PIOc_put_vara_short(ncid, varid, start, count, (const signed short*)buf);
    }
    else if (memtype == adios2::GetType<int16_t>())
    {
        ret = PIOc_put_vara_short(ncid, varid, start, count, (const signed short*)buf);
    }
    else if (memtype == adios2::GetType<int>())
    {
        ret = PIOc_put_vara_int(ncid, varid, start, count, (const signed int*)buf);
    }
    else if (memtype == adios2::GetType<int32_t>())
    {
        ret = PIOc_put_vara_int(ncid, varid, start, count, (const signed int*)buf);
    }
    else if (memtype == adios2::GetType<float>())
    {
        ret = PIOc_put_vara_float(ncid, varid, start, count, (const float *)buf);
    }
    else if (memtype == adios2::GetType<double>())
    {
        ret = PIOc_put_vara_double(ncid, varid, start, count, (const double *)buf);
    }
    else if (memtype == adios2::GetType<unsigned char>())
    {
        ret = PIOc_put_vara_uchar(ncid, varid, start, count, (const unsigned char *)buf);
    }
    else if (memtype == adios2::GetType<uint8_t>())
    {
        ret = PIOc_put_vara_uchar(ncid, varid, start, count, (const unsigned char *)buf);
    }
    else if (memtype == adios2::GetType<unsigned short>())
    {
        ret = PIOc_put_vara_ushort(ncid, varid, start, count, (const unsigned short *)buf);
    }
    else if (memtype == adios2::GetType<uint16_t>())
    {
        ret = PIOc_put_vara_ushort(ncid, varid, start, count, (const unsigned short *)buf);
    }
    else if (memtype == adios2::GetType<unsigned int>())
    {
        ret = PIOc_put_vara_uint(ncid, varid, start, count, (const unsigned int *)buf);
    }
    else if (memtype == adios2::GetType<uint32_t>())
    {
        ret = PIOc_put_vara_uint(ncid, varid, start, count, (const unsigned int *)buf);
    }
    else if (memtype == adios2::GetType<long long int>())
    {
        ret = PIOc_put_vara_longlong(ncid, varid, start, count, (const signed long long *)buf);
    }
    else if (memtype == adios2::GetType<int64_t>())
    {
        ret = PIOc_put_vara_longlong(ncid, varid, start, count, (const signed long long *)buf);
    }
    else if (memtype == adios2::GetType<unsigned long long int>())
    {
        ret = PIOc_put_vara_ulonglong(ncid, varid, start, count, (const unsigned long long *)buf);
    }
    else if (memtype == adios2::GetType<uint64_t>())
    {
        ret = PIOc_put_vara_ulonglong(ncid, varid, start, count, (const unsigned long long *)buf);
    }
    else if (memtype == adios2::GetType<std::string>())
    {
        ret = PIOc_put_vara_text(ncid, varid, start, count, (const char *)buf);
    }
    else
    {
        ret = PIOc_put_vara(ncid, varid, start, count, buf);
    }
#endif

    return ret;
}

template <class T>
int adios2_ConvertVariablePutVar(adios2::Variable<T>& v_base,
                                 adios2::IO &bpIO, adios2::Engine &bpReader,
								 int ncid, 
                                 const std::string& varname, Variable& var, 
								 uint64_t time_step,
								 MPI_Comm comm, int mpirank, int nproc,
								 int num_bp_writers)
{
    int ret = PIO_NOERR;

    v_base = bpIO.InquireVariable<T>(varname);

    int var_ndims = var.ndims;
    if (var_ndims == 0)
    {
        /* Scalar variable */
        std::vector<T> v_value;
        try
        {
            bpReader.Get(v_base, v_value, adios2::Mode::Sync);
        }
        catch (const std::exception &e)
        {
            return BP2PIO_ERROR;
        }
        catch (...)
        {
            return BP2PIO_ERROR;
        }
        ret = put_var_nm(ncid, var.nc_varid, var.nctype, v_base.Type(), v_value.data());
        if (ret != PIO_NOERR)
        {
            cout << "ERROR in PIOc_put_var(), code = " << ret
                 << " at " << __func__ << ":" << __LINE__ << endl;
            return BP2PIO_ERROR;
        }
    }
    else
    {
        /* An N-dimensional array that needs no rearrangement.
         * put_vara_nm() needs all processes participate */
        TimerStart(write);

        /* E3SM writes this array from I/O processor 0 */
        // PIOc_put_var may have been called multiple times with different start/count values
        // for a variable. We need to convert the output from each of those calls.

		/* Compute the total number of blocks */
		int total_num_blocks = 0;
		// adios2::Variable<T> v_base_tmp = bpIO.InquireVariable<T>(varname);
		const auto v_blocks = bpReader.BlocksInfo(v_base, time_step);
		total_num_blocks += v_blocks.size();

		int var_num_blocks = total_num_blocks/num_bp_writers;
		if (total_num_blocks != (var_num_blocks*num_bp_writers))
		{
			cout << "ERROR in PIOc_put_var(), code = " << ret
				 << " at " << __func__ << ":" << __LINE__ << endl;
			return BP2PIO_ERROR;
		}

		char start_varname[PIO_MAX_NAME];
		char count_varname[PIO_MAX_NAME];
		snprintf(start_varname, PIO_MAX_NAME, "start_id/%s", varname.c_str());
		snprintf(count_varname, PIO_MAX_NAME, "count_id/%s", varname.c_str());

		int elemsize = adios2_type_size_a2(v_base.Type());
		assert(elemsize > 0);

		for (size_t ii = 0; ii < var_num_blocks; ii++)
		{
			try
			{
#ifdef _MERGE_SCB
				std::vector<T> v_data;
				v_base.SetBlockSelection(ii);
				bpReader.Get(v_base, v_data, adios2::Mode::Sync);

				int64_t *pio_var_startp, *pio_var_countp;
				char *data_buf;
				pio_var_startp = (int64_t*)v_data.data();
				pio_var_countp = (int64_t*) ((char*)pio_var_startp + var_ndims*sizeof(int64_t));
				data_buf = (char*)pio_var_countp + var_ndims*sizeof(int64_t);

				PIO_Offset start[var_ndims], count[var_ndims];
				PIO_Offset *start_ptr, *count_ptr;

				if (pio_var_startp[0] < 0) /* NULL start */
				{
					start_ptr = NULL;
				}
				else
				{
					for (size_t d = 0; d < var_ndims; d++)
					{
						start[d] = (PIO_Offset) pio_var_startp[d];
					}
					start_ptr = start;
				}

				if (pio_var_countp[0] < 0) /* NULL count */
				{
					count_ptr = NULL;
					for (size_t d = 0; d < var_ndims; d++)
					{
						pio_var_countp[d] = -1 * (pio_var_countp[d] + 1);
					}
				}
				else
				{
					for (size_t d = 0; d < var_ndims; d++)
					{
						count[d] = (PIO_Offset) pio_var_countp[d];
					}
					count_ptr = count;
				}

				int64_t nelems = 1;
				for (size_t d = 0; d < var_ndims; d++)
				{
					nelems *= pio_var_countp[d];
				}

				ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base.Type(), start_ptr, count_ptr, data_buf);
				if (ret != PIO_NOERR)
				{
					cout << "rank " << mpirank << ":ERROR in PIOc_put_vara(), code = " << ret
						 << " at " << __func__ << ":" << __LINE__ << endl;
					return BP2PIO_ERROR;
				}
#else
				std::vector<int64_t> v_tmp;
				adios2::Variable<int64_t> v_var;
				v_var = bpIO.InquireVariable<int64_t>(start_varname);
				v_var.SetBlockSelection(ii);
				bpReader.Get(v_var, v_tmp, adios2::Mode::Sync);
				memcpy(pio_var_start, v_tmp.data(), sizeof(int64_t)*var_ndims);

				v_var = bpIO.InquireVariable<int64_t>(count_varname);
				v_var.SetBlockSelection(ii);
				bpReader.Get(v_var, v_tmp, adios2::Mode::Sync);
				memcpy(pio_var_count, v_tmp.data(), sizeof(int64_t)*var_ndims);

				PIO_Offset start[var_ndims], count[var_ndims];
				PIO_Offset *start_ptr, *count_ptr;

				/* Check if start was NULL */
				if (pio_var_start[0] < 0) /* NULL start */
				{
					start_ptr = NULL;
				}
				else
				{
					for (size_t d = 0; d < var_ndims; d++)
					{
						start[d] = (PIO_Offset) pio_var_start[d];
					}
					start_ptr = start;
				}

				/* Check if count was NULL */
				if (pio_var_count[0] < 0) /* NULL count */
				{
					count_ptr = NULL;
					for (size_t d = 0; d < var_ndims; d++)
					{
						pio_var_count[d] = -1 * (pio_var_count[d] + 1);
					}
				}
				else
				{
					for (size_t d = 0; d < var_ndims; d++)
					{
						count[d] = (PIO_Offset) pio_var_count[d];
					}
					count_ptr = count;
				}

				int64_t nelems = 1;
				for (size_t d = 0; d < var_ndims; d++)
				{
					nelems *= pio_var_count[d];
				}
				std::vector<char> d(nelems * elemsize);
				v_base.SetBlockSelection(ii);
				std::vector<T> v_data;
				bpReader.Get(v_base, v_data, adios2::Mode::Sync);
				memcpy(d.data(), v_data.data(), nelems * elemsize);

				ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base.Type(), start_ptr, count_ptr, d.data());
				if (ret != PIO_NOERR)
				{
					cout << "rank " << mpirank << ":ERROR in PIOc_put_vara(), code = " << ret
						 << " at " << __func__ << ":" << __LINE__ << endl;
					return BP2PIO_ERROR;
				}
#endif 
			}
			catch (const std::exception &e)
			{
				return BP2PIO_ERROR;
			}
			catch (...)
			{
				return BP2PIO_ERROR;
			}
		}

        TimerStop(write);
    }

    return BP2PIO_NOERR;
}

int ConvertVariablePutVar(adios2::IO &bpIO, 
						adios2::Engine &bpReader,
					  	int ncid, 
					  	const std::string &varname, Variable& var, 
						uint64_t time_step,
						MPI_Comm comm, int mpirank, int nproc,
						int num_bp_writers)
{
    std::string v_type = bpIO.VariableType(varname);
    if (v_type.empty())
    {
        return BP2PIO_ERROR;
    }

#define declare_template_instantiation(T) \
    else if (v_type == adios2::GetType<T>()) \
    { \
        adios2::Variable<T> v_base; \
        return adios2_ConvertVariablePutVar(v_base, bpIO, bpReader, ncid, \
                                            varname, var, time_step, \
											comm, mpirank, nproc, num_bp_writers); \
    }

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)

#undef declare_template_instantiation

    return BP2PIO_ERROR;
}

template <class T>
int adios2_ConvertVariableTimedPutVar(adios2::Variable<T> &v_base,
								adios2::IO &bpIO, adios2::Engine &bpReader, 
								int ncid,
								const std::string &varname, Variable& var, 
								int nblocks_per_step, int time_step,
								MPI_Comm comm, int mpirank, int nproc)
{
    int ret = PIO_NOERR;

    v_base = bpIO.InquireVariable<T>(varname);

    int var_ndims = var.ndims; 
    if (var_ndims==0) /* Scalar variable */
    {
        try
        {
			/* Scalar variable over time */
			/* Written by only one process, so steps = number of blocks in file */
			const auto v_blocks = bpReader.BlocksInfo(v_base, time_step);
			int nsteps = v_blocks.size();

			std::vector<T> v_mins(nsteps);
			for (int ts = 0; ts < nsteps; ts++)
				v_mins[ts] = v_base.Min(ts);

			PIO_Offset start[1], count[1];
			for (int ts = 0; ts < nsteps; ++ts)
			{
				start[0] = ts+var.start_time_step;
				count[0] = 1;
				ret = PIOc_put_vara(ncid, var.nc_varid, start, count, v_mins.data());
				if (ret != PIO_NOERR)
				{
					cout << "ERROR in PIOc_put_vara(), code = " << ret
						 << " at " << __func__ << ":" << __LINE__ << endl;
					return BP2PIO_ERROR;
				}
			}
			var.start_time_step += nsteps;  /* a timed variable may be stored across multiple adios time steps */
        }
        catch (const std::exception &e)
        {
            return BP2PIO_ERROR;
        }
        catch (...)
        {
            return BP2PIO_ERROR;
        }
    }
    else
    {
        /* Calculate how many records/steps we have for this variable */
        int nsteps = 1;
		int ierr = BP2PIO_NOERR;

        /* Compute the total number of blocks */
        int total_num_blocks = 0;
		try
		{
			// adios2::Variable<T> v_base_tmp = bpIO.InquireVariable<T>(varname);
			const auto vb_blocks = bpReader.BlocksInfo(v_base, time_step);
			total_num_blocks += vb_blocks.size();
		}
		catch (const std::exception &e)
		{
			ierr = BP2PIO_ERROR;
            return ierr;
		}
		catch (...)
		{
			ierr = BP2PIO_ERROR;
            return ierr;
		}

        if (var.is_timed)
        {
            nsteps = total_num_blocks / nblocks_per_step;
        }

        if (total_num_blocks != (nsteps * nblocks_per_step))
        {
            if (debug_out)
                cout << "rank " << mpirank << ":ERROR in processing variable '" << varname
                     << "'. Number of blocks = " << total_num_blocks
                     << " does not equal the number of steps * number of writers = "
                     << nsteps << " * " << nblocks_per_step << " = " << nsteps * nblocks_per_step
                     << endl;
        }

        char start_varname[PIO_MAX_NAME];
        char count_varname[PIO_MAX_NAME];
        snprintf(start_varname, PIO_MAX_NAME, "start_id/%s", varname.c_str());
        snprintf(count_varname, PIO_MAX_NAME, "count_id/%s", varname.c_str());

        /* Just read the arrays written by rank 0 (on every process here) and
         * write it collectively.
        */
        for (int ts = 0; ts < nsteps; ++ts)
        {
            TimerStart(read);

            try
            {
                int elemsize = adios2_type_size_a2(v_base.Type());
                assert(elemsize > 0);

#ifdef _MERGE_SCB
                std::vector<T> v_data;
				v_base.SetBlockSelection(ts);
                bpReader.Get(v_base, v_data, adios2::Mode::Sync);

				int64_t *pio_var_startp, *pio_var_countp;
				char *data_buf;
				pio_var_startp = (int64_t*)v_data.data();
				pio_var_countp = (int64_t*) ((char*)pio_var_startp + var_ndims*sizeof(int64_t));
				data_buf = (char*)pio_var_countp + var_ndims*sizeof(int64_t);

				PIO_Offset start[var_ndims], count[var_ndims];
				PIO_Offset *start_ptr, *count_ptr;

                if (pio_var_startp[0] < 0) /* NULL start */
                {
                    start_ptr = NULL;
                }
                else
                {
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        start[d] = (PIO_Offset) pio_var_startp[d];
                    }
                    start_ptr = start;
                }

                if (pio_var_countp[0] < 0) /* NULL count */
                {
                    count_ptr = NULL;
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        pio_var_countp[d] = -1 * (pio_var_countp[d] + 1);
                    }
                }
                else
                {
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        count[d] = (PIO_Offset) pio_var_countp[d];
                    }
                    count_ptr = count;
                }

				TimerStop(read);

                TimerStart(write);

                ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base.Type(), start_ptr, count_ptr, data_buf);
                if (ret != PIO_NOERR)
                {
                    cout << "ERROR in PIOc_put_vara(), code = " << ret
                         << " at " << __func__ << ":" << __LINE__ << endl;
                    return BP2PIO_ERROR;
                }

                TimerStop(write);
#else
                std::vector<int64_t> v_tmp;
                adios2::Variable<int64_t> v_var;
                v_var = bpIO.InquireVariable<int64_t>(start_varname);
                v_var.SetBlockSelection(ts);
                bpReader.Get(v_var, v_tmp, adios2::Mode::Sync);
                memcpy(pio_var_start, v_tmp.data(), sizeof(int64_t)*var_ndims);

                v_var = bpIO.InquireVariable<int64_t>(count_varname);
                v_var.SetBlockSelection(ts);
                bpReader.Get(v_var, v_tmp, adios2::Mode::Sync);
                memcpy(pio_var_count, v_tmp.data(), sizeof(int64_t)*var_ndims);

                PIO_Offset start[var_ndims], count[var_ndims];
                PIO_Offset *start_ptr, *count_ptr;
                if (pio_var_start[0] < 0) /* NULL start */
                {
                    start_ptr = NULL;
                }
                else
                {
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        start[d] = (PIO_Offset) pio_var_start[d];
                    }
                    start_ptr = start;
                }

                if (pio_var_count[0] < 0) /* NULL count */
                {
                    count_ptr = NULL;
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        pio_var_count[d] = -1 * (pio_var_count[d] + 1);
                    }
                }
                else
                {
                    for (size_t d = 0; d < var_ndims; d++)
                    {
                        count[d] = (PIO_Offset) pio_var_count[d];
                    }
                    count_ptr = count;
                }

                int64_t nelems = 1;
                for (size_t d = 1; d < var_ndims; d++)
                {
                    nelems *= pio_var_count[d];
                }
                std::vector<char> d(nelems * elemsize);
                v_base.SetBlockSelection(ts);
                std::vector<T> v_data;
                bpReader.Get(v_base, v_data, adios2::Mode::Sync);
                memcpy(d.data(), v_data.data(), nelems * elemsize);

                TimerStop(read);

                TimerStart(write);

                ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base.Type(), start_ptr, count_ptr, d.data());
                if (ret != PIO_NOERR)
                {
                    cout << "ERROR in PIOc_put_vara(), code = " << ret
                         << " at " << __func__ << ":" << __LINE__ << endl;
                    return BP2PIO_ERROR;
                }

                TimerStop(write);
#endif
            }
            catch (const std::exception &e)
            {
                return BP2PIO_ERROR;
            }
            catch (...)
            {
                return BP2PIO_ERROR;
            }
        }
		var.start_time_step += nsteps;  /* a timed variable may be stored across multiple adios time steps */
	}

    return BP2PIO_NOERR;
}

int ConvertVariableTimedPutVar(adios2::IO &bpIO, adios2::Engine &bpReader,
						   	int ncid, 
						   	const std::string &varname, Variable& var, 
							int nblocks_per_step, int time_step,
						   	MPI_Comm comm, int mpirank, int nproc)
{
    std::string v_type = bpIO.VariableType(varname);
    if (v_type.empty())
    {
        return BP2PIO_ERROR;
    }

#define declare_template_instantiation(T) \
    else if (v_type == adios2::GetType<T>()) \
    { \
        adios2::Variable<T> v_base; \
        return adios2_ConvertVariableTimedPutVar(v_base,bpIO, bpReader, ncid, \
                                                 varname, var, nblocks_per_step, \
                                                 time_step, comm, mpirank, nproc); \
    }

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)

#undef declare_template_instantiation

    return BP2PIO_ERROR;
}

template <class T>
int adios2_ConvertVariableDarray(adios2::Variable<T> *v_base, std::vector<T> v_value,
							 	IOVector &bpIO, EngineVector &bpReader, std::string varname,
							 	int ncid, Variable& var,
							 	const std::vector<int>& wfiles,
								DecompositionMap& decomp_map,
							 	int nblocks_per_step, int iosysid,
								std::string file0, adios2::ADIOS &adios, uint64_t time_step, 
							 	MPI_Comm comm, int mpirank, int nproc, int mem_opt, std::vector<int> &local_proc_blocks)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;

	/* Different decompositions at different frames */
    char decomp_varname[PIO_MAX_NAME];
    char frame_varname[PIO_MAX_NAME];
    char fillval_varname[PIO_MAX_NAME];
    sprintf(decomp_varname, "decomp_id/%s", varname.c_str());
    sprintf(frame_varname, "frame_id/%s", varname.c_str());
    sprintf(fillval_varname, "fillval_id/%s", varname.c_str());
    int  decomp_id, frame_id, fillval_exist;
    char decompname[PIO_MAX_NAME];
    char fillval_id[PIO_MAX_NAME];

    /* Calculate how many records/steps we have for this variable */
    /* Compute the number of time steps */
	int nsteps = 0;
	try
	{
		std::cout << "TIME STEP: " << time_step << " VARIABLE: " << frame_varname << std::endl;
		adios2::Variable<int32_t> f_var = bpIO[0].InquireVariable<int32_t>(frame_varname);
		const auto vb_blocks = bpReader[0].BlocksInfo(f_var, time_step);
		nsteps = (int)vb_blocks.size();
	}
	catch (const std::exception &e)
	{
		std::cout << "TIME STEPS Error: " << e.what() << std::endl;
		ierr = BP2PIO_ERROR;
	}
	catch (...)
	{
		ierr = BP2PIO_ERROR;
	}
    ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)

	/* Find the number of merge_blocks info per time step */
	int num_merge_blocks_per_step = 0;
	std::vector<int> merge_blocks;
	adios2::Variable<int> tmp_base = bpIO[0].InquireVariable<int>("num_data_block_writers/"+varname);
	{
		const auto v_blocks = bpReader[0].BlocksInfo(tmp_base, time_step);
		num_merge_blocks_per_step = v_blocks.size()/nsteps;
		merge_blocks.resize(num_merge_blocks_per_step);
	}

	printf("DARRAY_SIZE: %d\n",merge_blocks.size()); fflush(stdout);

	int start_merge_blocks_id = 0;
	int start_data_blocks_id  = 0;
	*v_base = bpIO[0].InquireVariable<T>(varname);
	int elemsize = adios2_type_size_a2(v_base->Type());
	assert(elemsize > 0);
    for (int ts = 0; ts < nsteps; ++ts)
    {
        try
        {
			if (mpirank==0) {
				printf("Converting: %s time: %d %d %d\n",varname.c_str(),ts,nsteps,nproc); 
				fflush(stdout);
			}

			/* get block merge status of the distributed array */
			for (int i = 0; i<local_proc_blocks.size(); i++) {
				std::vector<int> v_value;
				tmp_base.SetBlockSelection(i+start_merge_blocks_id);
				bpReader[0].Get(tmp_base, v_value, adios2::Mode::Sync);
				merge_blocks[i] = (int)v_value[0];
			}
			start_merge_blocks_id += num_merge_blocks_per_step;

			const auto vb_blocks = bpReader[0].BlocksInfo(*v_base, time_step);
            uint64_t nelems = 0;
			{
				/* Allocate +1 to prevent d_data.data() from returning NULL. Otherwise, read/write operations fail */
				/* nelems may be 0, when some processes do not have any data */
    			size_t blockid = 0;
				nelems = 0;
				int start_block = 0;
				for (int i = 0; i<local_proc_blocks[0]; i++) {
					start_block += merge_blocks[i];
				}
				for (int i = 0; i<local_proc_blocks.size(); i++) {
					for (int j=0; j<merge_blocks[i]; j++) {
						blockid = start_block+start_data_blocks_id;
						nelems += vb_blocks[blockid].Count[0];
						start_block++;
					}
				}
				printf("DARRAY NELEMS: %d\n",nelems); fflush(stdout); 
				std::vector<char> d_data((nelems + 1) * elemsize);

				/* read in the data array */
				uint64_t offset = 0;
				start_block = 0;
				for (int i = 0; i<local_proc_blocks[0]; i++) {
					start_block += merge_blocks[i];
				}
				for (int i = 0; i<local_proc_blocks.size(); i++) {
					for (int j=0; j<merge_blocks[i]; j++) {
						blockid = start_block+start_data_blocks_id;
						v_base->SetBlockSelection(blockid);
						bpReader[0].Get(*v_base, v_value, adios2::Mode::Sync);
						memcpy(d_data.data() + offset, v_value.data(), v_value.size()*elemsize);
						offset += vb_blocks[blockid].Count[0] * elemsize;
						start_block++;
					}
				}
				int sum_blocks = 0;
				for (int i = 0; i<merge_blocks.size(); i++) {
					sum_blocks += merge_blocks[i];
				}
				start_data_blocks_id += sum_blocks;
					

				std::vector<int32_t> v_tmp;
				adios2::Variable<int32_t> v_var;
				v_var = bpIO[0].InquireVariable<int32_t>(decomp_varname);
				v_var.SetBlockSelection(ts);
				bpReader[0].Get(v_var, v_tmp, adios2::Mode::Sync);
				memcpy(&decomp_id, v_tmp.data(), sizeof(int32_t));
	
				v_var = bpIO[0].InquireVariable<int32_t>(frame_varname);
				v_var.SetBlockSelection(ts);
				bpReader[0].Get(v_var, v_tmp, adios2::Mode::Sync);
				memcpy(&frame_id, v_tmp.data(), sizeof(int32_t));

				/* Fix for NUM_FRAMES */
				if (!var.is_timed && frame_id >= 0)
					var.is_timed = true;

				if (decomp_id > 0)
				{
					std::vector<T> v1_tmp;
					adios2::Variable<T> v1_var = bpIO[0].InquireVariable<T>(fillval_varname);
					v1_var.SetBlockSelection(ts);
					bpReader[0].Get(v1_var, v1_tmp, adios2::Mode::Sync);
					memcpy(fillval_id, v1_tmp.data(), v1_tmp.size()*sizeof(T));
					fillval_exist = 1;
				}
				else
				{
					decomp_id = -decomp_id;
					fillval_exist = 0;
				}

				Decomposition decomp;
				if (mem_opt) 
				{
					sprintf(decompname, "/__pio__/decomp/%d", decomp_id);
					decomp = LoadDecomposition(decomp_map, decompname, bpIO[1], bpReader[1],
						    				   ncid, wfiles, NC_NAT, iosysid, mpirank, nproc, comm, 
											   file0, adios, local_proc_blocks);
				}
				else
				{
					sprintf(decompname, "%d", decomp_id);
					decomp = decomp_map[decompname];
				}

				if (decomp.ioid == BP2PIO_ERROR)
				{
					ierr = BP2PIO_ERROR;
					break;
				}
						
				if (decomp.piotype != var.nctype)
				{
				   /* Type conversion may happened at writing. Now we make a new decomposition for this nctype */
				   if (mem_opt)
				   {
						ret = PIOc_freedecomp(iosysid, decomp.ioid);
						if (ret != PIO_NOERR)
						{
							ierr = BP2PIO_ERROR;
							break;
						}
						sprintf(decompname, "/__pio__/decomp/%d", decomp_id);
						decomp = LoadDecomposition(decomp_map, decompname, bpIO[1], bpReader[1],
													ncid, wfiles, var.nctype, iosysid, mpirank, nproc, comm, 
													file0,adios,local_proc_blocks);
					}
					else
					{
						sprintf(decompname, "%d", decomp_id);
						decomp = GetNewDecomposition(decomp_map, decompname, bpIO[1], bpReader[1],
													ncid, wfiles, var.nctype, iosysid, mpirank, nproc, comm, 
													file0,adios,local_proc_blocks);
					}

					if (decomp.ioid == BP2PIO_ERROR)
					{
						ierr = BP2PIO_ERROR;
						break;
					}
				}

				if (frame_id < 0)
					frame_id = 0;

				/* Different decompositions at different frames */
				/* Note: this variable can have an unlimited or limited time dimension */
				if (var.is_timed)
				{
					ret = PIOc_setframe(ncid, var.nc_varid, frame_id);
					if (ret != PIO_NOERR)
					{
						ierr = BP2PIO_ERROR;
						break;
					}
				}

				if (fillval_exist)
				{
					ret = PIOc_write_darray(ncid, var.nc_varid, decomp.ioid, (PIO_Offset)nelems,
											d_data.data(), fillval_id);
				}
				else
				{
					ret = PIOc_write_darray(ncid, var.nc_varid, decomp.ioid, (PIO_Offset)nelems,
											d_data.data(), NULL);
				}

				if (ret != PIO_NOERR)
				{
					ierr = BP2PIO_ERROR;
					break;
				}

				if (mem_opt)
				{
					ret = PIOc_sync(ncid);
					if (ret != PIO_NOERR)
					{
						ierr = BP2PIO_ERROR;
						break;
					}

					ret = PIOc_freedecomp(iosysid, decomp.ioid);
					if (ret != PIO_NOERR)
					{
						ierr = BP2PIO_ERROR;
						break;
					}
				}
			}
		}
		catch (const std::exception &e)
		{
			std::cout << "Error: " << e.what() << " Timestep: " << ts << " " << time_step<< std::endl;
			ierr = BP2PIO_ERROR;
			break;
		}
		catch (...)
		{
			ierr = BP2PIO_ERROR;
			break;
		}
	}
    ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm)
	
	if (mpirank==0) {
		printf("END Variable: %s nsteps: %d nblocks_per_step: %d time_step: %d\n",varname.c_str(),nsteps,nblocks_per_step,time_step);
		fflush(stdout);
	}

    return BP2PIO_NOERR;
}

int ConvertVariableDarray(IOVector &bpIO, EngineVector &bpReader,
					  	std::string varname,
					  	int ncid, Variable& var,
					  	const std::vector<int>& wfiles,
					  	DecompositionMap& decomp_map,
					  	int nblocks_per_step, int iosysid,
						std::string file0, adios2::ADIOS &adios, int time_step,
					  	MPI_Comm comm, int mpirank, int nproc, int mem_opt, std::vector<int> &local_proc_blocks)
{
    std::string v_type = bpIO[0].VariableType(varname);
    if (v_type.empty())
    {
        return BP2PIO_ERROR;
    }

#define declare_template_instantiation(T) \
    else if (v_type == adios2::GetType<T>()) \
    { \
        adios2::Variable<T> v_base; \
        std::vector<T> v_value; \
        return adios2_ConvertVariableDarray(&v_base, v_value, \
                                            bpIO, bpReader, varname, ncid, var, \
                                            wfiles, decomp_map, nblocks_per_step, iosysid, \
                                            file0, adios, time_step, comm, mpirank, nproc, mem_opt, local_proc_blocks); \
    }

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)

#undef declare_template_instantiation

    return BP2PIO_ERROR;
}

/*
 * Assumes a BP folder with name "infilename.dir" and
 * all the files in the folder are bp files. It also
 * assumes the file extensions are "infilename.bp.X"
 * where X is 0 to N-1.
 */
int GetNumOfFiles(const string &infilename)
{
    int file_count = 0;
    string foldername = infilename; // + ".dir/";

    DIR* dirp = opendir(foldername.c_str());
    if (!dirp)
    {
        fprintf(stderr, "Folder %s does not exist.\n", foldername.c_str());
        return BP2PIO_ERROR;
    }

    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL)
    {
        if (strstr(dp->d_name, "data.") != NULL)
            file_count++;
    }

    closedir(dirp);

    return file_count;
}

std::string ExtractFilename(const std::string &pathname)
{
    size_t pos = pathname.find_last_of("/\\");
    if (pos == std::string::npos)
    {
        return pathname;
    }
    else
    {
        return pathname.substr(pos + 1);
    }
}

std::string ExtractPathname(const std::string &pathname)
{
    size_t pos = pathname.find_last_of("/\\");
    if (pos == std::string::npos)
    {
        return "./";
    }
    else
    {
        return pathname.substr(0, pos);
    }
}

std::vector<int> FindProcessBlockGroupAssignments(std::vector<int> &block_procs,
							int mpirank, int nproc, MPI_Comm comm) 
{
	int num_blocks = block_procs.size();
	int nwb = num_blocks / nproc;
    int start_wb;

	if (mpirank < num_blocks % nproc)
		nwb++;

    if (mpirank < num_blocks % nproc)
        start_wb = mpirank * nwb;
    else
        start_wb = mpirank * nwb + num_blocks % nproc;

    std::vector<int> blocks;
	if (nwb>0) {
    	blocks.resize(nwb);
   		for (int i = 0; i < nwb; ++i)
       		blocks[i] = start_wb + i;
	}

    return blocks;
}

int ConvertBPFile(const string &infilepath, const string &outfilename,
                    int pio_iotype, int iosysid,
                    MPI_Comm comm, int mpirank, int nproc, int mem_opt)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    string err_msg = "No errors";
    int num_infiles = 0;
    int ncid = -1;
    int n_bp_writers;
    int ret = PIO_NOERR;

    MPI_Comm w_comm = comm;
    int w_mpirank, w_nproc;

    MPI_Comm_rank(w_comm, &w_mpirank);
    MPI_Comm_size(w_comm, &w_nproc);

    // Initialization of the class factory
    adios2::ADIOS adios(w_comm, adios2::DebugON);

    try
    {
        /*
         * Get the number of files in BP folder.
         * This operation assumes that the BP folder contains only the
         * BP files.
         */
        int n_bp_files = GetNumOfFiles(infilepath);
        if (n_bp_files == BP2PIO_ERROR)
        {
            err_msg = "Input folder doesn't exist.";
            ierr = BP2PIO_ERROR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg)

        if (nproc > n_bp_files)
        {
            if (debug_out)
                cout << "ERROR: nproc (" << nproc << ") is greater than #files ("
                     << n_bp_files << ")" << endl;
            ierr = BP2PIO_ERROR;
            err_msg =  "Use fewer processors.";
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg)

        /* Number of BP file writers != number of converter processes here */
        std::vector<int> wfiles;
        wfiles = AssignWriteRanks(n_bp_files, comm, mpirank, nproc);
        if (debug_out)
        {
            for (auto nb: wfiles)
                printf("Myrank: %d File id: %d\n", mpirank, nb);
        }
        num_infiles = wfiles.size() + 1; // infile[0] is opened by all processors

        /* Now allocate IO and Engine for each file */
        std::vector<adios2::IO> bpIO(num_infiles);
        std::vector<adios2::Engine> bpReader(num_infiles);

		// Open BP4 file
        string file0 = infilepath; 
        try
        {
			bpIO[0] = adios.DeclareIO(file0 + "_0");
			bpIO[0].SetParameter("StreamReader","ON");
			bpIO[0].SetEngine("FileStream");
			bpReader[0] = bpIO[0].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);

			/*
			 * Extra bpIO and bpReader are used for going over ADIOS steps 
			 * multiple times for decompositions, variables, etc.         
			 */
			bpIO[1] = adios.DeclareIO(file0 + "_1");
			bpIO[1].SetParameter("StreamReader","ON");
			bpIO[1].SetEngine("FileStream");
			bpReader[1] = bpIO[1].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
        }
        catch (const std::exception &e)
        {
            err_msg =  e.what();
            ierr = BP2PIO_ERROR;
        }
        catch (...)
        {
            err_msg = "Unknown exception.";
            ierr = BP2PIO_ERROR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg)

		/* Process nproc and block procs objects */
		std::vector<int> block_procs;
		uint64_t time_step = 0;
		while (bpReader[0].BeginStep()==adios2::StepStatus::OK) 
		{
    		adios2::Variable<int> bpNProc = bpIO[0].InquireVariable<int>("/__pio__/info/nproc");
			if (bpNProc) {
        		try {
            		bpReader[0].Get(bpNProc, &n_bp_writers, adios2::Mode::Sync);
        		} catch (const std::exception &e) {
					err_msg = e.what();
					ierr = BP2PIO_ERROR;
				} catch (...) {
					err_msg = "Unknown exception.";
					ierr = BP2PIO_ERROR;
				}
				ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg);
			}

			adios2::Variable<int> blockProcs = bpIO[0].InquireVariable<int>("/__pio__/info/block_nprocs");
			if (blockProcs) {
				const auto v_blocks = bpReader[0].BlocksInfo(blockProcs, time_step);
				block_procs.resize(v_blocks.size());
				for (int ii=0;ii<v_blocks.size();ii++) {
        			std::vector<int> v_value;
					blockProcs.SetBlockSelection(ii);
					bpReader[0].Get(blockProcs, v_value, adios2::Mode::Sync);
					block_procs[ii] = (int)v_value[0];
				}
			}
			bpReader[0].EndStep();
			time_step++;
		}

		printf("BLOCK: size: %d\n",block_procs.size());
    	std::cout << std::flush;
		MPI_Barrier(comm);

		std::vector<int> local_proc_blocks = FindProcessBlockGroupAssignments(block_procs,mpirank,nproc,comm);
		MPI_Barrier(comm);

		for (int i=0;i<local_proc_blocks.size();i++) {
			printf("LOCAL myrank: %d local_proc_block: %d\n",mpirank,local_proc_blocks[i]);
    		std::cout << std::flush;
		}

    	int io_proc = 0;
    	int split_comm = 0;
    	if (block_procs.size() < w_nproc)
    	{
        	if (w_mpirank == 0)
        	{
           	 	fprintf(stderr, "Warning: #files: %d < #procs: %d\n", block_procs.size(), w_nproc);
           	 	fflush(stderr);
        	}

			if (local_proc_blocks.size()>0) /* I/O process */ 
			{
           	 	io_proc = 1;
			}

        	MPI_Comm_split(w_comm, io_proc, w_mpirank, &comm);
        	split_comm = 1;
        	MPI_Comm_rank(comm, &mpirank);
        	MPI_Comm_size(comm, &nproc);
    	}
    	else
    	{
       	 	io_proc = 1;
       	 	comm    = w_comm;
       	 	mpirank = w_mpirank;
       	 	nproc   = w_nproc;
    	}

		if (io_proc==0) { /* Not I/O process */
			MPI_Barrier(w_comm);
			return BP2PIO_NOERR;
		}

		iosysid = InitPIO(comm, mpirank, nproc);
		if (iosysid == BP2PIO_ERROR)
		{
			ierr = BP2PIO_ERROR;
			ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg);
		}

		/* Create output file */
        /*
         *   Use NC_64BIT_DATA instead of PIO_64BIT_OFFSET. Some output files will have variables
         *   that require more than 4GB storage.
         */
        ret = PIOc_createfile(iosysid, &ncid, &pio_iotype, outfilename.c_str(), NC_64BIT_DATA);
        if (ret != PIO_NOERR)
            ierr = BP2PIO_ERROR;
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "Could not create output file.");

		MPI_Barrier(comm);

		/* Reopen the file to reset steps */
		ierr = BP2PIO_NOERR;
		try {
			bpReader[0].Close();
			std::string bpIO_name = bpIO[0].Name();
			adios.RemoveIO(bpIO_name);
            bpIO[0] = adios.DeclareIO(file0 + std::to_string(rand()));
			bpIO[0].SetParameter("StreamReader","ON");
			bpIO[0].SetEngine("FileStream");
            bpReader[0] = bpIO[0].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
        } catch (const std::exception &e) {
            err_msg = e.what();
            ierr = BP2PIO_ERROR;
        } catch (...) {
            err_msg = "Unknown exception.";
            ierr = BP2PIO_ERROR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg);

		/* Process one time objects */
		time_step = 0;
		while (bpReader[0].BeginStep()==adios2::StepStatus::OK) 
		{
	        /* Process the global fillmode */
        	ierr = ProcessGlobalFillmode(bpIO[0], ncid, comm);	
			bpReader[0].EndStep();
			time_step++;
		}

		/* Reopen the file to reset steps */
		ierr = BP2PIO_NOERR;
		try {
			bpReader[0].Close();
			std::string bpIO_name = bpIO[0].Name();
			adios.RemoveIO(bpIO_name);
            bpIO[0] = adios.DeclareIO(file0 + std::to_string(rand()));
			bpIO[0].SetParameter("StreamReader","ON");
			bpIO[0].SetEngine("FileStream");
            bpReader[0] = bpIO[0].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
        } catch (const std::exception &e) {
            err_msg = e.what();
            ierr = BP2PIO_ERROR;
        } catch (...) {
            err_msg = "Unknown exception.";
            ierr = BP2PIO_ERROR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg);

		/* Process dimensions, decomposition arrays, variable definitions, and global attributes */
        DimensionMap dimension_map;
        DecompositionMap decomp_map;
		VariableMap vars_map;
    	std::map<std::string, char> processed_attrs;
    	std::map<std::string, int> var_att_map;
		std::set<std::string> var_processed_set;
		time_step = 0;
		while (bpReader[0].BeginStep()==adios2::StepStatus::OK) 
		{
        	/* Process dimensions */
        	ProcessDimensions(bpIO[0], bpReader[0], ncid, comm, mpirank, nproc, dimension_map);

        	/* Process decompositions */
           	ProcessDecompositions(bpIO[0], bpReader[0], ncid, wfiles, iosysid, comm, mpirank, nproc, time_step, decomp_map, local_proc_blocks);

        	/* For each variable, define a variable with PIO */
        	ProcessVariableDefinitions(bpIO[0], bpReader[0], ncid, dimension_map, comm, mpirank, nproc, vars_map, 
									var_processed_set);

        	/* Process the global attributes */
        	ierr = ProcessGlobalAttributes(bpIO[0], bpReader[0], ncid, dimension_map, vars_map, comm, processed_attrs, var_att_map);

			time_step++;
			bpReader[0].EndStep();
		}

        ret = PIOc_enddef(ncid);
        if (ret != PIO_NOERR)
            ierr = BP2PIO_ERROR;
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "PIOc_enddef failed.");

		/* Reset time steps */
		try
        {
			bpReader[0].Close();
			std::string bpIO_name = bpIO[0].Name();
			adios.RemoveIO(bpIO_name);
            bpIO[0] = adios.DeclareIO(file0 + std::to_string(rand()));
			bpIO[0].SetParameter("StreamReader","ON");
			bpIO[0].SetEngine("FileStream");
            bpReader[0] = bpIO[0].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);

			bpReader[1].Close();
			bpIO_name = bpIO[1].Name();
			adios.RemoveIO(bpIO_name);
            bpIO[1] = adios.DeclareIO(file0 + std::to_string(rand()));
			bpIO[1].SetParameter("StreamReader","ON");
			bpIO[1].SetEngine("FileStream");
            bpReader[1] = bpIO[1].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);
        }
        catch (const std::exception &e)
        {
            err_msg =  e.what();
            ierr = BP2PIO_ERROR;
        }
        catch (...)
        {
            err_msg = "Unknown exception.";
            ierr = BP2PIO_ERROR;
        }
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg)

		/* Convert variables */
		time_step = 0;
		while (bpReader[0].BeginStep()==adios2::StepStatus::OK) 
		{
        	std::map<std::string, adios2::Params> a2_vi = bpIO[0].AvailableVariables();
			for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter != a2_vi.end(); ++a2_iter)
			{
				string v = a2_iter->first;
				if (v.find("/__") == string::npos)
				{
					printf("CONVERTING: %s\n",v.c_str()); fflush(stdout);

					if (v.find("decomp_id/") == string::npos &&
						v.find("frame_id/") == string::npos &&
						v.find("fillval_id/") == string::npos && 
						v.find("num_data_block_writers/") == string::npos &&
						v.find("num_decomp_block_writers/") == string::npos)
					{
						/* For each variable, read with ADIOS then write with PIO */
						if (!mpirank && debug_out)
							cout << "Convert variable: " << v << endl;

						Variable& var = vars_map[v];
						printf("NOW_CONVERTING: %s %s\n",v.c_str(),var.op.c_str()); fflush(stdout);

						if (var.op == "put_var")
						{
							if (var.is_timed)
							{
								if (debug_out)
								{
									printf("ConvertVariableTimedPutVar: %d\n", mpirank);
									fflush(stdout);
								}
								ierr = ConvertVariableTimedPutVar(bpIO[0], bpReader[0], ncid, v, var,
																  n_bp_writers, time_step, 
																  comm, mpirank, nproc);
								ERROR_CHECK_SINGLE_THROW(ierr, "ConvertVariableTimedPutVar failed.")
							} 
							else
							{
								if (debug_out)
								{
									printf("ConvertVariablePutVar: %d\n", mpirank);
									fflush(stdout);
								}
								ierr = ConvertVariablePutVar(bpIO[0], bpReader[0], ncid, v, var, time_step, 
															comm,  mpirank, nproc, n_bp_writers);
								ERROR_CHECK_SINGLE_THROW(ierr, "ConvertVariablePutVar failed.")
							}
						} 
						else if (var.op == "darray")
						{
							/* Variable was written with pio_write_darray() with a decomposition */
							if (debug_out)
							{
								printf("ConvertVariableDarray: %d\n", mpirank);
								fflush(stdout);
							}
							printf("CONVERT DARRAY: %s\n",v.c_str());
							fflush(stdout);
							ierr = ConvertVariableDarray(bpIO, bpReader, v, ncid, var, wfiles, decomp_map,
												n_bp_writers, iosysid, file0, adios, time_step, 
												comm, mpirank, nproc, mem_opt, local_proc_blocks);
							ERROR_CHECK_SINGLE_THROW(ierr, "ConvertVariableDarray failed.")
						}
					}
				}

				ret = PIOc_sync(ncid); /* FIXME: flush after each variable until development is done. Remove for efficiency */
				if (ret != PIO_NOERR)
					ierr = BP2PIO_ERROR;
				ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "PIOc_sync failed.");
			}

			MPI_Barrier(comm);
			printf("END TIME STEP: %d\n",time_step); fflush(stdout);
			MPI_Barrier(comm);

			time_step++;
			bpReader[0].EndStep();
		}

		printf("FINISHING %d ...\n",time_step); fflush(stdout);

#if 0
		/* Finish up */
        for (std::map<std::string, Decomposition>::iterator it = decomp_map.begin();
				it != decomp_map.end(); ++it)
        {
			Decomposition d = it->second;
			ret = PIOc_freedecomp(iosysid, d.ioid);
			if (ret != PIO_NOERR)
				ierr = BP2PIO_ERROR;
			ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "PIOc_freedecomp failed.");
		}
#endif 

        ret = PIOc_sync(ncid);
        if (ret != PIO_NOERR)
            ierr = BP2PIO_ERROR;
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "PIOc_sync failed.");

        ret = PIOc_closefile(ncid);
        if (ret != PIO_NOERR)
            ierr = BP2PIO_ERROR;
        ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, "PIOc_closefile failed.");

		ret = PIOc_finalize(iosysid);
        if (ret != PIO_NOERR)
        {
           throw std::runtime_error("PIOc_finalize error.");
        }
        TimerStop(write);
    }
    catch (const std::exception &e)
    {
		std::cout << "ADIOS ERROR: " << e.what() << std::endl;
        err_msg = e.what();
        ierr = BP2PIO_ERROR;
    }
    catch (...)
    {
        err_msg = "Unknown exception.";
        ierr = BP2PIO_ERROR;
    }
    ERROR_CHECK_THROW(ierr, err_val, err_cnt, comm, err_msg)



	MPI_Barrier(w_comm);
    return BP2PIO_NOERR;
}

enum PIO_IOTYPE GetIOType_nm(const string &t)
{
    enum PIO_IOTYPE iotype = PIO_IOTYPE_NETCDF;
    if (t == "pnetcdf" || t == "PNETCDF" || t == "1")
    {
        iotype = PIO_IOTYPE_PNETCDF;
    }
    else if (t == "netcdf" || t == "NETCDF" || t == "2")
    {
        iotype = PIO_IOTYPE_NETCDF;
    }
    else if (t == "netcdf4c" || t == "NETCDF4C" || t == "3")
    {
        iotype = PIO_IOTYPE_NETCDF4C;
    }
    else if (t == "netcdf4p" || t == "NETCDF4P" || t == "4")
    {
        iotype = PIO_IOTYPE_NETCDF4P;
    }
    else
    {
        throw invalid_argument("Invalid conversion type given: " + t + "\n");
    }

    return iotype;
}

int ConvertBPToNC(const string &infilepath, const string &outfilename,
                  const string &piotype, int mem_opt, MPI_Comm comm_in)
{
    int ierr = BP2PIO_NOERR, err_val = 0, err_cnt = 0;
    int ret = PIO_NOERR;
    int iosysid = 0;

    MPI_Comm comm   = comm_in;
    MPI_Comm w_comm = comm_in;
    int mpirank, w_mpirank;
    int nproc, w_nproc;

    MPI_Comm_set_errhandler(w_comm, MPI_ERRORS_RETURN);
    MPI_Comm_rank(w_comm, &w_mpirank);
    MPI_Comm_size(w_comm, &w_nproc);

    mpirank = w_mpirank;
    nproc   = w_nproc;

    /*
     * Check if the number of nodes is less than or equal to the number of BP files.
     * If not, create a new comm.
     *
     */
    int io_proc = 1;
#if 0 /* TAHSIN */
    int num_files = GetNumOfFiles(infilepath);
    int split_comm = 0;
    if (num_files < w_nproc)
    {
        if (w_mpirank == 0)
        {
            fprintf(stderr, "Warning: #files: %d < #procs: %d\n", num_files, w_nproc);
            fflush(stderr);
        }

        if (w_mpirank < num_files) /* I/O nodes */
            io_proc = 1;
		else
			io_proc = 0;

        MPI_Comm_split(w_comm, io_proc, w_mpirank, &comm);
        split_comm = 1;
        MPI_Comm_rank(comm, &mpirank);
        MPI_Comm_size(comm, &nproc);
    }
    else
    {
        io_proc = 1;
        comm    = w_comm;
        mpirank = w_mpirank;
        nproc   = w_nproc;
    }
#endif 

    TimerInitialize_nm();

    if (io_proc)
    {
        try
        {
#if 0
            iosysid = InitPIO(comm, mpirank, nproc);
            if (iosysid == BP2PIO_ERROR)
            {
                throw std::runtime_error("InitPIO error.");
            }
#endif 

            enum PIO_IOTYPE pio_iotype = GetIOType_nm(piotype);
            ierr = ConvertBPFile(infilepath, outfilename, pio_iotype, iosysid, comm, mpirank, nproc, mem_opt);
            if (ierr != BP2PIO_NOERR)
            {
                throw std::runtime_error("ConvertBPFile error.");
            }

            TimerReport_nm(comm);
        }
        catch (const std::exception &e)
        {
            fprintf(stderr, "mpirank = %d, exception: %s\n", mpirank, e.what());
            fflush(stderr);
            ierr = BP2PIO_ERROR;
        }
        catch (...)
        {
            fprintf(stderr, "mpirank = %d, exception: Unknown\n", mpirank);
            fflush(stderr);
            ierr = BP2PIO_ERROR;
        }
    }
    ERROR_CHECK_RETURN(ierr, err_val, err_cnt, comm_in)

#if 0
    if (split_comm)
        MPI_Comm_free(&comm);
#endif 

    TimerFinalize_nm();

    return BP2PIO_NOERR;
}

/* Find BP directories, named "*.bp.dir", in bppdir and the
 * corresponding file name prefixes to be used for converted
 * files
 */
static int FindBPDirs(const string &bppdir,
                      vector<string> &bpdirs,
                      vector<string> &conv_fname_prefixes)
{
    DIR *pdir = opendir(bppdir.c_str());
    if (!pdir)
    {
        fprintf(stderr, "Folder %s does not exist.\n", bppdir.c_str());
        return BP2PIO_ERROR;
    }

    const string BPDIR_NAME_RGX_STR("(.*)([.]nc)[.]bp[.]dir");
    regex bpdir_name_rgx(BPDIR_NAME_RGX_STR.c_str());
    struct dirent *pde = NULL;
    while ((pde = readdir(pdir)) != NULL)
    {
        smatch match;
        string dname(pde->d_name);
        assert(pde);
        /* Add dirs named "*.bp.dir" to bpdirs */
        if ((pde->d_type == DT_DIR) &&
            regex_search(dname, match, bpdir_name_rgx) &&
            (match.size() == 3))
        {
            conv_fname_prefixes.push_back(match.str(1));
            const std::string NC_SUFFIX(".nc");
            const std::string BP_SUFFIX(".bp");
            bpdirs.push_back(match.str(1) + match.str(2) + BP_SUFFIX);
        }
    }

    closedir(pdir);

    return BP2PIO_NOERR;
}

/* Convert all BP files in "bppdir" to NetCDF files
 * bppdir:  Directory containing multiple directories, named "*.bp.dir",
 *          each directory containing BP files corresponding to a single
 *          file. This is the "BP Parent Directory".
 * piotype: The PIO IO type used for converting BP files to NetCDF using PIO
 * comm:    The MPI communicator to be used for conversion
 *
 * The function looks for all directories in bppdir named "*.bp.dir"
 * and converts them, one at a time, to NetCDF files
 */
int MConvertBPToNC(const string &bppdir, const string &piotype, int mem_opt,
                    MPI_Comm comm)
{
    int ierr = BP2PIO_NOERR;
    vector<string> bpdirs;
    vector<string> conv_fname_prefixes;
    const std::string CONV_FNAME_SUFFIX(".nc");

    ierr = FindBPDirs(bppdir, bpdirs, conv_fname_prefixes);
    if (ierr != BP2PIO_NOERR)
    {
        fprintf(stderr, "Unable to read directory, %s\n", bppdir.c_str());
        return ierr;
    }

    assert(bpdirs.size() == conv_fname_prefixes.size());
    for (size_t i = 0; i < bpdirs.size(); i++)
    {
        MPI_Barrier(comm);
        ierr = ConvertBPToNC(bpdirs[i],
                conv_fname_prefixes[i] + CONV_FNAME_SUFFIX,
                piotype, mem_opt, comm);
        MPI_Barrier(comm);
        if (ierr != BP2PIO_NOERR)
        {
            fprintf(stderr, "Unable to convert BP file (%s) to NetCDF\n",
                    bpdirs[i].c_str());
            return ierr;
        }
    }

    return BP2PIO_NOERR;
}

#ifdef __cplusplus
extern "C" {
#endif

int C_API_ConvertBPToNC(const char *infilepath, const char *outfilename,
                        const char *piotype, int mem_opt, MPI_Comm comm_in)
{
    return ConvertBPToNC(string(infilepath), string(outfilename),
                         string(piotype), mem_opt, comm_in);
}

#ifdef __cplusplus
}
#endif
