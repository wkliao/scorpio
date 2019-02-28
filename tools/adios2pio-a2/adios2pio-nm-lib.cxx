#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <stdexcept>
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

/* debug output */
static int debug_out = 1;

nc_type PIOc_get_nctype_from_adios_type(std::string atype)
{
#define adios2_GET_TYPE(a_type,T,n_type) \
	if (a_type == adios2::GetType<T>()) { \
		return n_type; \
	}

	adios2_GET_TYPE(atype,char,NC_BYTE);
	adios2_GET_TYPE(atype,short,NC_SHORT);
	adios2_GET_TYPE(atype,int,NC_INT);
	adios2_GET_TYPE(atype,float,NC_FLOAT);
	adios2_GET_TYPE(atype,double,NC_DOUBLE);
	adios2_GET_TYPE(atype,unsigned char,NC_UBYTE);
	adios2_GET_TYPE(atype,unsigned short,NC_USHORT);
	adios2_GET_TYPE(atype,unsigned int,NC_UINT);
	adios2_GET_TYPE(atype,long,NC_INT64);
	adios2_GET_TYPE(atype,unsigned long,NC_UINT64);
	adios2_GET_TYPE(atype,std::string,NC_CHAR);
	return NC_BYTE;

#undef adios2_GET_TYPE
}

int adios2_type_size_a2(std::string atype)
{
#define adios2_GET_SIZE(a_type,T,a_size) \
	if (a_type == adios2::GetType<T>()) { \
		return a_size; \
	}

	adios2_GET_SIZE(atype,char,1);
	adios2_GET_SIZE(atype,unsigned char,1);
	adios2_GET_SIZE(atype,std::string,1);
	adios2_GET_SIZE(atype,short,2);
	adios2_GET_SIZE(atype,unsigned short,2);
	adios2_GET_SIZE(atype,int,4);
	adios2_GET_SIZE(atype,unsigned int,4);
	adios2_GET_SIZE(atype,long,8);
	adios2_GET_SIZE(atype,unsigned long,8);
	adios2_GET_SIZE(atype,float,4);
	adios2_GET_SIZE(atype,double,8);
	adios2_GET_SIZE(atype,std::complex<float>,8);
	adios2_GET_SIZE(atype,std::complex<double>,16);
	return -1;

#undef adios2_GET_SIZE
}

template <class T>
int adios2_adios_get_attr_a2(adios2::Attribute<T> *a_base, adios2::IO *bpIO, adios2::Engine *bpReader, char *aname, 
							std::string *atype, size_t *asize, void **adata) 
{
	std::string a_type = a_base->Type(); 
	*atype = a_type;
	const std::vector<T> a_data = a_base->Data();
	*asize = a_data.size(); 
	*adata = (T*)malloc(sizeof(T)*(*asize));
	memcpy(*adata,a_data.data(),a_data.size()*sizeof(T));

	return 0;
}

int adios_get_attr_a2(adios2::IO *bpIO, adios2::Engine *bpReader, char *aname, std::string *atype, size_t *asize, void **adata)
{
	std::string a_type = bpIO->AttributeType(aname);
	if (a_type.empty()) {
		return 1;
	} else if (a_type == adios2::GetType<std::string>()) {
		adios2::Attribute<std::string> a_base = bpIO->InquireAttribute<std::string>(aname);		
		const std::vector<std::string> a_data = a_base.Data();
		*asize = a_data.size();
		char **tmp_ptr = (char**)malloc(sizeof(char*)*(*asize));
        int total_size = 0;
        for (int ii=0;ii<a_data.size();ii++) {
        	tmp_ptr[ii] = (char*)malloc(sizeof(char)*(a_data[ii].length()+1));
            memcpy(tmp_ptr[ii],a_data[ii].c_str(),a_data[ii].length()+1);
		}
		*adata = (void*)tmp_ptr;
	}
#define my_declare_template_instantiation(T) \
    else if (a_type == adios2::GetType<T>()) \
    {                                        \
    	adios2::Attribute<T> a_base = bpIO->InquireAttribute<T>(aname); \
        int ret_val = adios2_adios_get_attr_a2(&a_base,bpIO, bpReader, aname, atype, asize, adata); \
    } 

    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(my_declare_template_instantiation)
#undef my_declare_template_instantiation
	
	return 0;
}

/* Timer functions to separate time for read (ADIOS) and write (PIO) */
#ifdef ADIOS_TIMING
static double time_read, time_write;
static double time_temp_read, time_temp_write;
void TimerInitialize_nm() { time_read = 0.0; time_write = 0.0; }
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
        cout << "ADIOS read time   = " << std::fixed << std::setw(8) << tr_max << "s " << std::setw(8) << tr_sum << "s\n";
        cout << "PIO  write time   = " << std::fixed << std::setw(8) << tw_max << "s " << std::setw(8) << tw_sum << "s\n";
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

struct Dimension {
    int dimid;
    PIO_Offset dimvalue;
};

using DimensionMap = std::map<std::string,Dimension>;

struct Variable {
    int nc_varid;
    bool is_timed;
    nc_type nctype;
};

using VariableMap = std::map<std::string,Variable>;

struct Decomposition {
    int ioid;
    int piotype;
};

using DecompositionMap = std::map<std::string,Decomposition>;

int InitPIO(MPI_Comm comm, int mpirank, int nproc)
{
	int iosysid; 

    if (PIOc_Init_Intracomm(comm, nproc, 1, 0, PIO_REARR_SUBSET, &iosysid))
        throw std::runtime_error("PIO initialization failed\n");

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
int GlobalMaxSteps_nm(int nsteps_in=0)
{
    static int nsteps_current = 1;
    if (nsteps_in > nsteps_current)
        nsteps_current = nsteps_in;
    return nsteps_current;
}

std::vector<int> AssignWriteRanks(int n_bp_writers, MPI_Comm comm, int mpirank, int nproc)
{
    if (!mpirank && debug_out) cout << "The BP file was written by " << n_bp_writers << " processes\n";
    int nwb = n_bp_writers / nproc; // number of blocks to process
    int start_wb; // starting wb (in [0..nwb-1])
    if (mpirank < n_bp_writers % nproc)
        nwb++;

    if (mpirank < n_bp_writers % nproc)
        start_wb = mpirank * nwb;
    else
        start_wb = mpirank * nwb + n_bp_writers % nproc;
	if (debug_out)
    	cout << "Process " << mpirank << " start block = " << start_wb <<
            	" number of blocks = " << nwb << std::endl;
    FlushStdout_nm(comm);
    std::vector<int> blocks(nwb);
    for (int i=0; i<nwb; ++i)
        blocks[i]=start_wb+i;
    return blocks;
}


void ProcessGlobalFillmode(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid)
{
    size_t  asize;
    void    *fillmode = NULL;

    std::string atype;
    adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)"/__pio__/fillmode", &atype, &asize, &fillmode);
	printf("Fillmode: %d\n",*(int*)fillmode);
    PIOc_set_fill(ncid, *(int*)fillmode, NULL);
    if (fillmode) free(fillmode);
}

void ProcessVarAttributes(adios2::IO *bpIO, adios2::Engine *bpReader, int adios_varid, std::string varname, int ncid, int nc_varid)
{
	std::map<std::string,adios2::Params> a2_vi = bpIO[0].AvailableAttributes(varname);

    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter!=a2_vi.end(); a2_iter++) {
		if (debug_out) cout << "    Attribute: " << a2_iter->first << std::endl;
		if (a2_iter->first.find("__pio__/")==string::npos) {
			std::string aa_name = varname + "/" + a2_iter->first;
			std::string atype = bpIO[0].AttributeType(aa_name.c_str());
			size_t asize;
			int len;
			if (atype == adios2::GetType<std::string>()) {
        		char   **adata = NULL;
        		adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)aa_name.c_str(), &atype, &asize, (void**)&adata);
        		nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
        		char *attname = (char*) a2_iter->first.c_str();
        		len = strlen(*adata);
				if (debug_out) cout << "        define PIO attribute: " << attname << ""
           		  		 			<< "  type=" << piotype << std::endl;
        		PIOc_put_att(ncid, nc_varid, attname, piotype, len, *adata);
        		if (adata) free(adata);	
			} else {
        		char   *adata = NULL;
        		adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)aa_name.c_str(), &atype, &asize, (void**)&adata);
        		nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
        		char *attname = (char*) a2_iter->first.c_str();
				if (debug_out) cout << "        define PIO attribute: " << attname << ""
           			  		 		<< "  type=" << piotype << std::endl;
        		len = 1;
        		PIOc_put_att(ncid, nc_varid, attname, piotype, len, adata);
        		if (adata) free(adata);	
			}
		}
    }
}

void ProcessGlobalAttributes(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, DimensionMap& dimension_map, VariableMap& vars_map)
{
	if (debug_out) cout << "Process Global Attributes: \n";

	std::string delimiter = "/";
	std::map<std::string,char> processed_attrs;
	std::map<std::string,int>  var_att_map;

	std::map<std::string,adios2::Params> a2_attr = bpIO[0].AvailableAttributes();
	size_t nattrs = a2_attr.size();

	std::map<std::string,adios2::Params>::iterator a2_iter = a2_attr.begin();

	int total_cnt = (int) nattrs, i = -1;
	char attr_name[512];
	size_t aname_size;
	while (total_cnt>0) {
		std::string a = a2_iter->first; 
		a2_iter++;
		if (a2_iter==a2_attr.end()) 
			a2_iter = a2_attr.begin();
		char *attr_namelist = (char*)a.c_str(); 

		std::string token = a.substr(0, a.find(delimiter));
		if (token=="") {  /* Not a variable with attributes */
			processed_attrs[a] = 1; total_cnt--;
		} else if (processed_attrs.find(a)==processed_attrs.end()) {
        	if (a.find("pio_global/") != string::npos) {
            	if (debug_out) cout << " GLOBAL Attribute: " << attr_namelist << std::endl;
            	size_t asize;
            	char *adata = NULL;
            	std::string atype;
            	adios_get_attr_a2(&bpIO[0], &bpReader[0], attr_namelist, &atype, &asize, (void**)&adata);
            	nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
            	char *attname = attr_namelist+strlen("pio_global/");
            	if (debug_out) cout << "        define PIO attribute: " << attname << ""
               		     			<< "  type=" << piotype << std::endl;
            	int len = 1;
            	if (atype == adios2::GetType<std::string>())
                	len = strlen(adata);
            	PIOc_put_att(ncid, PIO_GLOBAL, attname, piotype, len, adata);
            	if (adata) free(adata);
				processed_attrs[a] = 1; total_cnt--;
			} else {
				if (debug_out) cout << "    Attribute: " << attr_namelist << std::endl;
				if (vars_map.find(token)==vars_map.end()) { 
					if (var_att_map.find(token)==var_att_map.end()) {
						// first encounter 
           				string attname = token + "/__pio__/nctype";
						processed_attrs[attname] = 1; total_cnt--;
           				size_t asize;
           				int *nctype = NULL;
           				std::string atype;
           				adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&nctype);

						attname = token + "/__pio__/ndims";
						processed_attrs[attname] = 1; total_cnt--;
      		      		int *ndims = NULL;
      		      		adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&ndims);
	
       		     		char **dimnames = NULL;
       		     		int dimids[PIO_MAX_DIMS];
       		     		if (*ndims) {
							attname = token + "/__pio__/dims";
							processed_attrs[attname] = 1; total_cnt--;
							adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&dimnames);

							for (int d=0; d < *ndims; d++)
								dimids[d] = dimension_map[dimnames[d]].dimid;
						}
       		     		int varid;
       		     		PIOc_def_var(ncid, token.c_str(), *nctype, *ndims, dimids, &varid);
	   	        		var_att_map[token] = varid; 

       		     		if (nctype) free(nctype);
       		     		if (ndims) free(ndims);
       		     		if (dimnames) free(dimnames);
					} else {
						if (processed_attrs.find(a)==processed_attrs.end()) {
							processed_attrs[a] = 1; total_cnt--;
							std::string atype = bpIO[0].AttributeType(a.c_str());
							size_t asize;
							int len;
							if (atype == adios2::GetType<std::string>()) {
       	    					char **adata = NULL;
       	    					adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)a.c_str(), &atype, &asize, (void**)&adata);
								nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
       		 					char *attname = ((char*)a.c_str())+token.length()+1;;
       		 					len = strlen(*adata);
       		 					PIOc_put_att(ncid, var_att_map[token], attname, piotype, len, *adata);
       		 					if (adata) free(adata);
							} else {
       	    					char *adata = NULL;
       	    					adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)a.c_str(), &atype, &asize, (void**)&adata);
								nc_type piotype = PIOc_get_nctype_from_adios_type(atype);
       		 					char *attname = ((char*)a.c_str())+token.length()+1;;
       		 					len = 1;
       		 					PIOc_put_att(ncid, var_att_map[token], attname, piotype, len, adata);
       		 					if (adata) free(adata);
							}
						}
					}
				}
			}
		}
	}
}

template <class T>
Decomposition adios2_ProcessOneDecomposition(adios2::Variable<T> *v_base,  
		adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, 
        const char *varname, std::vector<int>& wfiles, int iosysid, int mpirank,
        int nproc, int forced_type=NC_NAT)
{
    /* 
 	 * Read all decomposition blocks assigned to this process,
     * create one big array from them and create a single big 
     * decomposition with PIO
     */

    /* Sum the sizes of blocks assigned to this process */
    TimerStart(read);
    uint64_t nelems = 0;
    for (int i=1;i<=wfiles.size();i++) { // iterate over all the files assigned to this process
		std::string v_type = bpIO[i].VariableType(varname);
       	*v_base = bpIO[i].InquireVariable<T>(varname);  
		const auto v_blocks = bpReader[i].BlocksInfo(*v_base,(const size_t)0);
		for (int j=0;j<v_blocks.size();j++) 
			nelems += v_blocks[j].Count[0];
	}

	/* allocate +1 to prevent d.data() from returning NULL. Otherwise, read/write operations fail */
	/* nelems may be 0, when some processes do not have any data */
	std::string d_type = bpIO[1].VariableType(varname);	
	char *d = NULL;
	int d_size;
	if (d_type == adios2::GetType<int>()) {
		d = (char*)malloc(sizeof(int)*(nelems+1));
		d_size = sizeof(int);
	} else if (d_type == adios2::GetType<long int>()) {
		d = (char*)malloc(sizeof(long int)*(nelems+1));
		d_size = sizeof(long int);
	} else {
		printf("ERROR: unsupported type: %s\n",d_type.c_str());
		return Decomposition{-1, -1};
	}

   	uint64_t offset = 0;
   	for (int i=1;i<=wfiles.size();i++) {
		std::string v_type = bpIO[i].VariableType(varname);
       	*v_base = bpIO[i].InquireVariable<T>(varname);  
		std::vector<T> v_data;
		const auto v_blocks = bpReader[i].BlocksInfo(*v_base,0); 
   		for (int j=0;j<v_blocks.size();j++) { 
   	 		if (debug_out) cout << " rank " << mpirank << ": read decomp wb = " << j << 
   	   			  				" start = " << offset << 
   	   			  				" elems = " << v_blocks[j].Count[0] << endl; 
			v_base->SetBlockSelection(j); 
			adios2::Dims start = {0}, count = {v_blocks[j].Count[0]};
			v_base->SetSelection({start,count});
			bpReader[i].Get(*v_base,v_data,adios2::Mode::Sync); 
			memcpy((d+offset),v_data.data(),v_data.size()*d_size); 
			offset += (v_blocks[j].Count[0]*d_size); 
		} 
   	}

   	std::string attname;
   	size_t asize;
   	int *piotype = NULL;
   	std::string atype;
   	if (forced_type == NC_NAT) {
   		 attname = string(varname) + "/piotype";
   		 adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&piotype);
   	} else {
   		 piotype = (int*) malloc(sizeof(int));
   		 *piotype = forced_type;
   	}

   	attname = string(varname) + "/ndims";
   	int *decomp_ndims = NULL;
   	adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&decomp_ndims);

   	int *decomp_dims = NULL;
   	attname = string(varname) + "/dimlen";
   	adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&decomp_dims);
   	TimerStop(read);

   	TimerStart(write);
   	int ioid;
   	PIOc_InitDecomp(iosysid, *piotype, *decomp_ndims, decomp_dims, (PIO_Offset)nelems, (PIO_Offset*)d, &ioid, NULL, NULL, NULL);
   	TimerStop(write);

	int decomp_piotype = *piotype;

   	if (piotype) free(piotype);
   	if (decomp_ndims) free(decomp_ndims);
   	if (decomp_dims) free(decomp_dims);
	if (d) free(d);

    return Decomposition{ioid, decomp_piotype};
}

Decomposition ProcessOneDecomposition(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, 
		const char *varname, std::vector<int>& wfiles, int iosysid, int mpirank,
		int nproc, int forced_type=NC_NAT)
{
    TimerStart(read);

	std::string v_type = bpIO[0].VariableType(varname);

	if (v_type.empty()) {
		return Decomposition{-1, -1};
	}
#define declare_template_instantiation(T)      \
   	else if (v_type == adios2::GetType<T>())   \
   	{                                          \
       	adios2::Variable<T> v_base;    \
		return adios2_ProcessOneDecomposition(&v_base, bpIO, bpReader, ncid, \
       			varname, wfiles, iosysid, mpirank, nproc, forced_type); \
   	}

   	ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation	
}


DecompositionMap ProcessDecompositions(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, 
							std::vector<int>& wfiles, int iosysid, 
							MPI_Comm comm, int mpirank, int nproc)
{
    DecompositionMap decomp_map;
	std::map<std::string, adios2::Params> a2_vi = bpIO[0].AvailableVariables();

    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter!=a2_vi.end(); a2_iter++)
    {
        string v = a2_iter->first;
        if (v.find("/__pio__/decomp/") != string::npos)
        {
            string decompname = v.substr(16);
            if (!mpirank && debug_out) cout << "Process decomposition " << decompname << endl;
            Decomposition d = ProcessOneDecomposition(bpIO, bpReader, ncid, (char*)a2_iter->first.c_str(), wfiles, iosysid, mpirank, nproc);
            decomp_map[decompname] = d;
        }
        FlushStdout_nm(comm);
    }

    return decomp_map;
}

Decomposition GetNewDecomposition(DecompositionMap& decompmap, string decompname,
        adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, std::vector<int>& wfiles, int nctype, int iosysid, int mpirank, int nproc)
{
	char ss[256];
	sprintf(ss,"%s_%d",decompname.c_str(),nctype);
	string key(ss);

    auto it = decompmap.find(key);
    Decomposition d;
    if (it == decompmap.end())
    {
        string varname = "/__pio__/decomp/" + decompname;
        d = ProcessOneDecomposition(bpIO, bpReader, ncid, (char*)varname.c_str(), wfiles, iosysid, mpirank, nproc, nctype);
        decompmap[key] = d;
    }
    else
    {
        d = it->second;
    }
    return d;
}

DimensionMap ProcessDimensions(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, MPI_Comm comm, int mpirank, int nproc)
{
    DimensionMap dimensions_map;

	std::map<std::string, adios2::Params> a2_vi = bpIO[0].AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter!=a2_vi.end(); a2_iter++)
    {
        string v = a2_iter->first;
		const char *varname = v.c_str();
        if (v.find("/__pio__/dim/") != string::npos)
        {
            /* For each dimension stored, define a dimension variable with PIO */
            string dimname = v.substr(13);
            if (!mpirank && debug_out) cout << "Process dimension " << dimname << endl;

			std::string v_type = bpIO[0].VariableType(a2_iter->first);
            TimerStart(read);
			if (v_type.empty()) {
				return dimensions_map;
			}
#define declare_template_instantiation(T)              \
    		else if (v_type == adios2::GetType<T>())   \
    		{                                                  \
        		adios2::Variable<T> v_base = bpIO[0].InquireVariable<T>(varname);       \
				std::vector<T> dimval; \
				bpReader[0].Get(v_base,dimval,adios2::Mode::Sync); \
				int dimid; \
            	TimerStart(write); \
				PIO_Offset *d_val = (PIO_Offset*)dimval.data(); \
            	PIOc_def_dim(ncid, dimname.c_str(), *d_val, &dimid); \
            	TimerStop(write); \
            	dimensions_map[dimname] = Dimension{dimid,(PIO_Offset)*d_val}; \
			} 

    		ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation	
		}
        FlushStdout_nm(comm);
    }
    return dimensions_map;
}

VariableMap ProcessVariableDefinitions(adios2::IO *bpIO, adios2::Engine *bpReader, int ncid, DimensionMap& dimension_map, MPI_Comm comm, int mpirank, int nproc)
{
    VariableMap vars_map;
	std::map<std::string, adios2::Params> a2_vi = bpIO[0].AvailableVariables();
    for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter!=a2_vi.end(); a2_iter++)
    {
        string v = a2_iter->first;
        if (!mpirank && debug_out) cout << "BEFORE Process variable " << v << endl;
        if (v.find("/__") == string::npos)
        {
            /* For each variable written define it with PIO */
            if (!mpirank && debug_out) cout << "Process variable " << v << endl;

			if (v.find("decomp_id/")==string::npos && 
				v.find("frame_id/")==string::npos && 
				v.find("fillval_id/")==string::npos) {

            	TimerStart(read);
            	string attname = v + "/__pio__/nctype";
            	size_t asize;
            	int *nctype = NULL;
            	std::string atype;
            	adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&nctype);

	            attname = v + "/__pio__/ndims";
	            int *ndims = NULL;
	            adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&ndims);

	            char **dimnames = NULL;
	            int  dimids[PIO_MAX_DIMS];
	            bool timed = false;
	            if (*ndims)
	            {
	                attname = v + "/__pio__/dims";
   					adios_get_attr_a2(&bpIO[0], &bpReader[0], (char*)attname.c_str(), &atype, &asize, (void**)&dimnames);

	                for (int d=0; d < *ndims; d++)
	                {
	                    dimids[d] = dimension_map[dimnames[d]].dimid;
	                    if (dimension_map[dimnames[d]].dimvalue == PIO_UNLIMITED) 
	                    {
							timed = true;
	                    }
	                }
	            }
	            TimerStop(read);

	            TimerStart(write);
	            int varid;
	            PIOc_def_var(ncid, (char*)v.c_str(), *nctype, *ndims, dimids, &varid);
	            TimerStop(write);
	            vars_map[v] = Variable{varid,timed,*nctype};

	            ProcessVarAttributes(bpIO,bpReader, 0, v, ncid, varid);

	            if (nctype) free(nctype);
	            if (ndims) free(ndims);
	            if (dimnames) free(dimnames);
			}
        }
        FlushStdout_nm(comm);
    }
    return vars_map;
}

int put_var_nm(int ncid, int varid, int nctype, std::string memtype, const void* buf)
{
    int ret = 0;

    if (memtype == adios2::GetType<char>()) { 
		if (nctype==PIO_CHAR)
			ret = PIOc_put_var_text(ncid, varid, (const char*)buf);
		else  
			ret = PIOc_put_var_schar(ncid, varid, (const signed char*)buf);
	} else if (memtype == adios2::GetType<short int>()) {
        ret = PIOc_put_var_short(ncid, varid, (const signed short*)buf);
	} else if (memtype == adios2::GetType<int>()) {
        ret = PIOc_put_var_int(ncid, varid, (const signed int*)buf);
	} else if (memtype == adios2::GetType<float>()) {
        ret = PIOc_put_var_float(ncid, varid, (const float *)buf);
	} else if (memtype == adios2::GetType<double>()) {
        ret = PIOc_put_var_double(ncid, varid, (const double *)buf);
	} else if (memtype == adios2::GetType<unsigned char>()) {
        ret = PIOc_put_var_uchar(ncid, varid, (const unsigned char *)buf);
	} else if (memtype == adios2::GetType<unsigned short>()) {
        ret = PIOc_put_var_ushort(ncid, varid, (const unsigned short *)buf);
	} else if (memtype == adios2::GetType<unsigned int>()) {
        ret = PIOc_put_var_uint(ncid, varid, (const unsigned int *)buf);
	} else if (memtype == adios2::GetType<long long int>()) {
        ret = PIOc_put_var_longlong(ncid, varid, (const signed long long *)buf);
	} else if (memtype == adios2::GetType<unsigned long long int>())  {
        ret = PIOc_put_var_ulonglong(ncid, varid, (const unsigned long long *)buf);
	} else if (memtype == adios2::GetType<std::string>()) {
        ret = PIOc_put_var_text(ncid, varid, (const char *)buf);
	} else {
        /* We can't do anything here, hope for the best, i.e. memtype equals to nctype */
        ret = PIOc_put_var(ncid, varid, buf);
    }
    return ret;
}

int put_vara_nm(int ncid, int varid, int nctype, std::string memtype, PIO_Offset *start, PIO_Offset *count, const void* buf)
{
    int ret = 0;
	if (memtype == adios2::GetType<char>()) {
		if (nctype == PIO_BYTE) 
            ret = PIOc_put_vara_schar(ncid, varid, start, count, (const signed char*)buf);
        else 
            ret = PIOc_put_vara_text(ncid, varid, start, count, (const char*)buf);
	} else if (memtype == adios2::GetType<short>()) {
		ret = PIOc_put_vara_short(ncid, varid, start, count, (const signed short*)buf);
	} else if (memtype == adios2::GetType<int>()) {
		ret = PIOc_put_vara_int(ncid, varid, start, count, (const signed int*)buf);
	} else if (memtype == adios2::GetType<float>()) {
		ret = PIOc_put_vara_float(ncid, varid, start, count, (const float *)buf);
	} else if (memtype == adios2::GetType<double>()) {
		ret = PIOc_put_vara_double(ncid, varid, start, count, (const double *)buf);
	} else if (memtype == adios2::GetType<unsigned char>()) {
		ret = PIOc_put_vara_uchar(ncid, varid, start, count, (const unsigned char *)buf);
	} else if (memtype == adios2::GetType<unsigned short>()) {
		ret = PIOc_put_vara_ushort(ncid, varid, start, count, (const unsigned short *)buf);
	} else if (memtype == adios2::GetType<unsigned int>()) {
		ret = PIOc_put_vara_uint(ncid, varid, start, count, (const unsigned int *)buf);
	} else if (memtype == adios2::GetType<long long int>()) {
		ret = PIOc_put_vara_longlong(ncid, varid, start, count, (const signed long long *)buf);
	} else if (memtype == adios2::GetType<unsigned long long int>()) {
		ret = PIOc_put_vara_ulonglong(ncid, varid, start, count, (const unsigned long long *)buf);
	} else if (memtype == adios2::GetType<std::string>()) {
		ret = PIOc_put_vara_text(ncid, varid, start, count, (const char *)buf);
	} else {
		ret = PIOc_put_vara(ncid, varid, start, count, buf);
	}
    return ret;
}

template <class T>
int adios2_ConvertVariablePutVar(adios2::Variable<T> *v_base, 
								adios2::IO *bpIO, adios2::Engine *bpReader, std::vector<int> wfiles, 
								std::string varname, int ncid, Variable& var, 
								int mpirank, int nproc)
{
    TimerStart(read);
    int ret = 0;

	*v_base = bpIO[0].InquireVariable<T>(varname);

	adios2::Dims v_dims = v_base->Shape(0); 
	if (v_dims.size()==0) {	
		/* Scalar variable */
        TimerStart(write);
		std::vector<T> v_value;
		bpReader[0].Get(*v_base,v_value,adios2::Mode::Sync);
       	ret = put_var_nm(ncid, var.nc_varid, var.nctype, v_base->Type(), v_value.data()); 
       	if (ret != PIO_NOERR) { 
			cout << "ERROR in PIOc_put_var(), code = " << ret 
				<< " at " << __func__ << ":" << __LINE__ << endl; 
		}
        TimerStop(write); 
	} else {
        /* An N-dimensional array that needs no rearrangement.
         * put_vara_nm() needs all processes participate */
        TimerStart(read);
		
		/* ACME writes this array from I/O processor 0 */
        PIO_Offset start[PIO_MAX_DIMS], count[PIO_MAX_DIMS];
		// PIOc_put_var may have been called multiple times with different start,count values 
		// for a variable. We need to convert the output from each of those calls.
	
		if (v_base) {
			const auto v_blocks = bpReader[0].BlocksInfo(*v_base,(const size_t)0);
			if (mpirank==0) {
				size_t mysize = 1, mycount;
				char   *buf   = NULL;
				for (int ii=0;ii<v_blocks.size();ii++) {
					mysize = 1;
					for (int d=0;d<v_dims.size();d++) 
						mysize *= (size_t) v_blocks[ii].Count[d];
					mycount = mysize;
					mysize = mysize*adios2_type_size_a2(v_base->Type());
   	   				if ((buf = (char *)malloc(mysize))==NULL) {
						printf("ERROR: cannot allocate memory: %ld\n",mysize);
						return 1;
					}
					v_base->SetBlockSelection(ii);
					v_base->SetSelection({v_blocks[ii].Start,v_blocks[ii].Count});
					std::vector<T> v_value;
					bpReader[0].Get(*v_base,v_value,adios2::Mode::Sync);
					memcpy(buf,v_value.data(),mysize);
					for (int d=0;d<v_dims.size();d++) {
						start[d] = (PIO_Offset) v_blocks[ii].Start[d]; 
       					count[d] = (PIO_Offset) v_blocks[ii].Count[d]; 
					}
 					ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base->Type(), start, count, buf);
       				if (ret != PIO_NOERR) {
       					cout << "rank " << mpirank << ":ERROR in PIOc_put_vara(), code = " << ret
       						 << " at " << __func__ << ":" << __LINE__ << endl;
						return 1;
					}
					if (buf) free(buf);
				}
			} else {
				char   temp_buf;
				for (int ii=0;ii<v_blocks.size();ii++) {
					for (int d=0;d<v_dims.size();d++) {
						start[d] = (PIO_Offset) 0;
       					count[d] = (PIO_Offset) 0;
					}
 					ret = put_vara_nm(ncid, var.nc_varid, var.nctype, v_base->Type(), start, count, &temp_buf);
       				if (ret != PIO_NOERR) {
       					cout << "rank " << mpirank << ":ERROR in PIOc_put_vara(), code = " << ret
       						 << " at " << __func__ << ":" << __LINE__ << endl;
						return 1;
					}
				}
			}
		}

        TimerStop(write);
    }

    return ret;
}

int ConvertVariablePutVar(adios2::IO *bpIO, adios2::Engine *bpReader, std::vector<int> wfiles, 
						std::string varname, int ncid, Variable& var, int mpirank, int nproc)
{
    TimerStart(read);
    int ret = 0;

	std::string v_type = bpIO->VariableType(varname);
	if (v_type.empty()) {
		return 1;
	}
#define declare_template_instantiation(T)              \
    else if (v_type == adios2::GetType<T>())   \
    {                                                  \
        adios2::Variable<T> v_base; \
		std::vector<T> v_value; \
		return adios2_ConvertVariablePutVar(&v_base, bpIO, bpReader, wfiles, \
								varname, ncid, var, mpirank, nproc); \
    }
    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation	
}

template <class T>
int adios2_ConvertVariableTimedPutVar(adios2::Variable<T> *v_base, std::vector<T> v_value,
                                adios2::IO *bpIO, adios2::Engine *bpReader, std::vector<int> wfiles, 
								std::string varname, int ncid, Variable& var, int nblocks_per_step,
                              	MPI_Comm comm, int mpirank, int nproc)
{
    TimerStart(read);
    int ret = 0;

	adios2::Dims v_dims = v_base->Shape(0); 
    if (v_dims.size()== 0)
    {
        /* Scalar variable over time */
        /* Written by only one process, so steps = number of blocks in file */
		const auto v_blocks = bpReader[0].BlocksInfo(*v_base,0);

        TimerStart(read);
        int nsteps = v_blocks.size(); 
		std::vector<T> v_mins(nsteps);
		for (int ts=0; ts<nsteps; ts++)
			v_mins[ts] = v_base->Min(ts);
        TimerStart(read);

        PIO_Offset start[1], count[1];
        for (int ts=0; ts < nsteps; ++ts)
        {
            TimerStart(write);
            start[0] = ts;
            count[0] = 1;
            int ret = PIOc_put_vara(ncid, var.nc_varid, start, count, v_mins.data()); 
            if (ret != PIO_NOERR)
				cout << "ERROR in PIOc_put_var(), code = " << ret
					<< " at " << __func__ << ":" << __LINE__ << endl;
            TimerStop(write);
        }
    }
    else
    {
        /* calculate how many records/steps we have for this variable */
        int nsteps = 1;

		/* compute the total number of blocks */
		int l_nblocks = 0;
		int g_nblocks = 0;
		for (int i=1;i<=wfiles.size();i++) {
			const auto v_blocks = bpReader[i].BlocksInfo(*v_base,0);
			l_nblocks += v_blocks.size();
		}
		MPI_Allreduce(&l_nblocks,&g_nblocks,1,MPI_INT,MPI_SUM,comm); 
		
        if (var.is_timed)
        {
            nsteps = g_nblocks / nblocks_per_step;
        }
        if (g_nblocks != nsteps * nblocks_per_step)
        {
            if (debug_out) cout << "rank " << mpirank << ":ERROR in processing variable '" << varname 
                 << "'. Number of blocks = " << g_nblocks
                 << " does not equal the number of steps * number of writers = "
                 << nsteps << " * " << nblocks_per_step << " = " << nsteps*nblocks_per_step
                 << endl;
        }

        /* Is this a local array written by each process, or a truly distributed global array */
        TimerStart(read);
		const auto v_blocks = bpReader[0].BlocksInfo(*v_base,0);
        TimerStart(read);
        bool local_array = true;
		typename adios2::Variable<T>::Info v_info = v_blocks[0];
		adios2::Dims b_dims = v_info.Count;
        for (int d = 0; d < v_dims.size(); d++)
        {
            if (b_dims[d] != v_dims[d])
            {
                local_array = false;
                break;
            }
        }
        if (var.nctype == PIO_CHAR && v_dims.size() == 1)
        {
            /* Character array over time may have longer dimension declaration than the actual content */
            local_array = true;
        }

        if (local_array)
        {
            /* Just read the arrays written by rank 0 (on every process here) and
             * write it collectively.
             */
            for (int ts=0; ts < nsteps; ++ts)
            {
                TimerStart(read);
                int elemsize = adios2_type_size_a2(v_base->Type());
                uint64_t nelems = 1;
                for (int d = 0; d < v_dims.size(); d++)
                {
                    nelems *= v_dims[d]; 
                }
                std::vector<char> d(nelems * elemsize);
				v_base->SetBlockSelection(ts);
				adios2::Dims vd_start(v_dims.size());
				adios2::Dims vd_count(v_dims.size());
				for (int d = 0; d < v_dims.size(); d++) {
					vd_start[d] = 0;
					vd_count[d] = v_dims[d];
				}
				v_base->SetSelection({vd_start,vd_count});
				std::vector<T> v_data; 
				bpReader[0].Get(*v_base,v_data,adios2::Mode::Sync);
				memcpy(d.data(),v_data.data(),nelems * elemsize);
                TimerStop(read);

                TimerStart(write);
                PIO_Offset start[v_dims.size()+1], count[v_dims.size()+1];
                start[0] = ts;
                count[0] = 1;
                for (int d = 0; d < v_dims.size(); d++)
                {
                    start[d+1] = 0;
                    count[d+1] = v_dims[d];
                }
                if ((ret = PIOc_put_vara(ncid, var.nc_varid, start, count, d.data())))
                if (ret != PIO_NOERR)
                    cout << "ERROR in PIOc_put_var(), code = " << ret
                    << " at " << __func__ << ":" << __LINE__ << endl;
                TimerStop(write);
            }
        }
        else
        {
			/* PIOc_put_vara_ writes out from processor 0       */
			/* Read in infile[0] and output using PIOc_put_vara */
			if (mpirank==0) {
				for (int ts=0; ts < nsteps; ++ts)
           		{
               		TimerStart(read);
               		int elemsize = adios2_type_size_a2(v_base->Type());
               		uint64_t nelems = 1;
               		for (int d = 0; d < v_dims.size(); d++) {
               	    	nelems *= v_dims[d];
               		}
               		std::vector<char> d(nelems * elemsize);
					v_base->SetBlockSelection(ts);
					adios2::Dims vd_start(v_dims.size());
					adios2::Dims vd_count(v_dims.size());
					for (int d = 0; d < v_dims.size(); d++) {
						vd_start[d] = 0;
						vd_count[d] = v_dims[d];
					}
					v_base->SetSelection({vd_start,vd_count});
					std::vector<T> v_data; 
					bpReader[0].Get(*v_base,v_data,adios2::Mode::Sync);
					memcpy(d.data(),v_data.data(),nelems * elemsize);
               		TimerStop(read);
	
   	            	TimerStart(write);
   	            	PIO_Offset start[v_dims.size()+1], count[v_dims.size()+1];
   	            	start[0] = ts;
   	            	count[0] = 1;
   	            	for (int d = 0; d < v_dims.size(); d++) {
						start[d+1] = 0;
   	                	count[d+1] = v_dims[d];
   	            	}
	
   	            	if ((ret = PIOc_put_vara(ncid, var.nc_varid, start, count, d.data())))
   	            	if (ret != PIO_NOERR)
   	                	cout << "ERROR in PIOc_put_var(), code = " << ret
   	                	<< " at " << __func__ << ":" << __LINE__ << endl;
   	            	TimerStop(write);
   	        	}
			} else {
				for (int ts=0; ts < nsteps; ++ts)
           		{
               		TimerStart(read);
               		int elemsize = adios2_type_size_a2(v_base->Type());
               		uint64_t nelems = 1;
               		std::vector<char> d(nelems * elemsize);

   	            	PIO_Offset start[v_dims.size()+1], count[v_dims.size()+1];
   	            	start[0] = 0;
   	            	count[0] = 0;
   	            	for (int d = 0; d < v_dims.size(); d++) {
						start[d+1] = 0;
   	                	count[d+1] = 0; 
   	            	}
	
   	            	if ((ret = PIOc_put_vara(ncid, var.nc_varid, start, count, d.data())))
   	            	if (ret != PIO_NOERR)
   	                	cout << "ERROR in PIOc_put_var(), code = " << ret
   	                	<< " at " << __func__ << ":" << __LINE__ << endl;
   	        	}
			}
        }
    }
    return ret;
}

int ConvertVariableTimedPutVar(adios2::IO *bpIO, adios2::Engine *bpReader, std::vector<int> wfiles, 
							std::string varname, int ncid, Variable& var, int nblocks_per_step,
							MPI_Comm comm, int mpirank, int nproc)
{
    TimerStart(read);
    int ret = 0;

	std::string v_type = bpIO->VariableType(varname);
	if (v_type.empty()) {
		return 1;
	}
#define declare_template_instantiation(T)              \
    else if (v_type == adios2::GetType<T>())   \
    {                                                  \
        adios2::Variable<T> v_base = bpIO[0].InquireVariable<T>(varname.c_str());       \
		std::vector<T> v_value; \
		return adios2_ConvertVariableTimedPutVar(&v_base, v_value, bpIO, bpReader, wfiles, \
                                varname, ncid, var, nblocks_per_step, \
                                comm, mpirank, nproc); \
    }
    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation	
}


template <class T>
int adios2_ConvertVariableDarray(adios2::Variable<T> *v_base, std::vector<T> v_value,
                            adios2::IO *bpIO, adios2::Engine *bpReader, 
							std::string varname, int ncid, Variable& var,
        					std::vector<int>& wfiles, DecompositionMap& decomp_map, int nblocks_per_step, int iosysid,
							MPI_Comm comm, int mpirank, int nproc, int mem_opt)
{
    int ret = 0;

    /* calculate how many records/steps we have for this variable */
    int nsteps = 1;
    int ts = 0; /* loop from ts to nsteps, below ts may become the last step */

	/* compute the total number of blocks */
	int l_nblocks = 0;
	int g_nblocks = 0;
	for (int i=1;i<=wfiles.size();i++) {
		*v_base = bpIO[i].InquireVariable<T>(varname);
		const auto vb_blocks = bpReader[i].BlocksInfo(*v_base,0);
		l_nblocks += vb_blocks.size();
	}
	MPI_Allreduce(&l_nblocks,&g_nblocks,1,MPI_INT,MPI_SUM,comm); 

	nsteps = g_nblocks / nblocks_per_step;
	if (g_nblocks != nsteps * nblocks_per_step) 
	{
		if (debug_out) cout << "rank " << mpirank << ":ERROR in processing darray '"  
				<< "'. Number of blocks = " << g_nblocks 
				<< " does not equal the number of steps * number of writers = "
				<< nsteps << " * " << nblocks_per_step << " = " << nsteps*nblocks_per_step
				<< endl;
	}

	/* different decompositions at different frames */
	char decomp_varname[128];
	char frame_varname[128];
	char fillval_varname[128];
	char decompname[64];
	sprintf(decomp_varname,"decomp_id/%s",varname.c_str());
	sprintf(frame_varname,"frame_id/%s",varname.c_str());
	sprintf(fillval_varname,"fillval_id/%s",varname.c_str());
	int  decomp_id, frame_id, fillval_exist; 
	char fillval_id[16];

	// TAHSIN -- THIS IS GETTING CONFUSING. NEED TO THINK ABOUT time steps. 
    for (; ts < nsteps; ++ts)
    {
        /* Sum the sizes of blocks assigned to this process */
		/* Compute the number of writers for each file from nsteps */
        uint64_t nelems = 0;
		int l_nwriters  = 0;
		for (int i=1;i<=wfiles.size();i++) {
			*v_base = bpIO[i].InquireVariable<T>(varname);
			const auto vb_blocks = bpReader[i].BlocksInfo(*v_base,0);
			l_nwriters = vb_blocks.size()/nsteps;
			for (int j=0;j<l_nwriters;j++) {	
				int blockid = j*nsteps+ts;
				if (blockid<vb_blocks.size())
					nelems += vb_blocks[blockid].Count[0];
			}
    	}

		/* Read local data for each file */
        int elemsize = adios2_type_size_a2(v_base->Type());
		/* allocate +1 to prevent d.data() from returning NULL. Otherwise, read/write operations fail */
		/* nelems may be 0, when some processes do not have any data */
        std::vector<char> d((nelems+1) * elemsize);
        uint64_t offset = 0;
		for (int i=1;i<=wfiles.size();i++) {
			*v_base = bpIO[i].InquireVariable<T>(varname);
			const auto vb_blocks = bpReader[i].BlocksInfo(*v_base,0);
            l_nwriters = vb_blocks.size()/nsteps;
            for (int j=0;j<l_nwriters;j++) {
				int blockid = j*nsteps+ts;
                if (blockid<vb_blocks.size()) {
					v_base->SetBlockSelection(blockid);
					v_base->SetSelection({vb_blocks[blockid].Start,vb_blocks[blockid].Count});
					bpReader[i].Get(*v_base,v_value,adios2::Mode::Sync);
					memcpy(d.data()+offset,v_value.data(),v_value.size()*elemsize);

					std::vector<int> v_tmp;
					adios2::Variable<int> v_var;
					v_var = bpIO[i].InquireVariable<int>(decomp_varname);
					v_var.SetBlockSelection(blockid);
					bpReader[i].Get(v_var,v_tmp,adios2::Mode::Sync);
					memcpy(&decomp_id,v_tmp.data(),sizeof(int));

					v_var = bpIO[i].InquireVariable<int>(frame_varname);
					v_var.SetBlockSelection(blockid);
					bpReader[i].Get(v_var,v_tmp,adios2::Mode::Sync);
					memcpy(&frame_id,v_tmp.data(),sizeof(int));

					/* Fix for NUM_FRAMES */
					if (!var.is_timed && frame_id>=0)
						var.is_timed = true;
						
					if (decomp_id>0) {
						adios2::Variable<T> v1_var = bpIO[i].InquireVariable<T>(fillval_varname);
						std::vector<T> v1_tmp;
						v1_var.SetBlockSelection(blockid);
                       	bpReader[i].Get(v1_var,v1_tmp,adios2::Mode::Sync);
						memcpy(fillval_id,v1_tmp.data(),v1_tmp.size()*sizeof(T));
						fillval_exist = 1;
					} else {
						decomp_id = -decomp_id;
						fillval_exist = 0;
					}
               		offset += vb_blocks[blockid].Count[0]*elemsize;
				}
            }   
        }	
        TimerStop(read);

        TimerStart(write);
		Decomposition decomp;
		if (mem_opt) {
			sprintf(decompname,"/__pio__/decomp/%d",decomp_id);
        	decomp = ProcessOneDecomposition(bpIO, bpReader, ncid, decompname, wfiles, iosysid, mpirank, nproc);
		} else {
			sprintf(decompname,"%d",decomp_id);
			decomp = decomp_map[decompname];
		} 
    	if (decomp.piotype != var.nctype) {
       		/* Type conversion may happened at writing. Now we make a new decomposition for this nctype */
			if (mem_opt) {
				PIOc_freedecomp(iosysid,decomp.ioid);
        		decomp = ProcessOneDecomposition(bpIO, bpReader, ncid, decompname, wfiles, iosysid, mpirank, nproc, var.nctype);
			} else {
        		decomp = GetNewDecomposition(decomp_map, decompname, bpIO, bpReader, 
											 ncid, wfiles, var.nctype, iosysid, mpirank, nproc);
			} 
		}
		if (frame_id<0) frame_id = 0;
        if (wfiles[0] < nblocks_per_step)
        {
			/* different decompositions at different frames */	
            if (var.is_timed)
                PIOc_setframe(ncid, var.nc_varid, frame_id);
			if (fillval_exist) {
            	ret = PIOc_write_darray(ncid, var.nc_varid, decomp.ioid, (PIO_Offset)nelems,
						d.data(), fillval_id); 
			} else {
            	ret = PIOc_write_darray(ncid, var.nc_varid, decomp.ioid, (PIO_Offset)nelems,
						d.data(), NULL); 
			}
        }
		if (mem_opt) {
			PIOc_sync(ncid); 
			PIOc_freedecomp(iosysid,decomp.ioid);
		}
        TimerStop(write);
    }

    return ret;
}

int ConvertVariableDarray(adios2::IO *bpIO, adios2::Engine *bpReader, 
		std::string varname, int ncid, Variable& var,
        std::vector<int>& wfiles, DecompositionMap& decomp_map, int nblocks_per_step, int iosysid,
		MPI_Comm comm, int mpirank, int nproc, int mem_opt)
{
	TimerStart(read);
    int ret = 0;

    std::string v_type = bpIO->VariableType(varname);
    if (v_type.empty()) {
        return 1;
    }
#define declare_template_instantiation(T)              \
    else if (v_type == adios2::GetType<T>())   \
    {                                                  \
        adios2::Variable<T> v_base; \
        std::vector<T> v_value; \
		return adios2_ConvertVariableDarray(&v_base,v_value, \
                            bpIO, bpReader, varname, ncid, var, \
                            wfiles, decomp_map, nblocks_per_step, iosysid, \
                            comm, mpirank, nproc, mem_opt); \
    }
    ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

}

/*
 * Assumes a BP folder with name "infilename.dir" and 
 * all the files in the folder are bp files. It also 
 * assumes the file extensions are "infilename.bp.X" 
 * where X is 0 to N-1.
 */
int GetNumOfFiles(string infilename) 
{
	int file_count = 0;
	string foldername = infilename + ".dir/";
	DIR* dirp = opendir(foldername.c_str());
	if (!dirp) {
		fprintf(stderr, "Folder %s does not exist.\n",foldername.c_str());
		return -1;
	}
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
		if (dp->d_type==DT_REG && strstr(dp->d_name,".bp.")!=NULL) file_count++; 
    }
    closedir(dirp);	
	return file_count;
}

std::string ExtractFilename(std::string pathname)
{
	size_t pos = pathname.find_last_of("/\\");
	if (pos==std::string::npos) {
		return pathname;
	} else {
		return pathname.substr(pos+1);
	}
}

std::string ExtractPathname(std::string pathname)
{
	size_t pos = pathname.find_last_of("/\\");
	if (pos==std::string::npos) {
		return "./";
	} else {
		return pathname.substr(0,pos);
	}
}

void ConvertBPFile(string infilepath, string outfilename, int pio_iotype, int iosysid, MPI_Comm comm, int mpirank, int nproc, int mem_opt)
{
	// ADIOS_FILE **infile = NULL;
	int num_infiles = 0; 
	int ncid = -1;
	int n_bp_writers;
	int ret = 0;

	// Initialization of the class factory	
   	adios2::ADIOS adios(comm, adios2::DebugON);

	try {

		int save_imax = pio_get_imax();

		/*
		 * This assumes we are running the program in the same folder where 
		 * the BP folder exists.
		 *
		 */
		std::string foldername   = ExtractPathname(infilepath);
		std::string basefilename = ExtractFilename(infilepath);

		/*
		 * Get the number of files in BP folder. 
		 * This operation assumes that the BP folder contains only the 
		 * BP files. 
		 */
		int n_bp_files = GetNumOfFiles(infilepath);
		if (n_bp_files<0) 
			throw std::runtime_error("Input folder doesn't exist.\n");

		if (nproc>n_bp_files) {
			if (debug_out) std::cout << "ERROR: nproc (" << nproc << ") is greater than #files (" << n_bp_files << ")" << std::endl;
			throw std::runtime_error("Use fewer processors.\n");
		}

		/* Number of BP file writers != number of converter processes here */
        std::vector<int> wfiles;
        wfiles = AssignWriteRanks(n_bp_files,comm, mpirank, nproc);
		if (debug_out) {
        	for (auto nb: wfiles)
            	printf("Myrank: %d File id: %d\n",mpirank,nb);
		}

		num_infiles = wfiles.size()+1; //  infile[0] is opened by all processors

		/* Now allocate IO and Engine for each file */
        adios2::IO *bpIO = new adios2::IO[num_infiles]; 
		adios2::Engine *bpReader = new adios2::Engine[num_infiles];

		// infile = (ADIOS_FILE **)malloc(num_infiles*sizeof(ADIOS_FILE *));
		// if (!infile) 
		// throw std::runtime_error("Cannot allocate space for infile array.\n");
		
		string file0 = infilepath + ".dir/" + basefilename + ".0";
		if (debug_out) printf("FILE0: %s\n",file0.c_str()); fflush(stdout);

		bpIO[0] = adios.DeclareIO(file0 + "_0");
		bpReader[0] = bpIO[0].Open(file0, adios2::Mode::Read, MPI_COMM_SELF);

	    adios2::Variable<int> bpNProc = bpIO[0].InquireVariable<int>("/__pio__/info/nproc");
		bpReader[0].Get(bpNProc,&n_bp_writers,adios2::Mode::Sync);

		// int ret = adios_schedule_read(infile[0], NULL, "/__pio__/info/nproc", 0, 1, &n_bp_writers);
		// if (ret)
		// 	throw std::runtime_error("Invalid BP file: missing '/__pio__/info/nproc' variable\n");
		// adios_perform_reads(infile[0], 1);
		
		if (debug_out) printf("AFTER AFTER FILE0: %s %d\n",file0.c_str(),n_bp_writers); fflush(stdout);
		if (n_bp_writers!=n_bp_files) {
			if (debug_out) std::cout  << "WARNING: #writers (" 
					<< n_bp_writers << ") != #files (" 
					<< n_bp_files << ")" << std::endl;
		} else {
			if (debug_out) std::cout << "n_bp_writers: " << n_bp_writers 
					  << " n_bp_files: " << n_bp_files << std::endl;
		}

		/*
		 * Open the BP files. 
		 * basefilename.bp.0 is opened by all the nodes. It contains all of the variables 
		 * and attributes. Each node then opens the files assigned to that node. 
		 */
		for (int i=1;i<=wfiles.size();i++) {
			// ostringstream ss; ss << wfiles[i-1];
			// string fileid_str = ss.str(); 
			char ss[64];
			sprintf(ss,"%d",wfiles[i-1]);
			string fileid_str(ss);

			string filei = infilepath + ".dir/" + basefilename + "." + fileid_str;

			if (debug_out) std::cout << "myrank " << mpirank << " file: " << filei << std::endl;

			try {
				bpIO[i] = adios.DeclareIO(filei + "_" + std::to_string(i));
				bpReader[i] = bpIO[i].Open(filei, adios2::Mode::Read, MPI_COMM_SELF);
			} catch (const std::exception &exc ) {
				std::cerr << exc.what();
			}

			// infile[i] = adios_read_open_file(filei.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_SELF); 
			// if (debug_out) std::cout << "myrank " << mpirank << " file: " << filei << std::endl;
		}

		/* Create output file */
		TimerStart(write);
		/* 
			Use NC_64BIT_DATA instead of PIO_64BIT_OFFSET. Some output files will have variables 
			that require more than 4GB storage. 
		*/
		ret = PIOc_createfile(iosysid, &ncid, &pio_iotype, outfilename.c_str(), NC_64BIT_DATA); 
		TimerStop(write);
		if (ret)
			throw std::runtime_error("Could not create output file " + outfilename + "\n");

		/* First process decompositions */
		DecompositionMap decomp_map;
		if (!mem_opt) decomp_map = ProcessDecompositions(bpIO, bpReader, ncid, wfiles,iosysid,comm, mpirank, nproc);

		/* Process the global fillmode */
		ProcessGlobalFillmode(bpIO, bpReader, ncid);

		/* Next process dimensions */
		DimensionMap dimension_map = ProcessDimensions(bpIO, bpReader, ncid,comm, mpirank, nproc);

		/* For each variable, define a variable with PIO */
		VariableMap vars_map = ProcessVariableDefinitions(bpIO, bpReader, ncid, dimension_map, comm, mpirank, nproc);

		/* Process the global attributes */
		ProcessGlobalAttributes(bpIO, bpReader, ncid, dimension_map, vars_map);

		PIOc_enddef(ncid);

		/* 
		 * For each variable,read in the data 
		 * with ADIOS then write it out with PIO 
		 */
	 	std::map<std::string, adios2::Params> a2_vi = bpIO[0].AvailableVariables();
    	for (std::map<std::string, adios2::Params>::iterator a2_iter = a2_vi.begin(); a2_iter!=a2_vi.end(); a2_iter++)
		{
			string v = a2_iter->first; // infile[0]->var_namelist[i];
			if (v.find("/__") == string::npos)
			{
				/* For each variable, read with ADIOS then write with PIO */
				if (!mpirank && debug_out) cout << "Convert variable: " << v << endl;

				if (v.find("decomp_id/")==string::npos && 
					v.find("frame_id/")==string::npos &&
					v.find("fillval_id/")==string::npos) {
					Variable& var = vars_map[v];

					TimerStart(read);
					string attname = string(a2_iter->first) + "/__pio__/ncop";
					size_t asize;
					char **ncop = NULL;
					std::string atype;
					adios_get_attr_a2(bpIO, bpReader, (char*)attname.c_str(), &atype, &asize, (void**)&ncop);
					TimerStop(read);

					std::string op(ncop[0]);
					if (op == "put_var") {
						if (var.is_timed) {
							if (debug_out) printf("ConvertVariableTimedPutVar: %d\n",mpirank); fflush(stdout);
							ConvertVariableTimedPutVar(bpIO, bpReader, wfiles, v, ncid, var, n_bp_writers, comm, mpirank, nproc);
						} else {
							if (debug_out) printf("ConvertVariablePutVar: %d\n",mpirank); fflush(stdout);
							ConvertVariablePutVar(bpIO, bpReader, wfiles, v, ncid, var, mpirank, nproc);
						}
					} else if (op == "darray") {
						/* Variable was written with pio_write_darray() with a decomposition */
						if (debug_out) printf("ConvertVariableDarray: %d\n",mpirank); fflush(stdout);
						ConvertVariableDarray(bpIO, bpReader, v, ncid, var, wfiles, decomp_map, n_bp_writers, iosysid, comm, mpirank, nproc, mem_opt);
					} else {
						if (!mpirank && debug_out)
							cout << "  WARNING: unknown operation " << op << ". Will not process this variable\n";
					}
					if (ncop) free(ncop);
				}
			}
			FlushStdout_nm(comm);
			PIOc_sync(ncid); /* FIXME: flush after each variable until development is done. Remove for efficiency */
		}
		TimerStart(write);
	
 		for (std::map<std::string,Decomposition>::iterator it=decomp_map.begin(); it!=decomp_map.end(); ++it) {
		     Decomposition d = it->second;
		     int err_code = PIOc_freedecomp(iosysid, d.ioid);
		     if (err_code!=0) {
		           printf("ERROR: PIOc_freedecomp: %d\n",err_code);
		             fflush(stdout);
		     }
		}
	
		pio_set_imax(save_imax);
	
		ret = PIOc_sync(ncid);
		ret = PIOc_closefile(ncid);
		TimerStop(write);
		TimerStart(read);

		MPI_Barrier(comm);
		// for (int i=0;i<num_infiles;i++) 
		// 	adios_read_close(infile[i]);
		TimerStop(read);

		return;
   } catch (std::exception &e) {
        if (ncid > -1)
            PIOc_closefile(ncid);
		// if (infile) {
		// 	for (int i=0;i<num_infiles;i++) 
		// 		adios_read_close(infile[i]);
		// }
        throw e;
    }
}

void usage_nm(string prgname)
{
        cout << "Usage: " << prgname << " bp_file  nc_file  pio_io_type memory_opt\n";
        cout << "   bp file   :  data produced by PIO with ADIOS format\n";
        cout << "   nc file   :  output file name after conversion\n";
        cout << "   memory opt:  (0/1) reduce memory usage. Reduces execution speed.\n";
        cout << "   pio format:  output PIO_IO_TYPE. Supported parameters:\n";
        cout << "                pnetcdf  netcdf  netcdf4c  netcdf4p   or:\n";
        cout << "                   1       2        3         4\n";
}

enum PIO_IOTYPE GetIOType_nm(string t)
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

int ConvertBPToNC(string infilepath, string outfilename, string piotype, int mem_opt, MPI_Comm comm_in)
{
	int ret = 0;
	int iosysid = 0;

	MPI_Comm comm   = comm_in;
	MPI_Comm w_comm = comm_in;
	int mpirank, w_mpirank;
	int nproc, w_nproc;

    MPI_Comm_set_errhandler(w_comm, MPI_ERRORS_RETURN);
    MPI_Comm_rank(w_comm, &w_mpirank);
    MPI_Comm_size(w_comm, &w_nproc);

	SetDebugOutput(1);

	if (mem_opt) 
		printf("INFO: Option selected to reduce memory usage. Execution time will likely increase.\n");
	else
		printf("INFO: Reduce memory usage option is set to 0.\n");
	fflush(stdout);

	/* 
 	 * Check if the number of nodes is less than or equal to the number of BP files.
 	 * If not, create a new comm.
 	 *
 	 */
	int num_files = GetNumOfFiles(infilepath); 
	int io_proc   = 0;
	if (num_files<w_nproc) {
		fprintf(stderr, "Warning: #files: %d < #procs: %d\n",num_files,w_nproc); 
		fflush(stderr);
		if (w_mpirank<num_files) /* I/O nodes */
			io_proc = 1; 
		MPI_Comm_split(w_comm,io_proc,w_mpirank,&comm);
    	MPI_Comm_rank(comm, &mpirank);
    	MPI_Comm_size(comm, &nproc);
	} else {
		io_proc = 1;
		comm    = w_comm;
		mpirank = w_mpirank;
		nproc   = w_nproc;
	}
	
    TimerInitialize_nm();

    try {
		if (io_proc) {
        	enum PIO_IOTYPE pio_iotype = GetIOType_nm(piotype);
        	iosysid = InitPIO(comm,mpirank,nproc);
        	ConvertBPFile(infilepath, outfilename, pio_iotype, iosysid, comm, mpirank, nproc, mem_opt);
        	PIOc_finalize(iosysid);
        	TimerReport_nm(comm);
		}
		MPI_Barrier(w_comm);
    } catch (std::invalid_argument &e) {
        std::cout << e.what() << std::endl;
        return 2;
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return 3;
    }

    TimerFinalize_nm();
    return ret;
}

#ifdef __cplusplus
extern "C" {
#endif

int C_API_ConvertBPToNC(const char *infilepath, const char *outfilename, const char *piotype, int mem_opt, MPI_Comm comm_in)
{
    return ConvertBPToNC(string(infilepath), string(outfilename), string(piotype), mem_opt, comm_in);
}

#ifdef __cplusplus
}
#endif

