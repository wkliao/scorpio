static int needs_to_write_decomp(file_desc_t *file, int ioid)
{
    int ret = 1; // Yes
    for (int i = 0; i < file->n_written_ioids; i++)
    {
        if (file->written_ioids[i] == ioid)
        {
            ret = 0; // No
            break;
        }
    }
    return ret;
}

static int register_decomp(file_desc_t *file, int ioid)
{
    if (file->n_written_ioids >= ADIOS_PIO_MAX_DECOMPS)
        return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__,"ADIOS_PIO_MAX_DECOMPS exceeded: %d",file->n_written_ioids);

    file->written_ioids[file->n_written_ioids] = ioid;
    ++file->n_written_ioids;

    return PIO_NOERR;
}

#define ADIOS_CONVERT_ARRAY(array,arraylen,from_type,to_type,ierr,buf) \
{ \
	from_type *d = (from_type*)array; \
	to_type *f   = (to_type*)malloc(arraylen*sizeof(to_type)); \
	if (f) { \
		for (int i=0; i<arraylen; ++i) \
			f[i] = (to_type)d[i]; \
		buf = f; \
	} else { \
		ierr = 1; \
	} \
}

#define ADIOS_CONVERT_FROM(FROM_TYPE_ID,from_type) \
{ \
   	if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_DOUBLE) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,double,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_FLOAT) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,float,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_REAL) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,float,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_INT) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,int,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UINT) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,unsigned int,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_SHORT) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,short int,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_USHORT) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,unsigned short int,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_INT64) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,int64_t,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UINT64) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,uint64_t,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_CHAR) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,char,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_BYTE) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,char,*ierr,buf); \
	} else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UBYTE) { \
		ADIOS_CONVERT_ARRAY(array,arraylen,from_type,unsigned char,*ierr,buf); \
	} \
}

static void *PIOc_convert_buffer_adios(file_desc_t *file, io_desc_t *iodesc, 
								adios_var_desc_t *av, void *array, int arraylen, 
								int *ierr) 
{
	void *buf = array;
	*ierr = 0;

	ADIOS_CONVERT_FROM(PIO_DOUBLE,double);
	ADIOS_CONVERT_FROM(PIO_FLOAT,float);
	ADIOS_CONVERT_FROM(PIO_REAL,float);
	ADIOS_CONVERT_FROM(PIO_INT,int);
	ADIOS_CONVERT_FROM(PIO_UINT,unsigned int);
	ADIOS_CONVERT_FROM(PIO_SHORT,short int);
	ADIOS_CONVERT_FROM(PIO_USHORT,unsigned short int);
	ADIOS_CONVERT_FROM(PIO_INT64,int64_t);
	ADIOS_CONVERT_FROM(PIO_UINT64,uint64_t);
	ADIOS_CONVERT_FROM(PIO_CHAR,char);
	ADIOS_CONVERT_FROM(PIO_BYTE,char);
	ADIOS_CONVERT_FROM(PIO_UBYTE,unsigned char);

	return buf;
}

#define ADIOS_COPY_ONE(temp_buf,array,var_type) \
{ \
	temp_buf = (var_type*)malloc(2*sizeof(var_type)); \
    memcpy(temp_buf,array,sizeof(var_type)); \
}

void *PIOc_copy_one_element_adios(void *array, io_desc_t *iodesc) 
{
	void *temp_buf = NULL;
	if (iodesc->piotype==PIO_DOUBLE) {
		ADIOS_COPY_ONE(temp_buf,array,double);
	} else if (iodesc->piotype==PIO_FLOAT || iodesc->piotype==PIO_REAL) {
		ADIOS_COPY_ONE(temp_buf,array,float);
	} else if (iodesc->piotype==PIO_INT || iodesc->piotype==PIO_UINT) {
		ADIOS_COPY_ONE(temp_buf,array,int);
	} else if (iodesc->piotype==PIO_SHORT || iodesc->piotype==PIO_USHORT) {
		ADIOS_COPY_ONE(temp_buf,array,short int);
	} else if (iodesc->piotype==PIO_INT64 || iodesc->piotype==PIO_UINT64) {
		ADIOS_COPY_ONE(temp_buf,array,int64_t);
	} else if (iodesc->piotype==PIO_CHAR || iodesc->piotype==PIO_BYTE || iodesc->piotype==PIO_UBYTE) {
		ADIOS_COPY_ONE(temp_buf,array,char);
	}
	return temp_buf;
}

static void PIOc_write_decomp_adios(file_desc_t *file, int ioid)
{
    io_desc_t *iodesc = pio_get_iodesc_from_id(ioid);
    char name[32];
    sprintf(name, "/__pio__/decomp/%d", ioid);

   	adios2_type type = adios2_type_int32_t; 
   	if (sizeof(PIO_Offset) == 8)
   		type = adios2_type_int64_t;

	size_t av_count[1];
	if (iodesc->maplen>1) {
		av_count[0]  = (size_t)iodesc->maplen;

		adios2_variable *variableH = adios2_inquire_variable(file->ioH,name);
		if (variableH==NULL) {
			variableH = adios2_define_variable(file->ioH, name, type,
												1, NULL, NULL, av_count, 
												adios2_constant_dims_true);
		}
   		adios2_put(file->engineH, variableH, iodesc->map, adios2_mode_sync);
	} else if (iodesc->maplen==0) { // Handle the case where maplen is 0
		long mapbuf[2];
		mapbuf[0] = 0; 
		mapbuf[1] = 0;
		av_count[0] = (size_t)2;

		adios2_variable *variableH = adios2_inquire_variable(file->ioH,name);
		if (variableH==NULL) {
			variableH = adios2_define_variable(file->ioH, name, type,
												1, NULL, NULL, av_count, 
												adios2_constant_dims_true);
		}
       	adios2_put(file->engineH, variableH, mapbuf, adios2_mode_sync);
	} else { // Handle the case where maplen is 1
		int maplen   = iodesc->maplen+1; 
		void *mapbuf = NULL;
		if (type==adios2_type_int32_t) {
			mapbuf = (int*)malloc(sizeof(int)*maplen);	
			((int*)mapbuf)[0] = iodesc->map[0];
			((int*)mapbuf)[1] = 0;
		} else {
			mapbuf = (long*)malloc(sizeof(long)*maplen);	
			((long*)mapbuf)[0] = iodesc->map[0];
			((long*)mapbuf)[1] = 0;
		}

		av_count[0] = (size_t)maplen;
		adios2_variable *variableH = adios2_inquire_variable(file->ioH,name);
		if (variableH==NULL) {
			variableH = adios2_define_variable(file->ioH, name, type,
												1, NULL, NULL, av_count, 
												adios2_constant_dims_true);
		}
   		adios2_put(file->engineH, variableH, mapbuf, adios2_mode_sync);
		free(mapbuf);
	}

	/* ADIOS: assume all procs are also IO tasks */
    if (file->adios_iomaster == MPI_ROOT)
   	{
		char att_name[PIO_MAX_NAME];

		sprintf(att_name,"%s/piotype",name);
		if (adios2_inquire_attribute(file->ioH,att_name)==NULL) 
			adios2_define_attribute(file->ioH,att_name,adios2_type_int32_t,&iodesc->piotype);

		sprintf(att_name,"%s/ndims",name);
		if (adios2_inquire_attribute(file->ioH,att_name)==NULL) 
			adios2_define_attribute(file->ioH,att_name,adios2_type_int32_t,&iodesc->ndims);

		sprintf(att_name,"%s/dimlen",name);
		if (adios2_inquire_attribute(file->ioH,att_name)==NULL) 
			adios2_define_attribute_array(file->ioH,att_name,adios2_type_int32_t,iodesc->dimlen,iodesc->ndims);
   	}
}

static int PIOc_write_darray_adios(file_desc_t *file, int varid, int ioid, io_desc_t *iodesc,
        PIO_Offset arraylen, void *array, void *fillvalue)
{
    int ierr = PIO_NOERR;
    if (varid < 0 || varid >= file->num_vars)
        return pio_err(file->iosystem, file, PIO_EBADID, __FILE__, __LINE__,"varid is wrong: %d %d",varid,file->num_vars);

    adios_var_desc_t * av = &(file->adios_vars[varid]);

	// Handle the case where there is zero or one array element 
	void *temp_buf = NULL;
	if (arraylen==0) {
		arraylen = 2;
		temp_buf = (int64_t*)malloc(sizeof(int64_t)*arraylen);
		memset(temp_buf,0,sizeof(int64_t)*arraylen);
		array = temp_buf;
	} else if (arraylen==1) { 
		arraylen = 2;
		temp_buf = PIOc_copy_one_element_adios(array,iodesc);
		array = temp_buf;
	} 

    if (av->adios_varid == NULL)
    {
        /* First we need to define the variable now that we know it's decomposition */
        adios2_type atype = av->adios_type;
		size_t av_count[1];
		av_count[0] = (size_t)arraylen;
        av->adios_varid = adios2_define_variable(file->ioH,av->name,atype,
												1,NULL,NULL,av_count,
												adios2_constant_dims_true);

		/* different decompositions at different frames */
		char name_varid[256];
		sprintf(name_varid,"decomp_id/%s",av->name);
		av_count[0] = 1;
		av->decomp_varid = adios2_inquire_variable(file->ioH,name_varid);
		if (av->decomp_varid==NULL) {
			av->decomp_varid = adios2_define_variable(file->ioH,name_varid,adios2_type_int32_t,
												1,NULL,NULL,av_count,
                                                adios2_constant_dims_true);
		}

		sprintf(name_varid,"frame_id/%s",av->name);
		av_count[0] = 1;
		av->frame_varid = adios2_inquire_variable(file->ioH,name_varid);
		if (av->frame_varid==NULL) {
			av->frame_varid = adios2_define_variable(file->ioH,name_varid,adios2_type_int32_t,
												1,NULL,NULL,av_count,
                                                adios2_constant_dims_true);
		}

		sprintf(name_varid,"fillval_id/%s",av->name);
		av_count[0] = 1;
		av->fillval_varid = adios2_inquire_variable(file->ioH,name_varid);
		if (av->fillval_varid==NULL) {
			av->fillval_varid = adios2_define_variable(file->ioH,name_varid,atype,
												1,NULL,NULL,av_count,
                                                adios2_constant_dims_true);
		}
		
        if (file->adios_iomaster == MPI_ROOT)
        { /* TAHSIN: Some of the codes were moved to pio_nc.c */
            char decompname[32], att_name[64];
            sprintf(decompname, "%d", ioid);

			sprintf(att_name,"%s/__pio__/decomp",av->name);
			if (adios2_inquire_attribute(file->ioH,att_name)==NULL) 
				adios2_define_attribute(file->ioH,att_name,adios2_type_string,decompname);

			sprintf(att_name,"%s/__pio__/ncop",av->name);
			if (adios2_inquire_attribute(file->ioH,att_name)==NULL) 
				adios2_define_attribute(file->ioH,att_name,adios2_type_string,"darray");
        }
    }

	/* Check if we need to write the decomposition. Write it */
	if (needs_to_write_decomp(file, ioid))
    {
        PIOc_write_decomp_adios(file,ioid);
        register_decomp(file, ioid);
    }

    /* ACME history data special handling: down-conversion from double to float */
    void *buf = array;
    int buf_needs_free = 0;
    if (iodesc->piotype != av->nc_type)
    {
		buf = PIOc_convert_buffer_adios(file,iodesc,av,array,arraylen,&ierr);
		if (ierr!=0) {
       		if (file->adios_iomaster == MPI_ROOT)
                    LOG((2,"Darray '%s' decomp type is double but the target type is float. "
                            "ADIOS cannot do type conversion because memory could not be allocated."
                            "Therefore the data will be corrupt for this variable in the .bp output\n",
                            file->adios_vars[varid].name));
        } else {
			buf_needs_free = 1;
		}
    }

	adios2_put(file->engineH, av->adios_varid, buf, adios2_mode_sync);

	/* NOTE: PIOc_setframe with different decompositions          */
	/* different decompositions at different frames and fillvalue */
	if (fillvalue!=NULL) { 
		adios2_put(file->engineH, av->fillval_varid, fillvalue, adios2_mode_sync);
	} else {
		ioid = -ioid;
	}
	adios2_put(file->engineH, av->decomp_varid, &ioid, adios2_mode_sync);
	adios2_put(file->engineH, av->frame_varid, &(file->varlist[varid].record), adios2_mode_sync);

    if (buf_needs_free) free(buf);
	if (temp_buf!=NULL) free(temp_buf);

    return PIO_NOERR;
}

