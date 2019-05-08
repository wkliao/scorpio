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
        return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);

    file->written_ioids[file->n_written_ioids] = ioid;
    ++file->n_written_ioids;

    return PIO_NOERR;
}

static int PIOc_write_decomp_adios(file_desc_t *file, int ioid)
{
    io_desc_t *iodesc = pio_get_iodesc_from_id(ioid);
    char name[PIO_MAX_NAME], ldim[PIO_MAX_NAME];
    sprintf(name, "/__pio__/decomp/%d", ioid);

    enum ADIOS_DATATYPES type = adios_integer;
    if (sizeof(PIO_Offset) == 8)
        type = adios_long;

    if (iodesc->maplen != 1)
    {
        sprintf(ldim, "%d", iodesc->maplen);
        int64_t vid = adios_define_var(file->adios_group, name, "", type, ldim, "", "");
        adios_write_byid(file->adios_fh, vid, iodesc->map);
    }
    else /* Handle the case where maplen is 1 */
    {
        int maplen = iodesc->maplen + 1;
        char *mapbuf = NULL;

        sprintf(ldim, "%d", maplen);

        if (type == adios_integer)
        {
            int *temp_mapbuf = (int*)malloc(sizeof(int)*maplen);
            if (temp_mapbuf == NULL)
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            temp_mapbuf[0] = iodesc->map[0];
            temp_mapbuf[1] = 0;
            mapbuf = (char*)temp_mapbuf;
        }
        else
        {
            long *temp_mapbuf = (long*)malloc(sizeof(long)*maplen);
            if (temp_mapbuf == NULL)
                return pio_err(NULL, NULL, PIO_ENOMEM, __FILE__, __LINE__);
            temp_mapbuf[0] = iodesc->map[0];
            temp_mapbuf[1] = 0;
            mapbuf = (char*)temp_mapbuf;
        }

        int64_t vid = adios_define_var(file->adios_group, name, "", type, ldim, "", "");
        adios_write_byid(file->adios_fh, vid, mapbuf);

        if (mapbuf != NULL)
            free(mapbuf);
    }

    /* ADIOS: assume all procs are also IO tasks */
    if (file->adios_iomaster == MPI_ROOT)
    {
        adios_define_attribute_byvalue(file->adios_group, "piotype", name, adios_integer, 1, &iodesc->piotype);
        adios_define_attribute_byvalue(file->adios_group, "ndims", name, adios_integer, 1, &iodesc->ndims);
        adios_define_attribute_byvalue(file->adios_group, "dimlen", name, adios_integer, iodesc->ndims, iodesc->dimlen);
    }

    return PIO_NOERR;
}

#define ADIOS_CONVERT_ARRAY(array, arraylen, from_type, to_type, ierr, buf) \
{ \
    from_type *d = (from_type*)array; \
    to_type *f = (to_type*)malloc(arraylen * sizeof(to_type)); \
    if (f) { \
        for (int i = 0; i < arraylen; ++i) \
            f[i] = (to_type)d[i]; \
        buf = f; \
    } \
    else \
    { \
        ierr = 1; \
    } \
}

#define ADIOS_CONVERT_FROM(FROM_TYPE_ID, from_type) \
{ \
    if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_DOUBLE) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, double, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_FLOAT) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, float, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_REAL) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, float, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_INT) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, int, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UINT) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, unsigned int, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_SHORT) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, short int, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_USHORT) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, unsigned short int, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_INT64) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, int64_t, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UINT64) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, uint64_t, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_CHAR) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, char, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_BYTE) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, char, *ierr, buf); \
    } \
    else if (iodesc->piotype == FROM_TYPE_ID && av->nc_type == PIO_UBYTE) \
    { \
        ADIOS_CONVERT_ARRAY(array, arraylen, from_type, unsigned char, *ierr, buf); \
    } \
}

static void *PIOc_convert_buffer_adios(file_desc_t *file, io_desc_t *iodesc,
                                       adios_var_desc_t *av, void *array, int arraylen,
                                       int *ierr)
{
    void *buf = array;
    *ierr = 0;

    ADIOS_CONVERT_FROM(PIO_DOUBLE, double);
    ADIOS_CONVERT_FROM(PIO_FLOAT, float);
    ADIOS_CONVERT_FROM(PIO_REAL, float);
    ADIOS_CONVERT_FROM(PIO_INT, int);
    ADIOS_CONVERT_FROM(PIO_UINT, unsigned int);
    ADIOS_CONVERT_FROM(PIO_SHORT, short int);
    ADIOS_CONVERT_FROM(PIO_USHORT, unsigned short int);
    ADIOS_CONVERT_FROM(PIO_INT64, int64_t);
    ADIOS_CONVERT_FROM(PIO_UINT64, uint64_t);
    ADIOS_CONVERT_FROM(PIO_CHAR, char);
    ADIOS_CONVERT_FROM(PIO_BYTE, char);
    ADIOS_CONVERT_FROM(PIO_UBYTE, unsigned char);

    return buf;
}

#define ADIOS_COPY_ONE(temp_buf, array, var_type) \
{ \
    temp_buf = (var_type*)malloc(2 * sizeof(var_type)); \
    assert(temp_buf != NULL); \
    memcpy(temp_buf, array, sizeof(var_type)); \
}

void *PIOc_copy_one_element_adios(void *array, adios_var_desc_t *av)
{
    void *temp_buf = NULL;
    if (av->nc_type == PIO_DOUBLE)
    {
        ADIOS_COPY_ONE(temp_buf, array, double);
    }
    else if (av->nc_type == PIO_FLOAT || av->nc_type == PIO_REAL)
    {
        ADIOS_COPY_ONE(temp_buf, array, float);
    }
    else if (av->nc_type == PIO_INT || av->nc_type == PIO_UINT)
    {
        ADIOS_COPY_ONE(temp_buf, array, int);
    }
    else if (av->nc_type == PIO_SHORT || av->nc_type == PIO_USHORT)
    {
        ADIOS_COPY_ONE(temp_buf, array, short int);
    }
    else if (av->nc_type == PIO_INT64 || av->nc_type == PIO_UINT64)
    {
        ADIOS_COPY_ONE(temp_buf, array, int64_t);
    }
    else if (av->nc_type == PIO_CHAR || av->nc_type == PIO_BYTE || av->nc_type == PIO_UBYTE)
    {
        ADIOS_COPY_ONE(temp_buf, array, char);
    }

    return temp_buf;
}

static int PIOc_write_darray_adios(file_desc_t *file, int varid, int ioid,
                                   io_desc_t *iodesc, PIO_Offset arraylen,
                                   void *array, void *fillvalue)
{
    int ierr = PIO_NOERR;
    if (varid < 0 || varid >= file->num_vars)
        return pio_err(file->iosystem, file, PIO_EBADID, __FILE__, __LINE__);

    adios_var_desc_t *av = &(file->adios_vars[varid]);

    void *temp_buf = NULL;
    if (arraylen == 1) /* Handle the case where there is one array element */
    {
        arraylen++;
        temp_buf = PIOc_copy_one_element_adios(array, av);
        array = temp_buf;
    }

    if (av->adios_varid == 0)
    {
        /* First we need to define the variable now that we know it's decomposition */
        char ldims[PIO_MAX_NAME];
        sprintf(ldims, "%lld", arraylen);
        enum ADIOS_DATATYPES atype = av->adios_type;

        av->adios_varid = adios_define_var(file->adios_group, av->name, "", atype, ldims, "", "");

        /* Different decompositions at different frames */
        char name_varid[PIO_MAX_NAME];
        sprintf(name_varid, "decomp_id/%s", av->name);
        av->decomp_varid = adios_define_var(file->adios_group, name_varid, "", adios_integer, "1", "", "");
        sprintf(name_varid, "frame_id/%s", av->name);
        av->frame_varid = adios_define_var(file->adios_group, name_varid, "", adios_integer, "1", "", "");
        sprintf(name_varid, "fillval_id/%s", av->name);
        av->fillval_varid = adios_define_var(file->adios_group, name_varid, "", atype, "1", "", "");

        if (file->adios_iomaster == MPI_ROOT)
        {
            /* Some of the codes were moved to pio_nc.c */
            char decompname[PIO_MAX_NAME];
            sprintf(decompname, "%d", ioid);
            adios_define_attribute(file->adios_group, "__pio__/decomp", av->name, adios_string, decompname, NULL);
            adios_define_attribute(file->adios_group, "__pio__/ncop", av->name, adios_string, "darray", NULL);
         }
    }

    /* PIOc_setframe with different decompositions */
    if (needs_to_write_decomp(file, ioid))
    {
        ierr = PIOc_write_decomp_adios(file, ioid);
        if (ierr != PIO_NOERR)
            return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);
        ierr = register_decomp(file, ioid);
        if (ierr != PIO_NOERR)
            return pio_err(NULL, NULL, ierr, __FILE__, __LINE__);
    }

    /* E3SM history data special handling: down-conversion from double to float */
    void *buf = array;
    int buf_needs_free = 0;
    if (iodesc->piotype != av->nc_type)
    {
        buf = PIOc_convert_buffer_adios(file, iodesc, av, array, arraylen, &ierr);
        if (ierr != 0)
        {
            if (file->adios_iomaster == MPI_ROOT)
                LOG((2, "Darray '%s' decomp type is double but the target type is float. "
                        "ADIOS cannot do type conversion because memory could not be allocated. "
                        "Therefore the data will be corrupt for this variable in the .bp output.",
                        file->adios_vars[varid].name));
        }
        else
        {
            buf_needs_free = 1;
        }
    }

    adios_write_byid(file->adios_fh, av->adios_varid, buf);

    /* Different decompositions at different frames and fillvalue */
    if (fillvalue != NULL)
    {
        adios_write_byid(file->adios_fh, av->fillval_varid, fillvalue);
    }
    else
    {
        ioid = -ioid;
    }
    adios_write_byid(file->adios_fh, av->decomp_varid, &ioid);
    adios_write_byid(file->adios_fh, av->frame_varid, &(file->varlist[varid].record));

    if (buf_needs_free)
        free(buf);

    if (temp_buf != NULL)
        free(temp_buf);

    return PIO_NOERR;
}

