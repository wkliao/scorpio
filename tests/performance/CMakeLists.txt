SET(SRC pioperformance_rearr.F90)
ADD_DEFINITIONS(${PIO_DEFINITIONS} -D_PIO1)

# Use to trace hdf5 errors with a special netcdf build
OPTION(HDF5_LOGGING "Turn on hdf5 logging (requires instrumented netcdf4)" OFF)
if(HDF5_LOGGING) 
  ADD_DEFINITIONS(-DLOGGING)
endif()

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS} )
LINK_DIRECTORIES(${PIO_LIB_DIR} )

ADD_EXECUTABLE(pioperf_rearr ${SRC})

if(${PIO_BUILD_TIMING})
  TARGET_LINK_LIBRARIES(pioperf_rearr pio timing)
ELSE()
  TARGET_LINK_LIBRARIES(pioperf_rearr pio)
ENDIF()

TARGET_LINK_LIBRARIES(pioperf_rearr ${PnetCDF_Fortran_LIBRARIES})
TARGET_LINK_LIBRARIES(pioperf_rearr ${NetCDF_Fortran_LIBRARIES})

#ADD_TEST(pioperf_rearr pioperf_rearr)
