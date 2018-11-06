###-------------------------------------------------------------------------###
### CMakeList.txt for examples using pio
###-------------------------------------------------------------------------###

if (WITH_ADIOS)
  find_package (ADIOS "1.10.0")
  if (ADIOS_FOUND)
    # ADD_SUBDIRECTORY(adios2pio)
    ADD_SUBDIRECTORY(adios2pio-nm)
  endif (ADIOS_FOUND)
endif(WITH_ADIOS)

if (WITH_ADIOS2)
  find_package (ADIOS "1.13.1")
  if (ADIOS_FOUND)
    # ADD_SUBDIRECTORY(adios2pio)
    ADD_SUBDIRECTORY(adios2pio-nm)
  endif (ADIOS_FOUND)
endif(WITH_ADIOS2)
