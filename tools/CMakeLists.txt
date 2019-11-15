###-------------------------------------------------------------------------###
### CMakeList.txt for examples using pio
###-------------------------------------------------------------------------###

if (WITH_ADIOS2)
   find_package (ADIOS2 "2.2.0")
   if (ADIOS2_FOUND)
     ADD_SUBDIRECTORY(adios2pio-a2)
   endif (ADIOS2_FOUND)
endif(WITH_ADIOS2)

