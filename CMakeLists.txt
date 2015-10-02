
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CRDModel)

# Add local CMake directory to CMake module path to load project-specific modules. 
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Find SUNDIALS
find_package(SUNDIALS REQUIRED)
# Add the location of the SUNDIALS headers.
include_directories(${SUNDIALS_INCLUDES})

# Unable to run without MPI.
find_package(MPI REQUIRED)
# Add the location of the MPI headers.
include_directories(${MPI_CXX_INCLUDE_PATH})

# Add variable targets.
SET(Targets
  FHNmodel
  FHNmodel_torus
  )

# Check whether we use MPI compiler scripts.
# If yes, then change the C compiler to the MPICC script.
# If not, then add the MPI include directory for MPI headers.

IF(MPI_MPICXX)
  # use MPI_MPICXX as the compiler
  SET(CMAKE_CXX_COMPILER ${MPI_MPICXX})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
ELSE(MPI_MPICXX)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICXX)

# Set-up linker flags and link libraries

# SET(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})
#IF(LAPACK_FOUND)
#  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
#ENDIF(LAPACK_FOUND)
#IF(KLU_FOUND)
#  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
#ENDIF(KLU_FOUND)
#IF(SUPERLUMT_FOUND)
#  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
#ENDIF(SUPERLUMT_FOUND)

# Add the build and install targets for each target.

FOREACH(target ${Targets})
  ADD_EXECUTABLE(${target} src/${target}.cpp)
  TARGET_LINK_LIBRARIES(${target} ${SUNDIALS_LIBS})
  IF(NOT MPI_MPICXX)
    TARGET_LINK_LIBRARIES(${target} ${MPI_LIBRARIES} ${SUNDIALS_LIBRARIES})
  ENDIF(NOT MPI_MPICXX)
ENDFOREACH(target ${Targets})


