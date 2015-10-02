# Find SUNDIALS libraries.

# Search hints for SUNDIALS installation in the home directory.
find_path(SUNDIALS_DIR include/sundials/sundials_config.h HINTS ENV SUNDIALS_DIR PATHS $ENV{HOME}/sundials DOC "Sundials Directory")

set(SUNDIALS_LIB_NAMES sundials_arkode sundials_nvecserial sundials_cvode sundials_cvodes sundials_ida sundials_idas sundials_kinsol sundials_nvecparallel)

set(SUNDIALS_LIBRARIES)

IF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)

	SET(SUNDIALS_FOUND YES)
	SET(SUNDIALS_INCLUDES ${SUNDIALS_DIR})
	find_path (SUNDIALS_INCLUDE_DIR sundials_config.h HINTS "${SUNDIALS_DIR}" PATH_SUFFIXES include/sundials NO_DEFAULT_PATH)
	list(APPEND SUNDIALS_INCLUDES ${SUNDIALS_INCLUDE_DIR})

	foreach(SUNDIALS_LIB ${SUNDIALS_LIB_NAMES})
		find_library(lib${SUNDIALS_LIB} ${SUNDIALS_LIB} "${SUNDIALS_DIR}/lib")
		if(lib${SUNDIALS_LIB})
			message(STATUS "Found ${SUNDIALS_LIB}: ${lib${SUNDIALS_LIB}}")
			set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES} ${lib${SUNDIALS_LIB}})
		else(lib${SUNDIALS_LIB})
			message(STATUS "NOTFOUND: ${SUNDIALS_LIB}")
		endif(lib${SUNDIALS_LIB})
	endforeach(SUNDIALS_LIB)

ELSE(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
	SET(SUNDIALS_FOUND NO)
	message(FATAL_ERROR "Can not find SUNDIALS!")
ENDIF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS DEFAULT_MSG SUNDIALS_LIBRARIES SUNDIALS_INCLUDES)

