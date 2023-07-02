# Check for ZLIB
FUNCTION(FIND_ZLIB)
    MESSAGE(STATUS "Checking for ZLIB...")
    FIND_PACKAGE(ZLIB)
    IF (ZLIB_FOUND)
        MESSAGE(STATUS "ZLIB lib: ${ZLIB_LIBRARY}")
        MESSAGE(STATUS "ZLIB include: ${ZLIB_INCLUDE_DIR}")
    ELSE ()
        MESSAGE(FATAL_ERROR "ZLIB not found. Please install zlib to continue.")
    ENDIF()
ENDFUNCTION(FIND_ZLIB)

# Check Perl minimum version
FUNCTION(PERL_REQ_VER VER)
    MESSAGE(STATUS "Checking Perl version...")
    FIND_PACKAGE(Perl)
    IF (PERL_FOUND)
        IF (PERL_VERSION_STRING VERSION_GREATER ${VER} OR PERL_VERSION_STRING VERSION_EQUAL ${VER})
            MESSAGE(STATUS "Perl >= ${VER} (${PERL_VERSION_STRING})")
        ELSE ()
            MESSAGE(FATAL_ERROR "Error: Perl version must be >= ${VER}. You have ${PERL_VERSION_STRING}.")
        ENDIF ()
    ELSE (PERL_FOUND)
        MESSAGE(FATAL_ERROR "Error: ${PACKAGE_NAME} requires Perl version ${VER} in order to run.")
    ENDIF(PERL_FOUND)
ENDFUNCTION(PERL_REQ_VER)

# Check GCC minimum version
FUNCTION(GCC_REQ_VER VER)
    MESSAGE(STATUS "Checking GCC version...")
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION)
    STRING(REGEX REPLACE "(\r?\n)+$" "" GCC_VERSION "${GCC_VERSION}")

    IF(GCC_VERSION VERSION_GREATER ${VER} OR GCC_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "GCC version >= ${VER} (${GCC_VERSION})")
    ELSE()
        MESSAGE(FATAL_ERROR "Error: GCC version must be >= ${VER}. You have ${GCC_VERSION}")
    ENDIF()
ENDFUNCTION(GCC_REQ_VER)

# Check SQLITE3 minimum version
FUNCTION(SQLITE3_REQ_VER VER)
    MESSAGE(STATUS "Checking sqlite3 version...")
    find_program(SQLITE3_CMD sqlite3)
    IF(NOT SQLITE3_CMD) 
        MESSAGE(FATAL_ERROR "sqlite3 not found!")
    ENDIF()
    EXECUTE_PROCESS(COMMAND ${SQLITE3_CMD} --version | head -1
        OUTPUT_VARIABLE SQLITE3_VERSION)
    STRING(REGEX MATCH "[0-9]\\.[0-9][0-9]?\\.[0-9]" SQLITE3_VERSION "${SQLITE3_VERSION}")
    IF(SQLITE3_VERSION VERSION_GREATER ${VER} OR SQLITE3_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "sqlite3 version >= ${VER} (${SQLITE3_VERSION})")
    ELSE()
        MESSAGE(FATAL_ERROR "Error: sqlite3 version must be >= ${VER}. You have ${SQLITE3_VERSION}")
    ENDIF()
ENDFUNCTION(SQLITE3_REQ_VER)


# Check GLIBC minimum version
FUNCTION(GLIBC_REQ_VER VER)
    MESSAGE(STATUS "Checking GLIBC version...")
    EXECUTE_PROCESS(COMMAND /lib/libc.so.6
        OUTPUT_VARIABLE GLIBC_VERSION)
    STRING(REGEX MATCH "[0-9]\\.[0-9][0-9]?" GLIBC_VERSION "${GLIBC_VERSION}")
    IF(NOT (GLIBC_VERSION STREQUAL ""))
        IF(GLIBC_VERSION VERSION_GREATER ${VER} OR GLIBC_VERSION VERSION_EQUAL ${VER})
            MESSAGE(STATUS "GLIBC version >= ${VER} (${GLIBC_VERSION})")
        ELSE()
            MESSAGE(FATAL_ERROR "Error: GLIBC version must be >= ${VER}. You have ${GLIBC_VERSION}")
        ENDIF()
    ELSE()
        MESSAGE(FATAL_ERROR "Error: No GLIBC version found.  GLIBC ${VER} or greater is required.")
    ENDIF()

    IF(GLIBC_VERSION VERSION_LESS ${VER} OR GLIBC_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "Downloading legacy build of TRF...")
        SET(TRFDLName "trf${TRFVer}.legacy${ARCH}" PARENT_SCOPE)
    ENDIF()
ENDFUNCTION(GLIBC_REQ_VER)

# Check for samtools
FUNCTION(SAMTOOLS_REQ_VER VER)
    MESSAGE(STATUS "Checking samtools version...")
    find_program(SAMTOOLS_CMD samtools)
    if(NOT SAMTOOLS_CMD) 
        message(FATAL_ERROR "samtools not found!")
    endif()

    EXECUTE_PROCESS(COMMAND ${SAMTOOLS_CMD} --version | head -1
        OUTPUT_VARIABLE SAMTOOLS_VERSION)
    STRING(REGEX MATCH "[0-9]\\.[0-9][0-9]?" SAMTOOLS_VERSION "${SAMTOOLS_VERSION}")

    IF(SAMTOOLS_VERSION VERSION_GREATER ${VER} OR SAMTOOLS_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "samtools version >= ${VER} (${SAMTOOLS_VERSION})")
    ELSE()
        MESSAGE(FATAL_ERROR "Error: samtools version must be >= ${VER}. You have ${SAMTOOLS_VERSION}")
    ENDIF()
ENDFUNCTION(SAMTOOLS_REQ_VER)

# Check for bedtools
FUNCTION(BEDTOOLS_REQ_VER VER)
    MESSAGE(STATUS "Checking bedtools version...")
    find_program(BEDTOOLS_CMD bedtools)
    if(NOT BEDTOOLS_CMD) 
        message(FATAL_ERROR "Bedtools not found!")
    endif()

    EXECUTE_PROCESS(COMMAND ${BEDTOOLS_CMD} --version | head -1
        OUTPUT_VARIABLE BEDTOOLS_VERSION)
    STRING(REGEX MATCH "[0-9]\\.[0-9][0-9]?\\.[0-9]" BEDTOOLS_VERSION "${BEDTOOLS_VERSION}")

    IF(BEDTOOLS_VERSION VERSION_GREATER ${VER} OR BEDTOOLS_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "bedtools version >= ${VER} (${BEDTOOLS_VERSION})")
    ELSE()
        MESSAGE(FATAL_ERROR "Error: bedtools version must be >= ${VER}. You have ${BEDTOOLS_VERSION}")
    ENDIF()
ENDFUNCTION(BEDTOOLS_REQ_VER)

function(find_cpanm)
    find_program(CPANM_CMD cpanm)
    if(NOT CPANM_CMD)
        # Use bundled cpanm
        set(CPANM_CMD "${CMAKE_SOURCE_DIR}/perl/cpanm")
    endif()
endfunction()