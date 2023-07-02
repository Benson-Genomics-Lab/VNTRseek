SET(PACKAGE VNTRseek)
SET(PACKAGE_NAME VNTRseek)

cmake_minimum_required(VERSION 3.2)
# Set current VNTRseek version. Not user configurable
SET(VNTRVer "2.0.3")
project(VNTRseek VERSION ${VNTRVer} LANGUAGES C CXX)
# project(VNTRseek VERSION ${VNTRVer} LANGUAGES C)

SET(InstallSuffix "vntrseek${VNTRVer}")

# Set current TRF version. Not user configurable
SET(TRFVer "409")

## Set install directory
#IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    SET(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix (eg, /opt, the default). All non-web components of the pipeline will be installed in a directory under this, including binaries and support scripts or files." FORCE)
#ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Set the build type -- Does nothing at the moment
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Flags for the different builds. So far, we only have a definition for flags
# in release mode
SET(CMAKE_C_FLAGS_RELEASE "-O2")
# Dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

## Detect arch and OS
INCLUDE("DetermineArch.cmake")
DETERMINE_ARCH(ARCH)

IF (ARCH STREQUAL "UNSUPPORTED")
    MESSAGE(FATAL_ERROR "Your current OS and architecture (${CMAKE_SYSTEM_NAME}, ${CMAKE_SYSTEM_PROCESSOR}) are unsupported. Installation cannot continue.")
ELSEIF (ARCH STREQUAL "i386")
    MESSAGE(FATAL_ERROR "32-bit architecture is not supported. Installation cannot continue.")
ENDIF()

## Check for required installed programs and versions
SET(TRFDLName "trf${TRFVer}.${ARCH}")
SET(GLIB2_LEGACY_VERSION "2.14")
INCLUDE("FindDeps.cmake")
FIND_ZLIB()
GCC_REQ_VER(4.1.2)
GLIBC_REQ_VER(2.14)
PERL_REQ_VER(5.24.0)
SQLITE3_REQ_VER(3.37.0)
SAMTOOLS_REQ_VER(1.8)
find_cpanm()

message(STATUS "Your processor is ${CMAKE_SYSTEM_PROCESSOR} and you are running ${CMAKE_SYSTEM_NAME}. This means we'll download ${TRFDLName}")
# May need to change if other platforms are supported
#SET(TRFBin "${TRFDLName}.exe")
SET(TRFBin "trf${TRFVer}-ngs.linux.exe")

## Check for required Perl modules
# find_package(PerlModules COMPONENTS DBI DBD::SQLite Try::Tiny)
# IF(NOT PERLMODULES_FOUND)
#     message(FATAL_ERROR "Not all required perl modules were found on your system")
# ENDIF(NOT PERLMODULES_FOUND)

# Build files in src
add_subdirectory(src)
add_subdirectory(perl)
# add_subdirectory(scripts)
#ADD_SUBDIRECTORY(etc)

# Install required perl modules locally
install(CODE "
    message(STATUS \"Installing required perl modules (using ${CPANM_CMD})\")
    execute_process(COMMAND ${CPANM_CMD} -L local --from \"${CMAKE_SOURCE_DIR}/vendor/cache\" --installdeps --notest --quiet .
        WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}
        )
")

# Post install message
add_subdirectory(postinst)

# Install config and other files
install(CODE "
    if (EXISTS \"${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/defaults.vs.cnf\")
        MESSAGE(STATUS \"defaults.vs.cnf file exists. Will install new version as 'defaults.vs.cnf.new'\")
        file(INSTALL \"${CMAKE_SOURCE_DIR}/defaults.vs.cnf\"
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}
            RENAME \"defaults.vs.cnf.new\"
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        )
    else()
        file(INSTALL \"${CMAKE_SOURCE_DIR}/defaults.vs.cnf\"
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        )
    endif(EXISTS \"${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/defaults.vs.cnf\")
")
SET(MISC_FILES
    README.md
    RUNNING.md
    RESULTS.md
    HISTORY.md
    COPYING
    REFERENCE_FORMATS
    reduce_db.sql
    sqlite_schema.sql
    cpanfile
    cpanfile.snapshot
    )
install(FILES ${MISC_FILES}
    DESTINATION ${InstallSuffix}
    )

# Download and install trf
file(DOWNLOAD https://github.com/Benson-Genomics-Lab/TRF/releases/latest/download/${TRFDLName}
    "${CMAKE_CURRENT_BINARY_DIR}/${TRFBin}"
    TIMEOUT 30
    STATUS TRFdlStat
    SHOW_PROGRESS)
list(GET TRFdlStat 0 TRFdlErrNo)
IF (TRFdlErrNo)
    message(WARNING "Error encountered downloading TRF. Download it manually from: https://github.com/Benson-Genomics-Lab/TRF/releases/latest.  Save it as ${CMAKE_CURRENT_BINARY_DIR}/${TRFBin}.")
ELSE (TRFdlErrNo)
    install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${TRFBin}" DESTINATION ${InstallSuffix})
ENDIF (TRFdlErrNo)
