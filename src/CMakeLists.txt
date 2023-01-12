# Install euclidean distance binary file
INSTALL(FILES eucledian.dst
    DESTINATION ${InstallSuffix}
    )

ADD_SUBDIRECTORY(libs)

#vntr
ADD_SUBDIRECTORY(pcr_dup)
ADD_SUBDIRECTORY(redund_code)
ADD_SUBDIRECTORY(newflankalign)
ADD_SUBDIRECTORY(join_clusters)
ADD_SUBDIRECTORY(trf2proclu-ngs)
ADD_SUBDIRECTORY(psearch1.91)
ADD_SUBDIRECTORY(newrefflankalign)
ADD_SUBDIRECTORY(edlib)

# seqtk
# update the hash in GIT_TAG with each new release
include(ExternalProject)
ExternalProject_Add(seqtk
  GIT_REPOSITORY    https://github.com/lh3/seqtk.git
  GIT_TAG           fe1cc9e8d926d40a88f6ff368c07466d33c3665d
  SOURCE_DIR        ${CMAKE_SOURCE_DIR}/src/seqtk
  BUILD_IN_SOURCE   1
  CONFIGURE_COMMAND ""
  INSTALL_DIR       "${CMAKE_INSTALL_PREFIX}/${InstallSuffix}"
  INSTALL_COMMAND   ""
)
ExternalProject_Get_Property(seqtk source_dir)
install(PROGRAMS ${source_dir}/seqtk
	COMPONENT RUNTIME
	DESTINATION ${InstallSuffix}
)