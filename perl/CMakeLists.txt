
FILE(GLOB perlscripts RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.pl")

# FILE(GLOB perlmods RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#    "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.pl")

# Install perl scripts.
INSTALL(PROGRAMS ${perlscripts}
    DESTINATION ${InstallSuffix}
    )

# Install perl modules
INSTALL(DIRECTORY lib/
    DESTINATION ${InstallSuffix}/lib
    )

# Modify vntrseek.pl
configure_file(vntrseek.pl vntrseek.pl @ONLY)
# Install vntrseek.pl
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/vntrseek.pl DESTINATION ${InstallSuffix})

# Create symlink to vntrseek.pl in bin
if(CMAKE_HOST_UNIX)
    INSTALL(CODE "
    EXECUTE_PROCESS(COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl vntrseek
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       )
    ")
    INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/vntrseek"
        DESTINATION bin
        OPTIONAL
        )
    # INSTALL(CODE "execute_process( \
    #         COMMAND ${CMAKE_COMMAND} -E create_symlink \
    #         ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl \
    #         ${CMAKE_INSTALL_PREFIX}/bin/vntrseek   \
    #     )"
    #     ERROR_VARIABLE VNTRSEEK_SYMLINK_ERROR
    # )
    # IF(VNTRSEEK_SYMLINK_ERROR)
    #     MESSAGE(WARNING "Problem creating symlink for vntrseek command. \
    #         If you'd like the command 'vntrseek' to be in your PATH either \
    #         add the directory '${CMAKE_INSTALL_PREFIX}/${InstallSuffix}' to \
    #         your PATH or create a symbolic link to \
    #         ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl named \
    #         'vntrseek' in a directory in your path (eg, /usr/local/bin).")
    # ENDIF()
endif(CMAKE_HOST_UNIX)
