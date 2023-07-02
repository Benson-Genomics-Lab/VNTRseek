INSTALL(CODE "MESSAGE(\"
Success: vntrseek is installed!

All executable code is stored at ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}.

A symbolic link has been created at ${CMAKE_INSTALL_PREFIX}/bin/vntrseek which points to the master program ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl.

Please check that ${CMAKE_INSTALL_PREFIX}/bin has been added to your PATH variable.

See the RUNNING.md instructions for extended operation instructions.\")")
