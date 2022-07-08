# source files for redund.exe
add_executable(redund.exe)
target_link_libraries(redund.exe easylife sqlite3 m)
target_sources(redund.exe
    PRIVATE redund2.c # redund.c
)

install(TARGETS redund.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
