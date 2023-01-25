set(CMAKE_CXX_STANDARD 14)
add_executable(edlib-align align.c edlib.cpp edlib.h)

install(TARGETS edlib-align
    RUNTIME DESTINATION ${InstallSuffix}
)