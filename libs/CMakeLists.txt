include(FetchContent)

FetchContent_Declare(genesis
    GIT_REPOSITORY https://github.com/lczech/genesis.git
    GIT_TAG v0.28.1)

FetchContent_GetProperties(genesis)
if(NOT genesis_POPULATED)
    FetchContent_Populate(genesis)
    add_subdirectory(${genesis_SOURCE_DIR} ${genesis_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
    GIT_TAG v3.11.2)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
