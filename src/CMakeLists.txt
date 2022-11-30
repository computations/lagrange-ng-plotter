add_executable(lagrange-ng-plotter main.cpp)

target_link_libraries(lagrange-ng-plotter PUBLIC genesis_lib_shared nlohmann_json::nlohmann_json)
