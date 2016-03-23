
# project (project_name)
project (Pendule_1D)

cmake_minimum_required(VERSION 2.6)

file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}" d)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${d}/CMake")
find_package(Imagine)

#ImagineAddExecutable(folder_name main.cpp example.cpp example.h)
 ImagineAddExecutable(Pendule pendule.cpp matrice.cpp matrice.hpp)

#ImagineUseModules(folder_name module_name)
 ImagineUseModules(Pendule Graphics)


#add_executable(folder_name main.cpp example.cpp example.h)