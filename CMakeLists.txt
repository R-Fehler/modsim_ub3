cmake_minimum_required(VERSION 3.10)
project(modsim_ub3)

find_package(GSL REQUIRED)
add_executable(gsl-ausgleichsproblem  gsl-ausgleichsproblem.c)
add_executable(ausgleichsproblem ausgleichsproblem.c)
include_directories(/usr/local/include)
target_link_libraries(gsl-ausgleichsproblem GSL::gsl GSL::gslcblas)
