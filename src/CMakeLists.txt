add_executable(Integrated_3PCF_theory main.cpp ClassEngine.cc Engine.cc cosmology_utils.cpp integration_utils.cpp interpolation_methods.cpp power_spectrum.cpp power_spectra_2D.cpp bispectrum.cpp bispectra_2D.cpp integrated_bispectra_2D.cpp real_space_2D.cpp halo_utils.cpp)

if (UNIX AND NOT APPLE)
    target_link_libraries(Integrated_3PCF_theory stdc++fs GSL::gsl GSL::gslcblas ${CMAKE_SOURCE_DIR}/external/class_files/class_lib/libclass.a ${CMAKE_SOURCE_DIR}/external/cubature_files/cubature_lib/libcubature.a)
else()
    target_link_libraries(Integrated_3PCF_theory GSL::gsl GSL::gslcblas ${CMAKE_SOURCE_DIR}/external/class_files/class_lib/libclass.a ${CMAKE_SOURCE_DIR}/external/cubature_files/cubature_lib/libcubature.a)
endif()