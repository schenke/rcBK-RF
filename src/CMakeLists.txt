add_library(
	amplitude
	STATIC
	amplitudelib.cpp
	datafile.cpp
	tools.cpp
	interpolation.cpp
)

target_link_libraries(
	amplitude
	GSL::gsl
	GSL::gslcblas
)
