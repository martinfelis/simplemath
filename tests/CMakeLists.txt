PROJECT (SIMPLEMATH_TESTS CXX)

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Enable proper C++11 flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Needed for UnitTest++
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ )

SET ( TESTS_SRCS
	test_main.cc
	#	TestSimpleMathFixed.cc
	# TestSimpleMathDynamic.cc
	# TestSimpleMathMixed.cc
	# TestSimpleMathQR.cc
	TestSimpleMathV2.cc
	# TestSimpleMathCholesky.cc
	# TestSimpleMathBlock.cc
	# TestSimpleMathQuaternion.cc
	# TestSimpleMathCommaInitializer.cc
	)

FIND_PACKAGE (UnitTest++ REQUIRED)

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/../include )

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
  LINKER_LANGUAGE CXX
)

IF ( UNITTEST++_FOUND )
	ADD_EXECUTABLE ( simplemath_tests ${TESTS_SRCS} )
	ADD_EXECUTABLE ( output_test output_test.cc )

  INCLUDE_DIRECTORIES ( ${UNITTEST++_INCLUDE_DIR} )

	SET_TARGET_PROPERTIES ( simplemath_tests PROPERTIES
		LINKER_LANGUAGE CXX
		OUTPUT_NAME runtests
		)

	TARGET_LINK_LIBRARIES ( simplemath_tests
			${UNITTEST++_LIBRARY}
		)
		
	OPTION (RUN_AUTOMATIC_TESTS "Perform automatic tests after compilation?" OFF)

	IF (RUN_AUTOMATIC_TESTS)
	ADD_CUSTOM_COMMAND (TARGET runtests
		POST_BUILD
		COMMAND ./runtests
		COMMENT "Running automated tests..."
		)
	ENDIF (RUN_AUTOMATIC_TESTS)

ENDIF ( UNITTEST++_FOUND )

