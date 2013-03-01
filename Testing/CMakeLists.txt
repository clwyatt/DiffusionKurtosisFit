# Testing/CMakeLists.txt

# test dkifit using small input size
ADD_TEST(test_small
	 ${EXECUTABLE_OUTPUT_PATH}/dkifit
	 -i
	 ${CMAKE_CURRENT_SOURCE_DIR}/test_small_dataset1.nrrd
	 ${CMAKE_CURRENT_SOURCE_DIR}/test_small_dataset2.nrrd
)

#set_tests_properties (test_name PROPERTIES TIMEOUT 3600)
