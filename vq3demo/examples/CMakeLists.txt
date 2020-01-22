# Make sure the compiler can find include files from our library.
include_directories (${CMAKE_SOURCE_DIR}/src)

# Define our examples to compile
file(
	GLOB 
	usage_examples
	*.cpp
)

# loop over the list to compile them
foreach(f ${usage_examples})
    get_filename_component(exampleName ${f} NAME_WE) 
    add_executable            (${exampleName}                            ${f}                     ) 
    set_target_properties     (${exampleName} PROPERTIES LINKER_LANGUAGE CXX                      )
    set_target_properties     (${exampleName} PROPERTIES COMPILE_FLAGS   "${PROJECT_CFLAGS} -Wall")
    target_include_directories(${exampleName} PUBLIC                     ${OPENCV_INCLUDE_DIRS}   )
    target_link_libraries     (${exampleName}                            ${OPENCV_LIBRARIES}      )
    target_link_libraries     (${exampleName}                            ${OPENCV_LDFLAGS}        )
    target_include_directories(${exampleName} PUBLIC                     ${VQ3_INCLUDE_DIRS}      )
    target_link_libraries     (${exampleName}                            ${VQ3_LIBRARIES}         )
    target_link_libraries     (${exampleName}                            ${VQ3_LDFLAGS}           )
    target_include_directories(${exampleName} PUBLIC                     ${DEMO_INCLUDE_DIRS}     )
    target_link_libraries     (${exampleName}                            ${DEMO_LIBRARIES}        )
    target_link_libraries     (${exampleName}                            ${DEMO_LDFLAGS}          )

    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${exampleName} 
            DESTINATION bin
            RENAME ${CMAKE_PROJECT_NAME}-${exampleName}
	    COMPONENT binary)
endforeach(f)
