add_executable(image-box-filter image_box_filter.cpp)
add_executable(image-transformer image_transformer.cpp)

set(app_targets image-box-filter image-transformer)

foreach(target ${app_targets})
  target_link_libraries(${target} PRIVATE ceng391)
  target_compile_features(${target} PRIVATE cxx_std_17)
  target_compile_options(${target} PRIVATE -Wall -Wextra)
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
endforeach()


