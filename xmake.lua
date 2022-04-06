add_requires("libsdl")
add_packages("libsdl")
add_rules("mode.debug","mode.release")

target("test_image_load")
    set_kind("binary")
    add_files("./tests/test_image_load.cpp")
target("test_scale")
    set_kind("binary")
    add_files("./tests/test_scale.cpp")
target("test_filter")
    set_kind("binary")
    add_files("./tests/test_filter.cpp")