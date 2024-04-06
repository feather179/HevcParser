add_rules("mode.debug", "mode.release")

target("HevcParser")
    set_kind("binary")
    -- set_toolchains("msvc")
    set_languages("c++20")
    add_cxxflags("/EHsc")
    -- add_cxxflags("/fsanitize=address")

    add_files("main.cpp")
    add_files("common/*.cpp")
    add_files("reader/*.cpp")
    add_files("decoder/*.cpp")

    add_includedirs("common")
    add_includedirs("reader")
    add_includedirs("decoder")
