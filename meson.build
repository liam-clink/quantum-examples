project(
    'quantum-examples',
    'cpp',
    version: '0.1',
    default_options: ['warning_level=3', 'cpp_standard=c++20'],
)


# Advice on adding precompiled library https://stackoverflow.com/questions/34663124/link-static-library-to-shared-library-or-to-a-binary
cpp = meson.get_compiler('cpp')

# (Meson requires an absolute path for find_library().)
cuba_libdir = meson.current_source_dir() + './cuba/lib/'
cuba = cpp.find_library('cuba', dirs : cuba_libdir) # ./lib/libA.lib
cuba_include = include_directories('cuba/include')



project_source_files = ['src/main.cpp']

project_dependencies = [cuba]

build_args = ['']

project_target = executable(
    meson.project_name(),
    project_source_files,
    include_directories: [cuba_include],
    dependencies: project_dependencies,
    install: true,
    c_args: build_args,
)

test('basic', project_target)
