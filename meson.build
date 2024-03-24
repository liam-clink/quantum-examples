project(
    'quantum-examples',
    'cpp',
    version: '0.1',
    default_options: ['warning_level=3', 'cpp_standard=c++20'],
)


# Advice on adding precompiled library https://stackoverflow.com/questions/34663124/link-static-library-to-shared-library-or-to-a-binary
cpp = meson.get_compiler('cpp')

# (Meson requires an absolute path for find_library().)
libdir = meson.current_source_dir() + '/third-party/build/lib/'
message(libdir)
# There should not be extension on library
gsl = cpp.find_library('gsl', dirs : libdir, static: true)
include = include_directories('third-party/build/include')



project_source_files = ['src/main.cpp']

project_dependencies = [gsl]

build_args = ['']

project_target = executable(
    meson.project_name(),
    project_source_files,
    include_directories: [include],
    dependencies: project_dependencies,
    install: true,
    c_args: build_args,
)

test('basic', project_target)
