project('ncl','cpp',
        version: '2.1',
        license: 'BSD'
       )

ncl_inc_dir = include_directories('.')

subdir('ncl')

subdir('example')

subdir('test')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries: ncl,
                 version: meson.project_version(),
                 name: 'libncl',
                 filebase: 'ncl',
                 description: 'Nexus class library'
                )
