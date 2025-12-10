ncl_headers = [
  'ncl.h',
  'nxsallocatematrix.h',
  'nxsassumptionsblock.h',
  'nxsblock.h',
  'nxscdiscretematrix.h',
  'nxscharactersblock.h',
  'nxscxxdiscretematrix.h',
  'nxsdatablock.h',
  'nxsdefs.h',
  'nxsdiscretedatum.h',
  'nxsdistancedatum.h',
  'nxsdistancesblock.h',
  'nxsexception.h',
  'nxsmultiformat.h',
  'nxspublicblocks.h',
  'nxsreader.h',
  'nxssetreader.h',
  'nxsstring.h',
  'nxstaxaassociationblock.h',
  'nxstaxablock.h',
  'nxstoken.h',
  'nxstreesblock.h',
  'nxsunalignedblock.h',
  'nxsutilcopy.h'
]

ncl_sources = [
  'nxsassumptionsblock.cpp',
  'nxscxxdiscretematrix.cpp',
  'nxsexception.cpp',
  'nxsreader.cpp',
  'nxstaxaassociationblock.cpp',
  'nxstreesblock.cpp',
  'nxsblock.cpp',
  'nxsdatablock.cpp',
  'nxsmultiformat.cpp',
  'nxssetreader.cpp',
  'nxstaxablock.cpp',
  'nxsunalignedblock.cpp',
  'nxscharactersblock.cpp',
  'nxsdistancesblock.cpp',
  'nxspublicblocks.cpp',
  'nxsstring.cpp',
  'nxstoken.cpp'
]

install_headers(ncl_headers, subdir: 'ncl')

ncl = both_libraries('ncl',
                     ncl_sources,
                     include_directories: ncl_inc_dir,
                     install: true)

ncl_dep = declare_dependency(
  link_with: ncl,
  include_directories: ncl_inc_dir
)

ncl_static_dep = declare_dependency(
  link_with: ncl.get_static_lib(),
  include_directories: ncl_inc_dir
)

ncl_shared_dep = declare_dependency(
  link_with: ncl.get_shared_lib(),
  include_directories: ncl_inc_dir
)
