autoreconf || exit
CPPFLAGS="-D__extern_always_inline=inline" CC=clang CXX=clang++ ./configure --prefix=$PWD/installed --with-constfuncs || exit
make -j4 || exit 
make check || exit
make install || exit
make installcheck || exit
