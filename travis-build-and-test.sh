autoreconf || exit
CPPFLAGS="-D__extern_always_inline=inline" CC=clang CXX=clang++ CXXFLAGS="-Weverything -pedantic -Wsign-conversion -Wpadded -Wconversion -Wweak-vtables -Wdocumentation-unknown-command -Wunused-exception-parameter -std=c++98" ./configure --prefix=$PWD/installed --with-constfuncs || exit
make -j4 || exit 
make check || exit
make install || exit
make installcheck || exit
