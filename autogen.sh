aclocal -I m4 --install
autoconf
automake --add-missing
#./configure --with-beta-optimize=together
bash external/buildexternal.bash
PATH=${PWD}/external/usr/bin/:${PATH} ./configure --with-beta-optimize=together --with-gsl-prefix=${PWD}/external/usr/ --with-nlopt-prefix=${PWD}/external/usr/
make clean
make
