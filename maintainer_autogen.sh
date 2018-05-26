git submodule init
git submodule update

if [ -x "$(command -v glibtoolize)" ]
then
	glibtoolize
else
	libtoolize
fi

aclocal -I m4 --install
autoconf
automake --add-missing
#Use this line for normal compilation
./configure --with-beta-optimize=together --enable-gtest
#use these lines instead to download and compile GSL and NLOPT, then statically link them into gemstat
#bash external/buildexternal.bash
#PATH=${PWD}/external/usr/bin/:${PATH} ./configure --with-beta-optimize=together --with-gsl-prefix=${PWD}/external/usr/ --with-nlopt-prefix=${PWD}/external/usr/
make clean
make
make check
