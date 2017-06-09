aclocal -I m4 --install
autoconf
automake --add-missing
#./configure --with-beta-optimize=together
bash external/buildexternal.bash
PATH=${PATH}:./external/usr/bin/ ./configure --with-beta-optimize=together --with-gsl-prefix=external/usr/ --with-nlopt-prefix=external/usr/
make
