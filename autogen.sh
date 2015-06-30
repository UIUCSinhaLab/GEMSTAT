aclocal -I m4 --install
autoconf
automake --add-missing
./configure --with-beta-optimize=together
make
