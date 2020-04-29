Docs from: https://majiq.biociphers.org/app_download/verify.php?email=weisburd@broadinstitute.org&hash=5e6fd3101d42b7fec52180e97da37bce

Before installation
MAJIQ instlalation requires the following lib/software to be installed in your system.

C++11 compiler with openMP. GCC includes that by default, but clang can be updated to include those (Clang/OMP). MAJIQ/VOILA has been tested to work with GNU GCC>=7.2, RedHat GCC>=4.7.2.
HTSlib library. This is a C library for reading/writing high-throughput sequencing data developed by Samtools organisation. MAJIQ installation assumes the library and its header files are present in the Unix default locations (/usr/local/lib, /usr/local/include). If that is not the case the appropiate locations can be specified setting the following environment variables.
$ export HTSLIB_LIBRARY_DIR=/path/to/htslib/lib
$ export HTSLIB_INCLUDE_DIR=/path/to/htslib/include
Installation
To download and install MAJIQ/Voila run the following commands:

$ python3 -m venv env
$ source env/bin/activate
$ pip install pip -U
$ pip install wheel setuptools -U
$ pip install cython numpy GitPython -U
$ pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
If there is an error during install please verify that you're installing these packages with python 3 and the installed version of pip is current. Also, check if you have git, a c compiler (gcc, clang, etc.), and zlib installed. If everything appears to be correct or you're in need of some help contact us at admin@biociphers.org.


-------
MacOSX installation commands:

sudo python3 -m pip install pip -U
sudo python3 -m pip install wheel setuptools -U
sudo python3 -m pip install cython numpy GitPython -U
CC=/usr/local/bin/gcc-9 CXX=/usr/local/bin/g++-9 python3 -m pip install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
