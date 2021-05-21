# GEMSTAT

This code implements the method from the paper

- Thermodynamics-based models of transcriptional regulation by enhancers: the roles of synergistic activation, cooperative binding and short-range repression
- Author: Xin He <xinhe2@illinois.edu>

## Authors

See the file [AUTHORS](./AUTHORS) for information about additional authors and maintainers of this software.

## Documentation

See the file named [README.txt](./README.txt) for documentation about building and running the program.

## building

After checking out the repo, do:

```bash
git submodule init
git submodule update --recursive
mkdir build
cd build
cmake ..
make
```

### Maintainers

See the file named [MAINTAINER_README](./MAINTAINER_README) for more hints on how to maintain this software.

### Builds

See [Travis-Ci](https://travis-ci.com/github/UIUCSinhaLab/GEMSTAT) for continuous integration.

- 'master': ![](https://api.travis-ci.com/UIUCSinhaLab/GEMSTAT.svg?branch=master)
- 'develop': ![](https://api.travis-ci.com/UIUCSinhaLab/GEMSTAT.svg?branch=develop)
- 'bleeding-edge': ![](https://api.travis-ci.com/UIUCSinhaLab/GEMSTAT.svg?branch=bleeding-edge)
