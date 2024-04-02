#!/usr/bin/env bash

set -e
set -x

if [[ "$RUNNER_OS" == "Linux" ]]; then
  # Download CUTEst and its dependencies
  mkdir "$GITHUB_WORKSPACE/cutest"
  git clone --depth 1 --branch v2.2.3 https://github.com/ralna/ARCHDefs.git "$GITHUB_WORKSPACE/cutest/archdefs"
  git clone --depth 1 --branch v2.1.5 https://github.com/ralna/SIFDecode.git "$GITHUB_WORKSPACE/cutest/sifdecode"
  git clone --depth 1 --branch v2.0.43 https://github.com/ralna/CUTEst.git "$GITHUB_WORKSPACE/cutest/cutest"
  git clone --depth 1 --branch v0.5 https://bitbucket.org/optrove/sif.git "$GITHUB_WORKSPACE/cutest/mastsif"

  # Set the environment variables
  export ARCHDEFS="$GITHUB_WORKSPACE/cutest/archdefs"
  export SIFDECODE="$GITHUB_WORKSPACE/cutest/sifdecode"
  export CUTEST="$GITHUB_WORKSPACE/cutest/cutest"
  export MASTSIF="$GITHUB_WORKSPACE/cutest/mastsif"
  export MYARCH=pc64.lnx.gfo
  {
    echo "ARCHDEFS=$ARCHDEFS"
    echo "SIFDECODE=$SIFDECODE"
    echo "CUTEST=$CUTEST"
    echo "MASTSIF=$MASTSIF"
    echo "MYARCH=$MYARCH"
  } >> "$GITHUB_ENV"

  # Build and install CUTEst
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/jfowkes/pycutest/master/.install_cutest.sh)"
elif [[ "$RUNNER_OS" == "macOS" ]]; then
  # Install gfortran
  sudo ln -fs /usr/local/bin/gfortran-12 /usr/local/bin/gfortran
  sudo ln -fs /usr/local/Cellar/gcc@12/*/lib/gcc/12/*.dylib /usr/local/lib/

  # Install CUTEst
  brew tap optimizers/cutest
  brew install cutest --without-single
  brew install mastsif

  # Set the environment variables
  for f in "archdefs" "sifdecode" "cutest" "mastsif"; do
    while IFS= read -r line; do
      echo "${line#export }" >> "$GITHUB_ENV"
    done <<< "$(cat "$(brew --prefix $f)/$f.bashrc")"
  done
fi
