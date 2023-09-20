#!/usr/bin/env bash

set -e
set -x

brew tap optimizers/cutest
brew install cutest --without-single
brew install mastsif
for f in "archdefs" "mastsif" "sifdecode" "cutest"; do
  echo ". $(brew --prefix $f)/$f.bashrc" >> ~/.bashrc;
done
