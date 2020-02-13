#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  python3 \
  python3-pip \
  python3-setuptools \
  wget


if [ ! -d $install_root ]; then
  mkdir $install_root
fi

#__________________________ ART _______________________________#
cd $install_root
wget -q https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
tar xf artbinmountrainier20160605linux64tgz.tgz
cp -s art_bin_MountRainier/art_illumina .


#__________________________ python tox ________________________#
# Why six>=1.14.0?
# See https://github.com/pypa/virtualenv/issues/1551
pip3 install tox "six>=1.14.0"

