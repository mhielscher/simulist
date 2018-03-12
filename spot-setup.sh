#!/bin/bash

if [[ $(whoami) != "root" ]]; then
    echo "Run as root."
    exit 1
fi

apt-get update
apt-get -y install python-pip python-matplotlib parallel
pip install scipy
