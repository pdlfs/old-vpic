#!/bin/bash
#
# Prepare a machine.conf for VPIC build

host=$(hostname)

rm -f machine.conf
if [ $(hostname | grep tt-fey) ]; then
    echo "Preparing machine.conf for Trinitite"
    cp trinitite.conf machine.conf
elif [ $(hostname | grep tr-fe) ]; then
    echo "Preparing machine.conf for Trinity"
    cp trinitite.conf machine.conf
elif [ $(hostname | grep narwhal) ]; then
    echo "Preparing machine.conf for Narwhal"
    cp narwhal.conf machine.conf
else
    echo "Preparing default machine.conf"
    cp narwhal.conf machine.conf
fi
