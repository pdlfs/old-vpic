#!/bin/bash
#
# Prepare a machine.conf for VPIC build

src_dir="$1"
host=$(hostname)

rm -f $src_dir/machine.conf
if [ $(hostname | grep tt-fey) ]; then
    echo "Preparing machine.conf for Trinitite"
    cp $src_dir/trinitite.conf $src_dir/machine.conf
elif [ $(hostname | grep tr-fe) ]; then
    echo "Preparing machine.conf for Trinity"
    cp $src_dir/trinitite.conf $src_dir/machine.conf
elif [ $(hostname | grep ga-fe) ]; then
    echo "Preparing machine.conf for Gadget"
    cp $src_dir/trinitite.conf $src_dir/machine.conf
elif [ $(hostname | grep narwhal) ] && [ `which mpirun.mpich` ]; then
    echo "Preparing machine.conf for Narwhal (MPICH)"
    cp $src_dir/narwhal-mpich.conf $src_dir/machine.conf
elif [ $(hostname | grep narwhal) ]; then
    echo "Preparing machine.conf for Narwhal (OpenMPI)"
    cp $src_dir/narwhal-openmpi.conf $src_dir/machine.conf
else
    echo "Preparing default machine.conf"
    cp $src_dir/narwhal-openmpi.conf $src_dir/machine.conf
fi
