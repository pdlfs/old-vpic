#!/bin/bash
#
# Prepare a machine.conf for VPIC build
#
# for CRAY targets, let the $CRAY_CPU_TARGET take priority (since
# the user should be setting the craype module target properly).
# otherwise fall back to hostname.
#

src_dir="$1"
host=$(hostname)

rm -f $src_dir/machine.conf

if [ x$CRAY_CPU_TARGET = xhaswell ]; then
    echo "Preparing machine.conf for CRAY Haswell target"
    cp $src_dir/cray-haswell.conf $src_dir/machine.conf
elif [ x$CRAY_CPU_TARGET = xmic-knl ]; then
    echo "Preparing machine.conf for CRAY MIC-KNL target"
    cp $src_dir/cray-mic-knl.conf $src_dir/machine.conf
elif [ $(hostname | grep tt-fey) ]; then
    echo "Preparing machine.conf for Trinitite (haswell)"
    cp $src_dir/cray-haswell.conf $src_dir/machine.conf
elif [ $(hostname | grep tr-fe) ]; then
    echo "Preparing machine.conf for Trinity (haswell)"
    cp $src_dir/cray-haswell.conf $src_dir/machine.conf
elif [ $(hostname | grep ga-fe) ]; then
    echo "Preparing machine.conf for Gadget (haswell)"
    cp $src_dir/cray-haswell.conf $src_dir/machine.conf
elif [ $(hostname | grep narwhal) ] && [ `which mpirun.mpich` ]; then
    echo "Preparing machine.conf for Narwhal (MPICH)"
    cp $src_dir/narwhal-mpich.conf $src_dir/machine.conf
elif [ $(hostname | grep narwhal) ]; then
    echo "Preparing machine.conf for Narwhal (OpenMPI)"
    cp $src_dir/narwhal-openmpi.conf $src_dir/machine.conf
else
    echo "prep_conf.sh: WARNING: unknown system, using default (mpich)"
    echo "Preparing default machine.conf"
    cp $src_dir/narwhal-mpich.conf $src_dir/machine.conf
fi

exit 0
