#!/bin/bash

NODES=4
CORES=16

build_dir="/users/gamvrosi/src/vpic"
deck_dir="/users/gamvrosi/src/vpic/decks/trecon-part"
output_dir="/users/gamvrosi/src/vpic/decks/dump"

die () { echo "Error: $@" 1>&2; exit 1;  }

# Generate Narwhal hosts file
fqdn_suffix="`hostname | sed 's/^[^\.]*././'`"
exp_hosts="`/share/testbed/bin/emulab-listall`"

echo $exp_hosts | awk -F, '{
    for (i=1; i<=NF; i++)
        print $i"'"$fqdn_suffix"'"
    }' > $output_dir/vpic.hosts || \
    die "failed to create vpic.hosts file"

dpoints=10
p=$CORES
while [ $dpoints -gt 0 ]
do
    echo "Running VPIC on $CORES cores, with $(( p * p * 200 )) particles."

    # Configure VPIC experiment
    cd $deck_dir || die "cd failed"
    mv $deck_dir/config.h $deck_dir/config.bkp || die "mv failed"
    cat $deck_dir/config.bkp | \
        sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$CORES'/' | \
        sed 's/VPIC_TOPOLOGY_Y.*/VPIC_TOPOLOGY_Y 1/' | \
        sed 's/VPIC_TOPOLOGY_Z.*/VPIC_TOPOLOGY_Z 1/' | \
        sed 's/VPIC_PARTICLE_X.*/VPIC_PARTICLE_X '$p'/' | \
        sed 's/VPIC_PARTICLE_Y.*/VPIC_PARTICLE_Y '$p'/' | \
        sed 's/VPIC_PARTICLE_Z.*/VPIC_PARTICLE_Z 1/' | > $deck_dir/config.h || \
        die "config.h editing failed"

    # Compile input deck
    $build_dir/build.op $deck_dir/turbulence.cxx || die "compilation failed"

    # Run VPIC experiment
    mkdir "$output_dir/run_$p" || die "mkdir failed"
    cd $output_dir/run_$p || die "cd failed"

    mpirun -np $CORES -npernode $(( CORES / NODES )) \
        --hostfile $output_dir/vpic.hosts $deck_dir/turbulence.op || \
        die "run failed"

    dpoints=$(( dpoints - 1 ))
done
