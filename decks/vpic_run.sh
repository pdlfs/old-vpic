#!/bin/bash

# TODO: infer parameters in emulab, set in cray
# TODO: deltafs server IP should not be hardcoded in script
# TODO: drop caches between baseline and deltafs
# TODO: add deltafs server code for mpich

NODES=16
CORES=64

umbrella_build_dir="$HOME/src/deltafs-umbrella/build"
output_dir="/panfs/probescratch/TableFS/vpic_test"
ip_subnet="10.92"
#output_dir="$HOME/src/vpic/decks/dump"

# Set internal variables
build_op_dir="$umbrella_build_dir/vpic-prefix/src/vpic-build"
deck_dir="$umbrella_build_dir/vpic-prefix/src/vpic/decks/trecon-part"

die () { echo "Error: $@" 1>&2; exit 1;  }

# Generate Narwhal hosts file
fqdn_suffix="`hostname | sed 's/^[^\.]*././'`"
exp_hosts="`/share/testbed/bin/emulab-listall`"

echo $exp_hosts | awk -F, '{
for (i=1; i<=NF; i++)
    print $i"'"$fqdn_suffix"'"
}' > "$output_dir/vpic.hosts" || die "failed to create vpic.hosts file"

dpoints=7
p=$CORES
while [ $dpoints -gt 0 ]
do
    /share/testbed/bin/emulab-mpirunall sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches && free -m'

    echo "Running VPIC on $CORES cores, with $(( p * p * 100 )) particles."

    # Configure VPIC experiment
    cd $deck_dir || die "cd failed"
    mv $deck_dir/config.h $deck_dir/config.bkp || die "mv failed"
    cat $deck_dir/config.bkp | \
        sed 's/^#define VPIC_FILE_PER_PARTICLE/\/\/#define VPIC_FILE_PER_PARTICLE/' | \
        sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$CORES'/' | \
        sed 's/VPIC_TOPOLOGY_Y.*/VPIC_TOPOLOGY_Y 1/' | \
        sed 's/VPIC_TOPOLOGY_Z.*/VPIC_TOPOLOGY_Z 1/' | \
        sed 's/VPIC_PARTICLE_X.*/VPIC_PARTICLE_X '$p'/' | \
        sed 's/VPIC_PARTICLE_Y.*/VPIC_PARTICLE_Y '$p'/' | \
        sed 's/VPIC_PARTICLE_Z.*/VPIC_PARTICLE_Z 1/' > $deck_dir/config.h || \
        die "config.h editing failed"

    # Compile input deck
    cd $deck_dir || die "cd failed"
    $build_op_dir/build.op ./turbulence.cxx || die "compilation failed"

    # Run VPIC experiment
    cd $output_dir || die "cd failed"
    mkdir "$output_dir/baseline_$p" || die "mkdir failed"
    cd $output_dir/baseline_$p || die "cd failed"

    # TODO: Fix for mpich, openmpi, aprun
    mpirun -np $CORES -npernode $(( CORES / NODES )) \
        --hostfile $output_dir/vpic.hosts \
        $deck_dir/turbulence.op 2>&1 | tee "$output_dir/baseline_$p.log" || \
        die "run failed"

    echo -n "Output size: " >> "$output_dir/baseline_$p.log"
    du -b $output_dir/baseline_$p | tail -1 | cut -f1 >> "$output_dir/baseline_$p.log"

    # Configure VPIC experiment
    cd $deck_dir || die "cd failed"
    mv $deck_dir/config.h $deck_dir/config.bkp || die "mv failed"
    cat $deck_dir/config.bkp | \
        sed 's/^\/\/#define VPIC_FILE_PER_PARTICLE/#define VPIC_FILE_PER_PARTICLE/' | \
        sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$CORES'/' | \
        sed 's/VPIC_TOPOLOGY_Y.*/VPIC_TOPOLOGY_Y 1/' | \
        sed 's/VPIC_TOPOLOGY_Z.*/VPIC_TOPOLOGY_Z 1/' | \
        sed 's/VPIC_PARTICLE_X.*/VPIC_PARTICLE_X '$p'/' | \
        sed 's/VPIC_PARTICLE_Y.*/VPIC_PARTICLE_Y '$p'/' | \
        sed 's/VPIC_PARTICLE_Z.*/VPIC_PARTICLE_Z 1/' > $deck_dir/config.h || \
        die "config.h editing failed"

    # Run VPIC with DeltaFS
    cd $output_dir || die "cd failed"
    mkdir "$output_dir/deltafs_$p" || die "mkdir failed"
    cd $output_dir/deltafs_$p || die "cd failed"

    # TODO: Optimize and move to top
    for mpi in openmpi mpich
    do
        which mpirun.${mpi}
        if [ $? -eq 0 ]; then
            MPI=$mpi
            break
        fi
    done

    mkdir -p $output_dir/deltafs_$p/metadata || die "deltafs metadata mkdir failed"
    mkdir -p $output_dir/deltafs_$p/data || die "deltafs data mkdir failed"

    preload_lib_path="$umbrella_build_dir/deltafs-vpic-preload-prefix/src/"\
"deltafs-vpic-preload-build/src/libdeltafs-preload.so"
    deltafs_srvr_path="$umbrella_build_dir/deltafs-prefix/src/"\
"deltafs-build/src/server/deltafs-srvr"
    deltafs_srvr_ip=`hostname -i`

    if [ x"$MPI" = xmpich ]; then
        mpirun.mpich -np $CORES --hostfile $output_dir/vpic.hosts -prepend-rank \
            -env LD_PRELOAD "$preload_lib_path" \
            -env PRELOAD_Deltafs_root "particle" \
            $deck_dir/turbulence.op 2>&1 | tee "$output_dir/deltafs_$p.log" || \
            die "mpich run failed"

    elif [ x"$MPI" = xopenmpi ]; then
        mpirun.openmpi -n 1 -tag-output \
            -x DELTAFS_MetadataSrvAddrs="$deltafs_srvr_ip:10101" \
            -x DELTAFS_FioName="posix" \
            -x DELTAFS_FioConf="root=$output_dir/deltafs_$p/data" \
            -x DELTAFS_Outputs="$output_dir/deltafs_$p/metadata" \
        	$deltafs_srvr_path | tee "$output_dir/deltafs_srvr_$p.log" || \
            die "openmpi run for deltafs server failed" &

        srvr_pid=$!

        mpirun.openmpi -np $CORES --hostfile $output_dir/vpic.hosts -tag-output \
            -x LD_PRELOAD=$preload_lib_path \
            -x PRELOAD_Deltafs_root=particle \
            -x DELTAFS_MetadataSrvAddrs="$deltafs_srvr_ip:10101" \
            -x SHUFFLE_Subnet=$ip_subnet \
            $deck_dir/turbulence.op 2>&1 | tee "$output_dir/deltafs_$p.log" || \
            (kill -KILL $srvr_pid && die "openmpi run failed")

        kill -KILL $srvr_pid
    else
        die "not running mpich or openmpi"
    fi

    echo -n "Output size: " >> "$output_dir/deltafs_$p.log"
    du -b $output_dir/deltafs_$p | tail -1 | cut -f1 >> "$output_dir/deltafs_$p.log"

    dpoints=$(( dpoints - 1 ))
    p=$(( p * 2 ))
done
