#!/bin/bash
#
#MSUB -N deltafs-test
#MSUB -l walltime=1:00:00
#MSUB -l nodes=5:haswell
#MSUB -o /users/$USER/joblogs/deltafs-test-$MOAB_JOBID.out
#MSUB -j oe
##MSUB -V
##MSUB -m b
##MSUB -m $USER@lanl.gov


######################
# Tunable parameters #
######################

# Reminder: ensure that the number of nodes here matches the preamble above.
# You'll get nicer particle counts if this is 2^n + 1 (see note below).
nodes=5

# N.B.: This is the number of cores allocated to VPIC. One node is reserved for
# the DeltaFS server, and the rest are used by VPIC. Because the number of VPIC
# particles across runs is a multiple of the number of cores allocated to it,
# you'll get nicer particle counts across iterations if this is a power of two.
cores=$(( (nodes - 1) * 8 ))

# Paths
umbrella_build_dir="$HOME/src/deltafs-umbrella/build"
output_dir="$HOME/src/vpic/decks/dump"


# DeltaFS vars
ip_subnet="10.92" # Subnet from which to pick IPs for Mercury address generation


###############
# Core script #
###############

message () { echo "$@" | tee $logfile; }
die () { message "Error $@"; exit 1; }

# Clear VPIC node caches
clear_caches() {
    message "Clearing node caches..."
    aprun -L $vpic_nodes -n $cores -N $((cores / (nodes-1))) sudo sh -c \
        'echo 3 > /proc/sys/vm/drop_caches'
}

# Internal variables (no need to touch)
build_op_dir="$umbrella_build_dir/vpic-prefix/src/vpic-build"
deck_dir="$umbrella_build_dir/vpic-prefix/src/vpic/decks/trecon-part"
logfile=""

# Configure config.h
# @1 in {"file-per-process", "file-per-particle"}
# @2 particles
build_deck() {
    p=$2

    cd $deck_dir || die "cd failed"
    mv $deck_dir/config.h $deck_dir/config.bkp || die "mv failed"

    case $1 in
    "file-per-process")
        cat $deck_dir/config.bkp | \
            sed 's/^#define VPIC_FILE_PER_PARTICLE/\/\/#define VPIC_FILE_PER_PARTICLE/' | \
            sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$cores'/' | \
            sed 's/VPIC_TOPOLOGY_Y.*/VPIC_TOPOLOGY_Y 1/' | \
            sed 's/VPIC_TOPOLOGY_Z.*/VPIC_TOPOLOGY_Z 1/' | \
            sed 's/VPIC_PARTICLE_X.*/VPIC_PARTICLE_X '$p'/' | \
            sed 's/VPIC_PARTICLE_Y.*/VPIC_PARTICLE_Y '$p'/' | \
            sed 's/VPIC_PARTICLE_Z.*/VPIC_PARTICLE_Z 1/' > $deck_dir/config.h || \
            die "config.h editing failed"
        ;;
    "file-per-particle")
        cat $deck_dir/config.bkp | \
            sed 's/^\/\/#define VPIC_FILE_PER_PARTICLE/#define VPIC_FILE_PER_PARTICLE/' | \
            sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$cores'/' | \
            sed 's/VPIC_TOPOLOGY_Y.*/VPIC_TOPOLOGY_Y 1/' | \
            sed 's/VPIC_TOPOLOGY_Z.*/VPIC_TOPOLOGY_Z 1/' | \
            sed 's/VPIC_PARTICLE_X.*/VPIC_PARTICLE_X '$p'/' | \
            sed 's/VPIC_PARTICLE_Y.*/VPIC_PARTICLE_Y '$p'/' | \
            sed 's/VPIC_PARTICLE_Z.*/VPIC_PARTICLE_Z 1/' > $deck_dir/config.h || \
            die "config.h editing failed"
        ;;
    *)
        die "build_deck: VPIC mode not supported"
        ;;
    esac

    # Compile input deck
    cd $deck_dir || die "cd failed"
    $build_op_dir/build.op ./turbulence.cxx || die "compilation failed"
}

# Run VPIC
# @1 in {"baseline", "deltafs"}
# @2 number of particles
do_run() {
    runtype=$1
    p=$2

    cd $output_dir || die "cd failed"
    mkdir "$output_dir/${runtype}_$p" || die "mkdir failed"
    cd $output_dir/${runtype}_$p || die "cd failed"

    # Define logfile before calling message()
    logfile="$output_dir/${runtype}_$p.log"

    clear_caches

    message ""
    message "=========================================================="
    message "Running VPIC ($runtype) with $(( p * p * 100 )) particles."
    message "=========================================================="
    message ""

    case $runtype in
    "baseline")
        aprun -L $vpic_nodes -n $cores -N $((cores / (nodes-1))) \
            "$deck_dir/turbulence.op" 2>&1 | tee $logfile
        if [ $? -ne 0 ]; then
            die "baseline: aprun failed"
        fi

        echo -n "Output size: " >> $logfile
        du -b $output_dir/baseline_$p | tail -1 | cut -f1 >> $logfile
        ;;

    "deltafs")
        mkdir -p $output_dir/deltafs_$p/metadata || \
            die "deltafs metadata mkdir failed"
        mkdir -p $output_dir/deltafs_$p/data || \
            die "deltafs data mkdir failed"

        preload_lib_path="$umbrella_build_dir/deltafs-vpic-preload-prefix/src/"\
"deltafs-vpic-preload-build/src/libdeltafs-preload.so"
        deltafs_srvr_path="$umbrella_build_dir/deltafs-prefix/src/"\
"deltafs-build/src/server/deltafs-srvr"

#        aprun -L $deltafs_node -n 1 -N 1 \
#            -e DELTAFS_MetadataSrvAddrs="$deltafs_srvr_ip:10101" \
#            -e DELTAFS_FioName="posix" \
#            -e DELTAFS_FioConf="root=$output_dir/deltafs_$p/data" \
#            -e DELTAFS_Outputs="$output_dir/deltafs_$p/metadata" \
#            $deltafs_srvr_path 2>&1 | tee $logfile
#        if [ $? -ne 0 ]; then
#            die "deltafs server: aprun failed"
#        fi
#
#        srvr_pid=$!

        aprun -L $vpic_nodes -n $cores -N $((cores / (nodes-1))) \
            -e LD_PRELOAD="$preload_lib_path" \
            -e PRELOAD_Deltafs_root="particle" \
            -e PRELOAD_Local_root="${output_dir}" \
            -e PRELOAD_Bypass_deltafs_namespace=1 \
            -e PRELOAD_Enable_verbose_error=1 \
            -e SHUFFLE_Virtual_factor=1024 \
            -e SHUFFLE_Mercury_proto="bmi+tcp" \
            -e SHUFFLE_Subnet="$ip_subnet" \
            "$deck_dir/turbulence.op" 2>&1 | tee $logfile
#            -e DELTAFS_MetadataSrvAddrs="$deltafs_srvr_ip:10101" \
        if [ $? -ne 0 ]; then
#            kill -KILL $srvr_pid
            die "deltafs: mpirun failed"
        fi

#        kill -KILL $srvr_pid

        echo -n "Output size: " >> $logfile
        du -b $output_dir/deltafs_$p | tail -1 | cut -f1 >> $logfile
        ;;
    esac
}

# Generate aprun -L comma-separated node lists
deltafs_node=$(cat $PBS_NODEFILE | uniq | sort | head -n 1 | tr '\n' ',')
vpic_nodes=$(cat $PBS_NODEFILE | uniq | sort -r | head -n $((nodes-1)) | \
             tr '\n' ',')

# Get Deltafs server node IP
aprun -L $deltafs_node -n 1 -N 1 hostname -i > $output_dir/hostname-i.txt
deltafs_srvr_ip=$(cat $output_dir/hostname-i.txt | head -1)
echo "DeltaFS Server IP: $deltafs_srvr_ip"

dpoints=7
parts=$cores
while [ $dpoints -gt 0 ]
do
    build_deck "file-per-process" $parts
    do_run "baseline" $parts

    build_deck "file-per-particle" $parts
    do_run "deltafs" $parts

    dpoints=$(( dpoints - 1 ))
    parts=$(( parts * 2 ))
done
