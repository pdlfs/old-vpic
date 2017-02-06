#!/bin/bash

# Tunable parameters: tweak to your liking
# Note: Use a number of cores that is a power of two plus 1.
#       The number of particles is a multiple of the number of cores,
#       but 1 core is reserved for DeltaFS server, so you want to be
#       left with a power of two to get better particle numbers.
CORES=5
umbrella_build_dir="$HOME/src/deltafs-umbrella/build"
output_dir="/panfs/probescratch/TableFS/vpic_test"
ip_subnet="10.92"

# Internal variables (no need to touch)
build_op_dir="$umbrella_build_dir/vpic-prefix/src/vpic-build"
deck_dir="$umbrella_build_dir/vpic-prefix/src/vpic/decks/trecon-part"
dpoints=7
logfile=""

message () { echo "$@" | tee $logfile; }
die () { message "Error $@"; exit 1; }

# Generate mpirun hostfile on Narwhal
gen_hosts() {
    message "Generating hostfile..."

    fqdn_suffix="`hostname | sed 's/^[^\.]*././'`"
    exp_hosts="`/share/testbed/bin/emulab-listall`"

    echo $exp_hosts | \
        awk -F, '{ for (i=1; i<=NF; i++) { print $i "'"$fqdn_suffix"'"}}' > \
        "$output_dir/vpic.hosts" || die "failed to create vpic.hosts file"
}

# Clear node caches on Narwhal
clear_caches() {
    message "Clearing caches..."

    /share/testbed/bin/emulab-mpirunall sudo sh -c \
        'echo 3 > /proc/sys/vm/drop_caches'
}

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
            sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$((CORES-1))'/' | \
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
            sed 's/VPIC_TOPOLOGY_X.*/VPIC_TOPOLOGY_X '$((CORES-1))'/' | \
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

# Run MPICH or OpenMPI mpirun command
# Arguments:
# @1 number of processes
# @2 array of env vars: ("name1", "val1", "name2", ... )
# @3 other options (to be appended)
# @4 executable
# @5 outfile
do_mpirun() {
    procs=$1
    declare -a envs=("${!2}")
    opts="$3"
    exe="$4"
    outfile="$5"

    if [ `which mpirun.mpich` ]; then
        if [ ${#envs[@]} -gt 0 ]; then
            envstr=`printf -- "-env %s \"%s\" " ${envs[@]}`
        else
            envstr=""
        fi

        mpirun.mpich -np $procs --hostfile $output_dir/vpic.hosts \
            $envstr -prepend-rank $opts $exe 2>&1 | tee "$outfile"

    elif [ `which mpirun.openmpi` ]; then
        if [ ${#envs[@]} -gt 0 ]; then
            envstr=`printf -- "-x %s=%s " ${envs[@]}`
        else
            envstr=""
        fi

        mpirun.openmpi -np $procs --hostfile $output_dir/vpic.hosts \
            $envstr -tag-output $opts $exe 2>&1 | tee "$outfile"

    else
        die "not running mpich or openmpi"
    fi
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
        vars=()

        do_mpirun $((CORES - 1)) vars[@] "" "$deck_dir/turbulence.op" $logfile
        if [ $? -ne 0 ]; then
            die "baseline: mpirun failed"
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
        deltafs_srvr_ip=`hostname -i`

        vars=("DELTAFS_MetadataSrvAddrs" "$deltafs_srvr_ip:10101"
              "DELTAFS_FioName" "posix"
              "DELTAFS_FioConf" "root=$output_dir/deltafs_$p/data"
              "DELTAFS_Outputs" "$output_dir/deltafs_$p/metadata")

        do_mpirun 1 vars[@] "" "$deltafs_srvr_path" $logfile
        if [ $? -ne 0 ]; then
            die "deltafs server: mpirun failed"
        fi

        srvr_pid=$!

        vars=("LD_PRELOAD" "$preload_lib_path"
              "PRELOAD_Deltafs_root" "particle"
              "DELTAFS_MetadataSrvAddrs" "$deltafs_srvr_ip:10101"
              "SHUFFLE_Subnet" "$ip_subnet")

        do_mpirun $((CORES - 1)) vars[@] "" "$deck_dir/turbulence.op" $logfile
        if [ $? -ne 0 ]; then
            kill -KILL $srvr_pid
            die "deltafs: mpirun failed"
        fi

        kill -KILL $srvr_pid

        echo -n "Output size: " >> $logfile
        du -b $output_dir/deltafs_$p | tail -1 | cut -f1 >> $logfile
        ;;
    esac
}

gen_hosts

parts=$((CORES - 1))
while [ $dpoints -gt 0 ]
do
    build_deck "file-per-process" $parts
    do_run "baseline" $parts

    build_deck "file-per-particle" $parts
    do_run "deltafs" $parts

    dpoints=$(( dpoints - 1 ))
    parts=$(( parts * 2 ))
done
