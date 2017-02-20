#!/bin/bash -x
#
#MSUB -N deltafs-test
#MSUB -l walltime=1:00:00
#MSUB -l nodes=5:haswell
#MSUB -o /users/$USER/joblogs/deltafs-test-$MOAB_JOBID.out
#MSUB -j oe
##MSUB -V
##MSUB -m b
##MSUB -m $USER@lanl.gov

# TODO:
# - Convert node lists to ranges on CRAY

######################
# Tunable parameters #
######################

# Notes:
# ------
#
# nodes: Use an odd number of nodes. The number of particles is a multiple of
# the number of cores, but 1 node is reserved for DeltaFS server, so you want
# to be left with power-of-2 cores to get better particle numbers.
#
# bbos_buddies: An additional number of node dedicated for burst buffer
# communication. Should be set to the same number of nodes as the burst buffer
# nodes.

# Node topology
cores_per_node=4
nodes=3
bbos_buddies=2

# Paths
umbrella_build_dir="$HOME/src/deltafs-umbrella/build"
output_dir="$HOME/src/vpic/decks/dump"

# DeltaFS config
ip_subnet="10.92"


###############
# Core script #
###############

cores=$(((nodes-1) * cores_per_node))
build_op_dir="$umbrella_build_dir/vpic-prefix/src/vpic-build"
deck_dir="$umbrella_build_dir/vpic-prefix/src/vpic/decks/trecon-part"
dpoints=1
logfile=""

bb_sst_size=2           # BBOS SST table size in MB
bb_log_size=1024        # BBOS max per-core log size in MB

bb_clients=$cores
bb_client_cfg="$umbrella_build_dir/deltafs-bb-prefix/src/deltafs-bb/config/narwhal_8_client.conf"
bb_client="$umbrella_build_dir/deltafs-bb-prefix/src/deltafs-bb-build/src/bbos_client"

bb_servers=$bbos_buddies
bb_server_cfg="$umbrella_build_dir/deltafs-bb-prefix/src/deltafs-bb/config/narwhal_2_server.conf"
bb_server="$umbrella_build_dir/deltafs-bb-prefix/src/deltafs-bb-build/src/bbos_server"

message () { echo "$@" | tee -a $logfile; }
die () { message "Error $@"; exit 1; }

# Generate host lists
gen_hosts() {
    message "Generating host lists..."

    if [ `which aprun` ]; then
        # Generate host lists on CRAY and store them on disk
        cat $PBS_NODEFILE | uniq | sort | head -n 1 | \
            tr '\n' ',' > $output_dir/deltafs.hosts || \
            die "failed to create deltafs.hosts file"
        cat $PBS_NODEFILE | uniq | sort | head -n $nodes | tail -n $((nodes-1)) | \
            tr '\n' ',' > $output_dir/vpic.hosts || \
            die "failed to create vpic.hosts file"
        cat $PBS_NODEFILE | uniq | sort | tail -n $bbos_buddies | \
            tr '\n' ',' > $output_dir/bbos.hosts || \
            die "failed to create bbos.hosts file"

    else
        # Generate host lists on Emulab and store them on disk
        fqdn_suffix="`hostname | sed 's/^[^\.]*././'`"
        exp_hosts="`/share/testbed/bin/emulab-listall | tr ',' '\n'`"

        echo $exp_hosts | head -n 1 | \
            tr '\n' ',' > $output_dir/deltafs.hosts || \
            die "failed to create deltafs.hosts file"
        echo $exp_hosts | head -n $nodes | tail -n $((nodes-1)) | \
            tr '\n' ',' > $output_dir/vpic.hosts || \
            die "failed to create vpic.hosts file"
        echo $exp_hosts | tail -n $bbos_buddies | \
            tr '\n' ',' > $output_dir/bbos.hosts || \
            die "failed to create bbos.hosts file"
    fi

    # Populate host list variables
    deltafs_nodes=$(cat $output_dir/deltafs.hosts)
    vpic_nodes=$(cat $output_dir/vpic.hosts)
    bbos_nodes=$(cat $output_dir/bbos.hosts)
}

# Clear node caches on Narwhal
clear_caches() {
    message "Clearing node caches..."

    if [ `which aprun` ]; then
        aprun -L $vpic_nodes -n $cores -N $((nodes - 1)) sudo sh -c \
            'echo 3 > /proc/sys/vm/drop_caches'
    else 
        /share/testbed/bin/emulab-mpirunall sudo sh -c \
            'echo 3 > /proc/sys/vm/drop_caches'
    fi
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

# Run CRAY MPICH, ANL MPICH, or OpenMPI run command
# Arguments:
# @1 number of processes
# @2 array of env vars: ("name1", "val1", "name2", ... )
# @3 host list (comma-separated)
# @4 executable (and any options that don't fit elsewhere)
# @5 outfile (used to log output)
do_mpirun() {
    procs=$1
    declare -a envs=("${!2}")
    hosts="$3"
    exe="$4"
    outfile="$5"

    if [ `which aprun` ]; then
        # This is likely a CRAY machine. Deploy an aprun job.
        if [ ${#envs[@]} -gt 0 ]; then
            envstr=`printf -- "-e %s=\"%s\" " ${envs[@]}`
        else
            envstr=""
        fi

        aprun -L $hosts -n $procs $envstr $exe 2>&1 | \
            tee -a $outfile

    elif [ `which mpirun.mpich` ]; then
        if [ ${#envs[@]} -gt 0 ]; then
            envstr=`printf -- "-env %s \"%s\" " ${envs[@]}`
        else
            envstr=""
        fi

        mpirun.mpich -np $procs --host $hosts $envstr -prepend-rank $exe 2>&1 | \
            tee -a $outfile

    elif [ `which mpirun.openmpi` ]; then
        if [ ${#envs[@]} -gt 0 ]; then
            envstr=`printf -- "-x %s=%s " ${envs[@]}`
        else
            envstr=""
        fi

        mpirun.openmpi -np $procs --host $hosts $envstr -tag-output $exe 2>&1 | \
            tee -a "$outfile"

    else
        die "could not find a supported mpirun or aprun command"
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

        do_mpirun $cores vars[@] "$vpic_nodes" "$deck_dir/turbulence.op" $logfile
        if [ $? -ne 0 ]; then
            die "baseline: mpirun failed"
        fi

        echo -n "Output size: " >> $logfile
        du -b $output_dir/baseline_$p | tail -1 | cut -f1 >> $logfile
        ;;

    "deltafs")
        # Start BBOS servers and clients

        message "BBOS Per-core log size: ${bb_log_size}MB"
        
        bb_server_list=$(cat $output_dir/bbos.hosts | tr '\n' ' ')
        n=1
        for s in $bb_server_list; do
            container_dir=$output_dir/bbos/containers.$n
            mkdir -p $container_dir

            # Copying config files for every server
            new_server_config=$output_dir/bbos/server.$n
            cp $bb_server_cfg $new_server_config
            echo $s >> $new_server_config
            echo $container_dir >> $new_server_config

            do_mpirun 1 "" "$s" "$bb_server $new_server_config" "$logfile" &
            
            message "BBOS server started at $s"

            sleep 5

            # Copying config files for every client of this server
            cp $bb_client_cfg $output_dir/bbos/client.$n
            echo $s >> $output_dir/bbos/client.$n

            n=$((n + 1))
        done

        c=1
        while [ $c -le $bb_clients ]; do
            s=$(((c % bb_servers) + 1))
            cfg_file=$output_dir/bbos/client.$s
            do_mpirun 1 "" "$bbos_nodes" \
                "$bb_client $c.obj $cfg_file $((bb_log_size * (2 ** 20))) $((bb_sst_size * (2 ** 20)))" \
                "$logfile" &

            message "BBOS client #$c started bound to server #$s"

            sleep 1

            c=$((c + 1))
        done

        # Start DeltaFS processes
        mkdir -p $output_dir/deltafs_$p/metadata || \
            die "deltafs metadata mkdir failed"
        mkdir -p $output_dir/deltafs_$p/data || \
            die "deltafs data mkdir failed"

        preload_lib_path="$umbrella_build_dir/deltafs-vpic-preload-prefix/src/"\
"deltafs-vpic-preload-build/src/libdeltafs-preload.so"
        deltafs_srvr_path="$umbrella_build_dir/deltafs-prefix/src/"\
"deltafs-build/src/server/deltafs-srvr"
        deltafs_srvr_ip=`hostname -i`

#        vars=("DELTAFS_MetadataSrvAddrs" "$deltafs_srvr_ip:10101"
#              "DELTAFS_FioName" "posix"
#              "DELTAFS_FioConf" "root=$output_dir/deltafs_$p/data"
#              "DELTAFS_Outputs" "$output_dir/deltafs_$p/metadata")
#
#        do_mpirun 1 vars[@] "$deltafs_nodes" "$deltafs_srvr_path" $logfile
#        if [ $? -ne 0 ]; then
#            die "deltafs server: mpirun failed"
#        fi
#
#        srvr_pid=$!

        vars=("LD_PRELOAD" "$preload_lib_path"
              "PRELOAD_Deltafs_root" "particle"
              "PRELOAD_Local_root" "${output_dir}"
              "PRELOAD_Bypass_deltafs_namespace" "1"
              "PRELOAD_Enable_verbose_error" "1"
              "SHUFFLE_Virtual_factor" "1024"
              "SHUFFLE_Mercury_proto" "bmi+tcp"
              "SHUFFLE_Subnet" "$ip_subnet")
#              "DELTAFS_MetadataSrvAddrs" "$deltafs_srvr_ip:10101"

        do_mpirun $cores vars[@] "$vpic_nodes" "$deck_dir/turbulence.op" $logfile
        if [ $? -ne 0 ]; then
#            kill -KILL $srvr_pid
            die "deltafs: mpirun failed"
        fi

#        kill -KILL $srvr_pid

        echo -n "Output size: " >> $logfile
        du -b $output_dir/deltafs_$p | tail -1 | cut -f1 >> $logfile

        # Kill BBOS clients and servers
        message ""
        message "Killing BBOS servers"
        for s in $bb_server_list; do
            message "- Killing BBOS server: $s"
            do_mpirun 1 "" "$s" "pkill -SIGINT bbos_server" "$logfile"
        done

        ;;
    esac
}

rm $logfile
gen_hosts

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
