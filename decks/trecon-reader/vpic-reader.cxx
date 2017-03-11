#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <dirent.h>
#include <fcntl.h>
#include <sys/time.h>
#include <mpi.h>

#include <list>
#include <map>
#include <set>

using namespace std;

typedef map<int64_t,int64_t> ParticleMap;
typedef map<int64_t,int64_t> RevParticleMap;
typedef list<int> EpochList;

char *me;
int rank, worldsz;

void usage(int ret)
{
    printf("\n"
           "usage: %s [options] -i input_dir\n"
           "  options:\n"
           "    -o dir    Output directory, /dev/null if unspecified\n"
           "    -n num    Number of particles to read (reading ID 1 to num)\n"
           "              If unspecified we explore 10**i range, i in {0, k}\n"
           "              with 10**k approaching total num of particles\n"
           "    -r num    Number of query retries (before averaging, def. 3)\n"
           "    -h        This usage info\n"
           "\n",
           me);

    exit(ret);
}

/*
 * Particle structure (species_advance.h):
 * - float dx, dy, dz; // Particle position, cell coordinates ([-1,1])
 * - int32_t i;
 * - float ux, uy, uz; // Particle normalized momentum
 * - float q;          // Particle charge
 * - int64_t tag, tag2; // particle identification tags
 *
 * On Emulab, the VPIC output file preamble is 115B.
 * The particle output per frame is 48B.
 */
#define DATA_LEN (7*sizeof(float) + sizeof(int32_t) + 2*sizeof(int64_t))
#define TAG_OFFT (7*sizeof(float) + sizeof(int32_t))

int process_file_metadata(FILE *fp, int *wsize, int *wnum)
{
    int x = 0;

    /* Verify V0 header */
    if (fseek(fp, 5 * sizeof(char), SEEK_CUR)) return 1;
    if (!fread(&x, sizeof(short int), 1, fp) || x != 0xcafe) return 1;
    if (!fread(&x, sizeof(int), 1, fp) || x != 0xdeadbeef) return 1;
    if (fseek(fp, sizeof(float) + sizeof(double), SEEK_CUR)) return 1;
    if (!fread(&x, sizeof(int), 1, fp) || x != 0) return 1;
    if (fseek(fp, (11*sizeof(float)) + (8*sizeof(int)), SEEK_CUR)) return 1;

    /* Read array header */
    if (!fread(wsize, sizeof(int), 1, fp)) return 1;
    if (!fread(&x, sizeof(int), 1, fp) || x != 1) return 1;
    if (!fread(wnum, sizeof(int), 1, fp)) return 1;

    return 0;
}

int pick_particles(char *ppath, int epoch, int64_t num, ParticleMap *ids,
                   RevParticleMap *rids)
{
    DIR *d;
    struct dirent *dp;
    char epath[PATH_MAX];
    char fprefix[PATH_MAX];
    char fpath[PATH_MAX];
    int64_t cur = 1;

    if (snprintf(epath, PATH_MAX, "%s/T.%d", ppath, epoch) <= 0) {
        fprintf(stderr, "Error: snprintf for epath failed\n");
        return 1;
    }

    if (snprintf(fprefix, PATH_MAX, "particle.%d.", epoch) <= 0 ) {
        fprintf(stderr, "Error: snprintf for fprefix failed\n");
        return 1;
    }

    if ((d = opendir(epath)) == NULL) {
        perror("Error: cannot open epoch directory");
        return 1;
    }

    /* Open each per-process file and process it */
    while (dp = readdir(d)) {
        FILE *fp;
        int x = 0, wsize, wnum;
        char data[DATA_LEN];
        int64_t tag;

        if (dp->d_type != DT_REG)
            continue;

        if (strncmp(dp->d_name+2, fprefix, strnlen(fprefix, PATH_MAX))) {
            fprintf(stderr, "Warning: unexpected file %s in %s\n",
                    dp->d_name, epath);
            continue;
        }

        if (snprintf(fpath, PATH_MAX, "%s/%s", epath, dp->d_name) <= 0) {
            fprintf(stderr, "Error: snprintf for fpath failed\n");
            goto err;
        }

        if (!(fp = fopen(fpath, "rb"))) {
            perror("Error: fopen epoch file failed");
            goto err;
        }

        if (process_file_metadata(fp, &wsize, &wnum)) {
            fclose(fp);
            goto err;
        }

        //printf("Array: %d elements, %db each\n", wnum, wsize);

        for (int i = 1; i <= wnum; i++) {
            if (fread(data, 1, DATA_LEN, fp) != DATA_LEN) {
                fclose(fp);
                goto err;
            }

            memcpy(&tag, data + TAG_OFFT, sizeof(int64_t));

            if ((*rids).find(tag) == (*rids).end()) {
                (*ids)[cur] = tag;
                (*rids)[tag] = cur;
                //printf("Particle #%ld: ID 0x%016lx\n", cur, tag);
                cur++;

                if (cur > num) {
                    fclose(fp);
                    closedir(d);
                    return 0;
                }
            }
        }

        fclose(fp);
    }

    closedir(d);
    return 0;

err:
    closedir(d);
    return 1;
}

int process_epoch(char *ppath, char *outdir, int it,
                  ParticleMap ids, RevParticleMap rids)
{
    DIR *d;
    FILE *fp, *fd;
    struct dirent *dp;
    int wsize, wnum;
    int64_t idx, tag, num = ids.size();
    char epath[PATH_MAX], fprefix[PATH_MAX], fpath[PATH_MAX];
    char wpath[PATH_MAX], data[DATA_LEN], preamble[64];

    if (snprintf(epath, PATH_MAX, "%s/T.%d", ppath, it) <= 0) {
        fprintf(stderr, "Error: snprintf for epath failed\n");
        return 1;
    }

    if (snprintf(fprefix, PATH_MAX, "particle.%d.", it) <= 0) {
        fprintf(stderr, "Error: snprintf for fprefix failed\n");
        return 1;
    }

    if ((d = opendir(epath)) == NULL) {
        perror("Error: cannot open epoch directory");
        return 1;
    }

    /* Open each per-process file and process it */
    while (dp = readdir(d)) {
        if (dp->d_type != DT_REG)
            continue;

        if (strncmp(dp->d_name+2, fprefix, strnlen(fprefix, PATH_MAX))) {
            fprintf(stderr, "Warning: unexpected file %s in %s\n",
                    dp->d_name, epath);
            continue;
        }

        //printf("Found file %s, epoch %d.\n", dp->d_name, it);

        if (snprintf(fpath, PATH_MAX, "%s/%s", epath, dp->d_name) <= 0) {
            fprintf(stderr, "Error: snprintf for fpath failed\n");
            goto err;
        }

        if (!(fp = fopen(fpath, "rb"))) {
            perror("Error: fopen epoch file failed");
            goto err;
        }

        if (process_file_metadata(fp, &wsize, &wnum))
            goto err_fd;

        for (int i = 1; i <= wnum; i++) {
            if (!fread(data, 1, DATA_LEN, fp) == DATA_LEN) {
                perror("Error: fread failed");
                goto err_fd;
            }

            memcpy(&tag, data + TAG_OFFT, sizeof(int64_t));

            if ((rids.find(tag) != rids.end())) {
                idx = rids[tag];

                if (!outdir[0])
                    continue;

                if (!snprintf(wpath, PATH_MAX, "%s/particle%ld.txt", outdir, idx)) {
                    perror("Error: snprintf for wpath failed");
                    return 1;
                }

                if (!(fd = fopen(wpath, "w"))) {
                    perror("Error: fopen failed for output");
                    return 1;
                }

                /* Write out particle data */
                if (sprintf(preamble, "Epoch: %d\nTag: 0x%016lx\nData:", it, tag) <= 0) {
                    fprintf(stderr, "Error: sprintf for preamble failed\n");
                    fclose(fd);
                    goto err_fd;
                }

                if (!fwrite(preamble, 1, strlen(preamble), fd) ||
                    !fwrite(data, 1, DATA_LEN, fd) ||
                    !fwrite("\n\n", 1, 2, fd)) {
                    perror("Error: fwrite failed");
                    fclose(fd);
                    goto err_fd;
                }

                //printf("Found 0x%016lx in epoch %d.\n", tag, it);

                fclose(fd);
            }
        }

        fclose(fp);
    }

    closedir(d);
    return 0;

err_fd:
    fclose(fp);
err:
    closedir(d);
    return 1;
}

int read_particles(int64_t num, char *indir, char *outdir)
{
    DIR *in;
    struct dirent *dp;
    char ppath[PATH_MAX];
    EpochList epochs;
    EpochList::iterator it;
    ParticleMap ids;
    RevParticleMap rids;

    //printf("Reading particles from %s.\n", indir);
    //printf("Storing trajectories in %s.\n", outdir);

    /* Open particle directory and sort epoch directories */
    if (snprintf(ppath, PATH_MAX, "%s/particle", indir) <= 0) {
        fprintf(stderr, "Error: snprintf for ppath failed\n");
        return 1;
    }

    if ((in = opendir(ppath)) == NULL) {
        perror("Error: cannot open input directory");
        return 1;
    }

    while (dp = readdir(in)) {
        int epoch;
        char *end;

        if (dp->d_type != DT_DIR)
            continue;

        if (!strcmp(dp->d_name, ".") ||
            !strcmp(dp->d_name, "..") ||
            dp->d_name[0] != 'T')
            continue;

        epoch = strtoll(dp->d_name+2, &end, 10);
        if (*end) {
            perror("Error: strtoll failed");
            closedir(in);
            return 1;
        }

        //printf("Found subdir %s, epoch %d.\n", dp->d_name, epoch);

        epochs.push_back(epoch);
    }

    closedir(in);
    epochs.sort();

    /* Pick the particle IDs to query */
    if (pick_particles(ppath, *(epochs.begin()), num, &ids, &rids))
        return 1;

    for (it = epochs.begin(); it != epochs.end(); ++it) {
        if ((*it) / (*epochs.begin()) == (rank + 1)) {
            //printf("Rank %d processing epoch %d.\n", rank, *it);
            if (process_epoch(ppath, outdir, *it, ids, rids)) {
                fprintf(stderr, "Error: epoch data processing failed\n");
                return 1;
            }
        }
    }

    return 0;
}

int64_t get_total_particles(char *indir)
{
    int64_t total = 0;
    char infop[PATH_MAX], buf[256];
    FILE *fd;

    if (snprintf(infop, PATH_MAX, "%s/info", indir) <= 0) {
        fprintf(stderr, "Error: snprintf for infop failed\n");
        return 1;
    }

    if (!(fd = fopen(infop, "r"))) {
        perror("Error: fopen failed for info");
        return 1;
    }

    while (fgets(buf, sizeof(buf), fd) != NULL) {
        char *str;

        if ((str = strstr(buf, "total # of particles = ")) != NULL) {
            total = (int64_t) strtof(str + 23, NULL);
            break;
        }
    }

    fclose(fd);
    return total;
}

int query_particles(int64_t retries, int64_t num, char *indir, char *outdir)
{
    int ret = 0;
    struct timeval ts, te;
    int64_t elapsed_sum = 0, max_elapsed_avg = 0;

    if (rank == 0)
        printf("Querying %ld particles (%ld retries)\n", num, retries);

    for (int64_t i = 1; i <= retries; i++) {
        int64_t elapsed, max_elapsed;
        int64_t *elapsed_all;

        if (rank == 0 &&
            !(elapsed_all = (int64_t *)malloc(sizeof(int64_t) * worldsz))) {
            perror("Error: malloc failed");
            exit(1);
        }

        gettimeofday(&ts, 0);
        ret = read_particles(num, indir, outdir);
        gettimeofday(&te, 0);

        elapsed = (te.tv_sec-ts.tv_sec)*1000 + (te.tv_usec-ts.tv_usec)/1000;

        MPI_Gather(&elapsed, 1, MPI_LONG_LONG_INT, elapsed_all, 1,
                   MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //printf("(Rank %d, %ld) %ldms / query, %ld ms / particle\n",
        //       rank, i, elapsed, elapsed / num);
        if (rank == 0) {
            elapsed = max_elapsed = 0;
            for (int j = 0; j < worldsz; j++) {
                elapsed += elapsed_all[j];
                if (max_elapsed < elapsed_all[j])
                    max_elapsed = elapsed_all[j];
            }
            elapsed /= worldsz;
            printf("Overall: %ldms / query, %ld ms / particle\n",
                   max_elapsed, elapsed / num);
            free(elapsed_all);
        }

        elapsed_sum += elapsed;
        max_elapsed_avg += max_elapsed;
    }

    if (rank == 0)
        printf("Querying results: %ld ms / query, %ld ms / particle\n\n",
                elapsed_sum / retries, elapsed_sum / num / retries);

    return ret;
}

int main(int argc, char **argv)
{
    int ret, c;
    int64_t num = 0, retries = 3, total = 0;
    char indir[PATH_MAX], outdir[PATH_MAX];
    char *end;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "Error: MPI_Init failed\n");
        return 1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &worldsz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    me = argv[0];
    indir[0] = outdir[0] = '\0';

    while ((c = getopt(argc, argv, "hi:n:o:p:r:")) != -1) {
        switch(c) {
        case 'h': /* print help */
            usage(0);
        case 'i': /* input directory (VPIC output) */
            if (!strncpy(indir, optarg, PATH_MAX)) {
                perror("Error: invalid input dir");
                usage(1);
            }
            break;
        case 'n': /* number of particles to fetch */
            num = strtoll(optarg, &end, 10);
            if (*end) {
                perror("Error: invalid num argument");
                usage(1);
            }
            break;
        case 'r': /* number of query retries */
            retries = strtoll(optarg, &end, 10);
            if (*end) {
                perror("Error: invalid retry argument");
                usage(1);
            }
        case 'o': /* output directory (trajectory files) */
            if (!strncpy(outdir, optarg, PATH_MAX)) {
                perror("Error: invalid output dir");
                usage(1);
            }
            break;
        default:
            usage(1);
        }
    }

    if (!indir[0]) {
        fprintf(stderr, "Error: input directory unspecified\n");
        usage(1);
    }

    /* Get total number of particles */
    total = get_total_particles(indir);
    if (!total) {
        fprintf(stderr, "Error: failed to read the total number of particles\n");
        return 1;
    }

    if (rank == 0)
        printf("\nNumber of particles: %ld\n", total);
    /* XXX: The following is only until we figure out caching */
    if (total > 1e5) {
        total = 1e5;
        if (rank == 0)
            printf("Warning: will stop querying at 100K particles\n");
    }
    if (rank == 0)
        printf("\n");

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Go through the query dance: increment num from 1 to total particles
     * multiplying by 10 each time, unless num is specified.
     */
    if (num) {
        ret = query_particles(retries, num, indir, outdir);
    } else {
        num = 1;

        while (num <= total) {
            ret = query_particles(retries, num, indir, outdir);
            if (ret)
                return ret;

            num *= 10;
        }
    }

    MPI_Finalize();

    return ret;
}
