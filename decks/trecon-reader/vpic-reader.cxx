#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <dirent.h>
#include <fcntl.h>
#include <list>
#include <map>
#include <assert.h>

#include <deltafs/deltafs_api.h>

using namespace std;

typedef map<int,FILE*> FileMap;
typedef map<int,int64_t> ParticleMap;
typedef list<int> EpochList;

char *me;

void usage(int ret)
{
    printf("\n"
           "usage: %s [options] -i input_dir -o output_dir\n"
           "  options:\n"
           "    -d        Run in DeltaFS mode\n"
           "    -n num    Number of particles to read (reading ID 1 to num)\n"
           "    -h        This usage info\n"
           "\n",
           me);

    exit(ret);
}

void close_files(FileMap *out)
{
    FileMap::iterator it;

    for (it = out->begin(); it != out->end(); ++it) {
        if (it->second)
            fclose(it->second);
    }
}

int generate_files(char *outdir, long long int num, FileMap *out)
{
    char fpath[PATH_MAX];

    for (long long int i = 1; i <= num; i++) {
        if (!snprintf(fpath, PATH_MAX, "%s/particle%lld.txt", outdir, i)) {
            perror("Error: snprintf failed");
            return 1;
        }

        if (!((*out)[i] = fopen(fpath, "w"))) {
            perror("Error: fopen failed");
            close_files(out);
            return 1;
        }
    }

    return 0;
}

int deltafs_read_particles(long long int num, char *indir, char *outdir)
{
    int ret;
    deltafs_plfsdir_t *dir;
    char *file_data, fname[PATH_MAX];
    FileMap out;
    size_t len;

    /* Iterate through epoch frames */
    if (generate_files(outdir, num, &out)) {
        fprintf(stderr, "Error: particle trajectory file creation failed\n");
        return 1;
    }

    dir = deltafs_plfsdir_create_handle(O_RDONLY);

    if ((ret = deltafs_plfsdir_open(dir, indir, NULL))) {
        perror("Error: cannot open input PLFS directory");
        deltafs_plfsdir_free_handle(dir);
        goto err;
    }

    for (int i=1; i<=num; i++) {
        /* Determine fname for particle */
        if (!snprintf(fname, PATH_MAX, "eparticle.%016lx", (long int) i)) {
            perror("Error: snprintf failed");
            goto err;
        }

        if (!(file_data = (char*) deltafs_plfsdir_readall(dir, fname, &len))) {
            perror("Error: failed to read particle data");
            deltafs_plfsdir_free_handle(dir);
            goto err;
        }

        /* Write out particle trajectory data */
        if (fwrite(file_data, 1, len, out[i]) != len) {
            perror("Error: fwrite failed");
            goto err;
        }

        free(file_data);
    }

    deltafs_plfsdir_free_handle(dir);
    close_files(&out);
    return 0;

err:
    close_files(&out);
    return 1;
}

int pick_particles(char *ppath, int epoch, long long int num, ParticleMap *ids)
{
    DIR *d;
    struct dirent *dp;
    char epath[PATH_MAX];
    char fprefix[PATH_MAX];
    char fpath[PATH_MAX];
    long long int cur = 1;

    if (!snprintf(epath, PATH_MAX, "%s/T.%d", ppath, epoch)) {
        perror("Error: snprintf for epath failed");
        return 1;
    }

    if (!snprintf(fprefix, PATH_MAX, "eparticle.%d.", epoch)) {
        perror("Error: snprintf for fprefix failed");
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

        if (dp->d_type != DT_REG)
            continue;

        if (strncmp(dp->d_name, fprefix, strnlen(fprefix, PATH_MAX))) {
            fprintf(stderr, "Warning: unexpected file %s in %s\n",
                    dp->d_name, epath);
            continue;
        }

        if (!snprintf(fpath, PATH_MAX, "%s/%s", epath, dp->d_name)) {
            perror("Error: snprintf for fpath failed");
            goto err;
        }

        if (!(fp = fopen(fpath, "rb"))) {
            perror("Error: fopen epoch file failed");
            goto err;
        }

        /* Verify V0 header */
        assert(!fseek(fp, 5 * sizeof(char), SEEK_CUR));
        assert(fread(&x, sizeof(short int), 1, fp) && x == 0xcafe);
        assert(fread(&x, sizeof(int), 1, fp) && x == 0xdeadbeef);
        assert(!fseek(fp, sizeof(float) + sizeof(double), SEEK_CUR));
        assert(fread(&x, sizeof(int), 1, fp) && x == 0);
        assert(!fseek(fp, (11*sizeof(float)) + (8*sizeof(int)), SEEK_CUR));

        /* Read array header */
        assert(fread(&wsize, sizeof(int), 1, fp));
        assert(fread(&x, sizeof(int), 1, fp) && x == 1);
        assert(fread(&wnum, sizeof(int), 1, fp));

        /*
         * Particle structure (species_advance.h):
         * - float dx, dy, dz; // Particle position, cell coordinates ([-1,1])
         * - int32_t i;
         * - float ux, uy, uz; // Particle normalized momentum
         * - float q;          // Particle charge
         * - int64_t tag, tag2; // particle identification tags
         */
#define DATA_LEN (7*sizeof(float) + sizeof(int32_t) + 2*sizeof(int64_t))
#define TAG_OFFT (7*sizeof(float) + sizeof(int32_t))
        char data[DATA_LEN];
        int64_t tag;

        for (int i = 1; i <= wnum; i++) {
            assert(fread(data, 1, DATA_LEN, fp) == DATA_LEN);
            memcpy(&tag, data + TAG_OFFT, sizeof(int64_t));

            if ((tag & 0x3ffffffffff) == cur) {
                (*ids)[cur] = tag;
                printf("Particle #%lld: ID 0x%016llx\n",
                       cur, (long long unsigned int) tag);
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

int process_epoch(char *ppath, int it, ParticleMap ids, FileMap out)
{
    DIR *d;
    struct dirent *dp;
    char epath[PATH_MAX];
    char fprefix[PATH_MAX];
    char fpath[PATH_MAX];
    long long int num = out.size();

    if (!snprintf(epath, PATH_MAX, "%s/T.%d", ppath, it)) {
        perror("Error: snprintf for epath failed");
        return 1;
    }

    if (!snprintf(fprefix, PATH_MAX, "eparticle.%d.", it)) {
        perror("Error: snprintf for fprefix failed");
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

        if (dp->d_type != DT_REG)
            continue;

        if (strncmp(dp->d_name, fprefix, strnlen(fprefix, PATH_MAX))) {
            fprintf(stderr, "Warning: unexpected file %s in %s\n",
                    dp->d_name, epath);
            continue;
        }

        //printf("Found file %s, epoch %d.\n", dp->d_name, it);

        if (!snprintf(fpath, PATH_MAX, "%s/%s", epath, dp->d_name)) {
            perror("Error: snprintf for fpath failed");
            goto err;
        }

        if (!(fp = fopen(fpath, "rb"))) {
            perror("Error: fopen epoch file failed");
            goto err;
        }

        /* Verify V0 header */
        assert(!fseek(fp, 5 * sizeof(char), SEEK_CUR));
        assert(fread(&x, sizeof(short int), 1, fp) && x == 0xcafe);
        assert(fread(&x, sizeof(int), 1, fp) && x == 0xdeadbeef);
        assert(!fseek(fp, sizeof(float) + sizeof(double), SEEK_CUR));
        assert(fread(&x, sizeof(int), 1, fp) && x == 0);
        assert(!fseek(fp, (11*sizeof(float)) + (8*sizeof(int)), SEEK_CUR));

        /* Read array header */
        assert(fread(&wsize, sizeof(int), 1, fp));
        assert(fread(&x, sizeof(int), 1, fp) && x == 1);
        assert(fread(&wnum, sizeof(int), 1, fp));
        //printf("Array: %d elements, %db each\n", wnum, wsize);

        /*
         * Particle structure (species_advance.h):
         * - float dx, dy, dz; // Particle position, cell coordinates ([-1,1])
         * - int32_t i;
         * - float ux, uy, uz; // Particle normalized momentum
         * - float q;          // Particle charge
         * - int64_t tag, tag2; // particle identification tags
         */
#define DATA_LEN (7*sizeof(float) + sizeof(int32_t) + 2*sizeof(int64_t))
#define TAG_OFFT (7*sizeof(float) + sizeof(int32_t))
        char data[DATA_LEN];
        int64_t tag;

        for (int i = 1; i <= wnum; i++) {
            assert(fread(data, 1, DATA_LEN, fp) == DATA_LEN);
            memcpy(&tag, data + TAG_OFFT, sizeof(int64_t));

            int idx = (int)(tag & 0x3ffffffffff);
            if (ids[idx] && ids[idx] == tag) {
                char preamble[64];

                /* Write out particle data */
                assert(sprintf(preamble, "Epoch: %d\nTag: 0x%016llx\nData:", it,
                               (long long int) tag));
                assert(fwrite(preamble, 1, strlen(preamble), out[idx]));
                assert(fwrite(data, 1, DATA_LEN, out[idx]));
                assert(fwrite("\n\n", 1, 2, out[idx]));
                //printf("Found 0x%016llx in epoch %d.\n", (long long unsigned int) tag, it);
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

int read_particles(long long int num, char *indir, char *outdir)
{
    DIR *in;
    struct dirent *dp;
    char ppath[PATH_MAX];
    EpochList epochs;
    EpochList::iterator it;
    FileMap out;
    ParticleMap ids;

    printf("Reading particles from %s.\n", indir);
    printf("Storing trajectories in %s.\n", outdir);

    /* Open particle directory and sort epoch directories */
    if (!snprintf(ppath, PATH_MAX, "%s/particle", indir)) {
        perror("Error: snprintf failed");
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

    /* Iterate through epoch frames */
    if (generate_files(outdir, num, &out)) {
        fprintf(stderr, "Error: particle trajectory file creation failed\n");
        return 1;
    }

    /* Pick the particle IDs to query */
    if (pick_particles(ppath, *(epochs.begin()), num, &ids)) {
        fprintf(stderr, "Error: failed to pick particles to query\n");
        close_files(&out);
        return 1;
    }

    for (it = epochs.begin(); it != epochs.end(); ++it) {
        printf("Processing epoch %d.\n", *it);
        if (process_epoch(ppath, *it, ids, out)) {
            fprintf(stderr, "Error: epoch data processing failed\n");
            close_files(&out);
            return 1;
        }
    }

    close_files(&out);
    return 0;
}

int main(int argc, char **argv)
{
    int ret, c, d = 0;
    long long int num = 1;
    char indir[PATH_MAX], outdir[PATH_MAX];

    me = argv[0];
    indir[0] = outdir[0] = '\0';

    while ((c = getopt(argc, argv, "dhi:n:o:p:")) != -1) {
        switch(c) {
        case 'd': /* run in DeltaFS mode */
            d = 1;
            break;
        case 'h': /* print help */
            usage(0);
        case 'i': /* input directory (VPIC output) */
            if (!strncpy(indir, optarg, PATH_MAX)) {
                perror("Error: invalid input dir");
                usage(1);
            }
            break;
        case 'n': /* number of particles to fetch */
            char *end;
            num = strtoll(optarg, &end, 10);
            if (*end) {
                perror("Error: invalid num argument");
                usage(1);
            }
            break;
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

    if (!indir[0] || !outdir[0]) {
        fprintf(stderr, "Error: input and output directories are mandatory\n");
        usage(1);
    }

    /* Do particle things */
    if (d)
        ret = deltafs_read_particles(num, indir, outdir);
    else
        ret = read_particles(num, indir, outdir);

    return ret;
}
