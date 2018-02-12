/* This file contains all the necessary routines for tagging
 * and advancing tracer particles in VPIC. The methedology here
 * is to copy all tagged particles into a new tracer species.
 * This tracer species does not back react on the fields/plasma
 * so that the evolution of this tracer species is identical
 * to that of the actual simulation particles. The reason for
 * copying is that we can then use the q subfield in the
 * particle structure as a tag for the tracers, but by altering
 * q we cannot accuratly accumulate moments and thus there
 * can be no back-reaction and so we need to copy. The particle
 * push is unaffected since this depends on q/m which is
 * unmodified and specified external to the particle structure.
 *
 *
 *********************************************************************
 * Example of how these routines should be called form the input deck:
 *********************************************************************
 * begin_globals{
 *      ...
 *      species_t * tracers_list ;
 *      ...
 * }
 *
 * begin_initilization{
 *  ...
 *      define_species(,...);       // Because I'm lazy, tracer species
 *      ...                             // must be defined after all other
 *      ...                 // species !!
 *      tag_tracer(...);        // Call to tag/create a new tracer particle
 *      ...
 *      hijack_tracers(...);        // MUST BE CALLED AFTER DEFINING TRACERS
 *  ...
 * }
 *
 * begin_diagnostics{
 *      ...
 *      tag_tracer(...);                // Call to tag/create a new tracer particle
 *      ...
 *      dump_tracers(...);              // Call to dump tracers when required
 *      ...
 *      dump_tracer_restart();          // Call to dump tracers to a restart file...
 *      ...
 * }
 *
 * begin_particle_injection{
 *      advance_tracers();              // MUST BE CALLED AT EVERY TIME STEP
 *      ...
 *  }
 *
 */

// Here, p is the particle to copy and tag and tracer is the species
// itno which we will inject the tracer. If only one species is
// defined then tracer = tracers_list. tag shold be a unique
// 4 byte long identifier for this tracer; it is the only way to
// distinguish tracers! Also note that inject_particle_raw doesn't
// check for free storage space, so we're silently assuming all nodes
// have enough room to hold all tracers. Not a bad assumption since
// typically number of tracers will be small due to output size constaraints.
#define tag_tracer(p, tracer, tag) BEGIN_PRIMITIVE{           \
    float q = * reinterpret_cast<float *>( &tag ) ;           \
    inject_particle_raw(tracer, p->dx, p->dy, p->dz, p->i,    \
                        p->ux, p->uy, p->uz, q);              \
} END_PRIMITIVE

// I'm lazy, so here we are assuming tracers are at the head
// of the species list (defined last), and that we know
// exactly how many tracer species there are.
#define hijack_tracers( num_tracer_species )   BEGIN_PRIMITIVE{ \
    species_t *s         = species_list ;                       \
    global->tracers_list = species_list ;                       \
                                                                \
    for(int i = 1; i < num_tracer_species ; ++i) s = s->next;   \
                                                                \
    species_list = s->next;                                     \
    s->next      = NULL ;                                       \
                                                                \
} END_PRIMITIVE

// advance_p takes care of all local moves, but to resolve cross-domain
// moves and boundary interactions we need to call boundary_p as well.
#define advance_tracers() BEGIN_PRIMITIVE{          \
    static int a_initted = 0;                       \
    static accumulator_t *dummy_a;                  \
    if (a_initted == 0){                            \
        dummy_a = new_accumulators(grid);           \
        if(step) read_tracer_restart();             \
        a_initted = 1;                              \
    }                                               \
                                                    \
    species_t *s = global->tracers_list ;           \
    while( s ){                                     \
        s->nm += advance_p(s->p, s->np, s->q_m, s->pm, s->max_nm,   \
                   dummy_a, interpolator, grid) ;           \
        boundary_p(s, field, dummy_a, grid, rng);           \
        boundary_p(s, field, dummy_a, grid, rng);           \
        boundary_p(s, field, dummy_a, grid, rng);           \
        if (step%s->sort_interval == 0) sort_p(s,grid);     \
        s = s->next ;                                       \
    }                                                       \
} END_PRIMITIVE

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template <typename T> int is_negative(T val) {
    return (val < T(0));
}

#define nxg (grid->nx + 2)
#define nyg (grid->ny + 2)
#define nzg (grid->nz + 2)
#define i0 (ii%nxg)
#define j0 ((ii/nxg)%nyg)
#define k0 (ii/(nxg*nyg))
#define tracer_x ((i0 + (dx0-1)*0.5) * grid->dx + grid->x0)
#define tracer_y ((j0 + (dy0-1)*0.5) * grid->dy + grid->y0)
#define tracer_z ((k0 + (dz0-1)*0.5) * grid->dz + grid->z0)

// This dump routine probably isn't very fast, but it slightly reduces the
// output size of files and takes care of all the messy calculations we
// would have to do afterwards (local->global coordinate conversion and
// momentum->velocity)
#define dump_tracers(fbase) BEGIN_PRIMITIVE{          \
    char dname[256], fname[256] ;                     \
    float ex, ey, ez, bx, by, bz;                     \
    float dx0, dy0, dz0;                              \
    float ux, uy, uz;                                 \
    species_t *s = global->tracers_list ;             \
    const particle_t     * ALIGNED(32) p;             \
    const particle_t     * ALIGNED(32) p0;            \
    const interpolator_t * ALIGNED(16) f;             \
    const grid_t * g;                                 \
    int ii, j, n, nvar;                               \
    float *pout;                                      \
                                                      \
    FileIO fh;                                        \
    sprintf(dname, "%s/T.%d", fbase, step );          \
    dump_mkdir(dname);                                \
                                                      \
    nvar = 13;                                        \
    g = grid;                                         \
                                                      \
    while( s ){                                       \
        n = s->np;                                    \
        if ( n > 0 ){                                 \
            p0 = s->p;                                \
            pout = new float[n*nvar];                 \
            j = 0;                                    \
            for ( p=p0; n; n--, j++, p++ ){           \
                dx0 = p->dx;                          \
                dy0 = p->dy;                          \
                dz0 = p->dz;                          \
                ii = p->i;                            \
                ux = p->ux;                           \
                uy = p->uy;                           \
                uz = p->uz;                           \
                f = interpolator + ii;                \
                ex = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz); \
                ey = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx); \
                ez = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy); \
                bx = f->cbx + dx0*f->dcbxdx;               \
                by = f->cby + dy0*f->dcbydy;               \
                bz = f->cbz + dz0*f->dcbzdz;               \
                pout[j*nvar + 0] = (float) p->q;           \
                pout[j*nvar + 1] = (float) tracer_x ;      \
                pout[j*nvar + 2] = (float) tracer_y ;      \
                pout[j*nvar + 3] = (float) tracer_z ;      \
                pout[j*nvar + 4] = ux;                     \
                pout[j*nvar + 5] = uy;                     \
                pout[j*nvar + 6] = uz;                     \
                pout[j*nvar + 7] = ex;                     \
                pout[j*nvar + 8] = ey;                     \
                pout[j*nvar + 9] = ez;                     \
                pout[j*nvar + 10] = bx;                    \
                pout[j*nvar + 11] = by;                    \
                pout[j*nvar + 12] = bz;                    \
            }                                              \
            sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank()); \
            fh.open(fname, io_write);                      \
            WRITE(int,    s->np ,         fh ) ;           \
            WRITE(float,  s->q_m ,        fh ) ;           \
            WRITE(double, step*grid->dt , fh ) ;           \
            WRITE(float,  1.0 ,           fh ) ;           \
            WRITE(double, 1.0 ,           fh ) ;           \
            fh.write(pout, nvar*j);                        \
            fh.close();                                    \
            delete [] pout;                                \
        }                                                  \
        s = s->next;                                       \
    }                                                      \
} END_PRIMITIVE


// This dump routine save tracer information every time but only dump the
// tracer information every a few steps. This reduces I/O operations.
#define dump_tracers_multi(fbase) BEGIN_PRIMITIVE{          \
    char dname[256], fname[256] ;                           \
    float ex, ey, ez, bx, by, bz;                           \
    float dx0, dy0, dz0;                                    \
    static int a_initted = 0;                               \
    species_t *s = global->tracers_list ;                   \
    const particle_t     * ALIGNED(32) p;                   \
    const particle_t     * ALIGNED(32) p0;                  \
    const interpolator_t * ALIGNED(16) f;                   \
    const grid_t * g;                                       \
    int ii, j, n, nvar;                                     \
    static int numspecies = 0;                              \
    static float  *pout;                                    \
    static int    *np;                                      \
    static double *t;                                       \
    int Ntracer = global->Ntracer;                          \
    int tri = global->tracer_interval;                      \
    int ntr, step_tmp;                                      \
    static int *ntraj_point;                                \
    int isp = 0;                                            \
                                                            \
    ntr = Ntracer / nproc();                                \
    nvar = 13;                                              \
    g = grid;                                               \
                                                            \
    if (a_initted == 0) {                                   \
        while( s ){                                         \
            ++numspecies;                                   \
            s = s->next;                                    \
        }                                                   \
        pout = new float[2*numspecies*ntr*nvar*tri];        \
        np = new int[tri*numspecies];                       \
        t  = new double[tri];                               \
        ntraj_point  = new int[numspecies];                 \
        if (step != 0) read_tracer_multi_restart();         \
        s = global->tracers_list ;                          \
        a_initted = 1;                                      \
    }                                                       \
    step_tmp = step % tri;                                  \
    if (step_tmp == 0) {                                    \
        if (step != 0) tracer_output(fbase);                \
        memset(pout, 0, sizeof(float) * 2*numspecies*ntr*nvar*tri);  \
        memset(np,   0, sizeof(int) * tri*numspecies);          \
        memset(t,    0, sizeof(double) * tri);                  \
        memset(ntraj_point, 0, sizeof(int) * numspecies);       \
        s = global->tracers_list ;                              \
    }                                                           \
    t[step_tmp] = step * g->dt;                                 \
    int index = 0;                                              \
    isp = 0;                                                    \
    while( s ){                                                 \
        n = s->np;                                              \
        np[isp*tri + step_tmp] = n;                             \
        if ( n > 0 ){                                           \
            p0 = s->p;                                          \
            j = 0;                                              \
            for ( p=p0; n; n--, j++, p++ ){                     \
                dx0 = p->dx;                                    \
                dy0 = p->dy;                                    \
                dz0 = p->dz;                                    \
                ii = p->i;                                      \
                f  = interpolator + ii;                         \
                ex = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz); \
                ey = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx); \
                ez = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy); \
                bx = f->cbx + dx0*f->dcbxdx;                    \
                by = f->cby + dy0*f->dcbydy;                    \
                bz = f->cbz + dz0*f->dcbzdz;                    \
                index = ntraj_point[isp]*nvar + 2*isp*ntr*nvar*tri; \
                pout[index + 0] = (float) p->q;                 \
                pout[index + 1] = (float) tracer_x ;            \
                pout[index + 2] = (float) tracer_y ;            \
                pout[index + 3] = (float) tracer_z ;            \
                pout[index + 4] = p->ux;                        \
                pout[index + 5] = p->uy;                        \
                pout[index + 6] = p->uz;                        \
                pout[index + 7] = ex;                           \
                pout[index + 8] = ey;                           \
                pout[index + 9] = ez;                           \
                pout[index + 10] = bx;                          \
                pout[index + 11] = by;                          \
                pout[index + 12] = bz;                          \
                ntraj_point[isp]++;                             \
            }                                                   \
        }                                                       \
        s = s->next;                                            \
        isp++;                                                  \
    }                                                           \
    if (global->write_restart) {                                \
        dump_tracer_multi_restart();                            \
    }                                                           \
    if (global->write_end_restart || (step==num_step)) {        \
        if (global->write_end_restart) {                        \
            dump_tracer_multi_restart();                        \
        }                                                       \
        if (step == num_step) tracer_output(fbase);             \
        delete [] ntraj_point;                                  \
        delete [] pout;                                         \
        delete []   np;                                         \
        delete []    t;                                         \
    }                                                           \
} END_PRIMITIVE


#define tracer_output(fbase) BEGIN_PRIMITIVE{                       \
    FileIO fh;                                                      \
    sprintf(dname, "%s/T.%d", fbase, step );                        \
    dump_mkdir(dname);                                              \
    isp = 0;                                                        \
    while ( s ) {                                                   \
        sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank());   \
        fh.open(fname, io_write);                                   \
        WRITE(float, s->q_m, fh);                                   \
        WRITE(int,     nvar, fh);                                   \
        WRITE(int,      tri, fh);                                   \
        fh.write(np + isp*tri, tri);                                \
        fh.write(t, tri);                                           \
        fh.write(pout + 2*isp*ntr*nvar*tri, ntraj_point[isp]*nvar); \
        fh.close();                                                 \
        s = s->next;                                                \
        isp++;                                                      \
    }                                                               \
} END_PRIMITIVE


#define dump_tracer_multi_restart() BEGIN_PRIMITIVE{                \
    FileIO fh;                                                      \
    char fname[256];                                                \
    sprintf(fname, "restart%d/tracer_multi%d.%d",                   \
            global->restart_set, global->particle_tracing,          \
            (int)rank() ) ;                                         \
    fh.open(fname, io_write);                                       \
    fh.write(np, tri * numspecies);                                 \
    fh.write(t, tri);                                               \
    fh.write(ntraj_point, numspecies);                              \
    fh.write(pout, 2*numspecies*ntr*nvar*tri);                      \
    fh.close();                                                     \
} END_PRIMITIVE


#define read_tracer_multi_restart() BEGIN_PRIMITIVE{                \
    FileIO fh;                                                      \
    char fname[256];                                                \
    sprintf(fname, "restart%d/tracer_multi%d.%d",                   \
            global->restart_set, global->particle_tracing,          \
            (int)rank() ) ;                                         \
    fh.open(fname, io_read);                                        \
    fh.read(np, tri * numspecies);                                  \
    fh.read(t, tri);                                                \
    fh.read(ntraj_point, numspecies);                               \
    fh.read(pout, 2*numspecies*ntr*nvar*tri);                       \
    fh.close();                                                     \
} END_PRIMITIVE


#define dump_tracer_restart() BEGIN_PRIMITIVE{      \
    species_t *s = global->tracers_list ;           \
    char fname[256] ;                               \
    int numspecies = 0;                             \
    FileIO f;                                       \
    sprintf(fname, "restart%d/tracer%d.%d", global->restart_set,    \
            global->particle_tracing, (int)rank() ) ;               \
    f.open(fname,io_write);                     \
                                                \
    while( s ){                                 \
        ++numspecies;                           \
        s = s->next;                            \
    }                                           \
    WRITE(int, numspecies, f);                  \
    s = global->tracers_list;                   \
    while( s ){                                 \
        WRITE_STRING(s->name, f);               \
        WRITE(species_id, s->id, f);            \
        WRITE(float, s->q_m, f);                \
        WRITE(int, s->max_np, f);               \
        WRITE(int, s->max_nm, f);               \
        WRITE(int, s->sort_interval, f);        \
        WRITE(int, s->sort_out_of_place, f);    \
        WRITE(int, s->np, f);                   \
        if (s->np > 0 ) f.write(s->p, s->np);   \
        s = s->next;                            \
    }                                           \
    f.close();                                  \
} END_PRIMITIVE


#define read_tracer_restart()  BEGIN_PRIMITIVE{ \
    char fname[256] ;                           \
    FileIO f;                                   \
    int namelen, maxnp, maxnm, sort ;           \
    int sorttype, np, numspecies, i ;           \
    float qm ;                                  \
    species_id id;                              \
    char *name;                                 \
    species_t * s;                              \
                                                \
    global->tracers_list = NULL ;               \
    sprintf(fname, "restart%d/tracer%d.%d", global->restart_set,    \
            global->particle_tracing, (int)rank() ) ;               \
    f.open(fname,io_read);                      \
                                                \
    READ(int, numspecies, f);                   \
    for( i=0 ; i < numspecies ; ++i){           \
        READ(int, namelen, f) ;                 \
        name = (char *)malloc(namelen+1) ;      \
        f.read(name, namelen) ;                 \
        name[namelen] = '\0' ;                  \
        READ(species_id,    id,   f);           \
        READ(float,     qm,   f);               \
        READ(int,       maxnp,    f);           \
        READ(int,       maxnm,    f);           \
        READ(int,       sort,     f);           \
        READ(int,       sorttype, f);           \
        READ(int,       np,   f);               \
        s = new_species(name, qm, maxnp, maxnm, sort, sorttype, \
                &(global->tracers_list));                       \
        s->id = id ;                                            \
        s->np = np ;                                            \
        if ( np > 0 ) f.read(s->p, np);                         \
    }                                                           \
    f.close();                                                  \
} END_PRIMITIVE

// dump tracer by particle trajectory
#define dump_traj(fbase) BEGIN_PRIMITIVE{   \
    char dname[256], fname[256] ;           \
    species_t *s = global->tracers_list ;   \
    float ex, ey, ez, bx, by, bz;           \
    float dx0, dy0, dz0;                    \
    float ux, uy, uz, q;                    \
    int ii, n, nvar;                        \
    const particle_t     * ALIGNED(32) p;   \
    const particle_t     * ALIGNED(32) p0;  \
    const interpolator_t * ALIGNED(16) f;   \
    const grid_t * g;                       \
    FileIO fh;                              \
                                            \
    sprintf(dname, "%s", fbase );           \
    dump_mkdir(dname);                      \
    nvar = 13;                              \
    float pout[nvar];                       \
    g = grid;                               \
                                            \
    while( s ){                             \
        n = s->np;                          \
        if ( n > 0 ){                       \
            p0 = s->p;                      \
            for ( p=p0; n; n--, p++ ){                      \
                dx0 = p->dx;                                \
                dy0 = p->dy;                                \
                dz0 = p->dz;                                \
                ii = p->i;                                  \
                ux = p->ux;                                 \
                uy = p->uy;                                 \
                uz = p->uz;                                 \
                q = p->q;                                   \
                int tag = *reinterpret_cast<int*>(&q);      \
                f = interpolator + ii;                      \
                if (tag != 0) {                             \
                    ex = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz); \
                    ey = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx); \
                    ez = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy); \
                    bx = f->cbx + dx0*f->dcbxdx;                \
                    by = f->cby + dy0*f->dcbydy;                \
                    bz = f->cbz + dz0*f->dcbzdz;                \
                    sprintf(fname, "%s/%s.%i", dname , s->name, tag);   \
                    fh.open(fname,io_append);                   \
                    pout[0] = step*grid->dt ;                   \
                    pout[1] = (float) tracer_x ;                \
                    pout[2] = (float) tracer_y ;                \
                    pout[3] = (float) tracer_z ;                \
                    pout[4] = ux;                               \
                    pout[5] = uy;                               \
                    pout[6] = uz;                               \
                    pout[7] = ex;                               \
                    pout[8] = ey;                               \
                    pout[9] = ez;                               \
                    pout[10] = bx;                              \
                    pout[11] = by;                              \
                    pout[12] = bz;                              \
                    fh.write(pout,nvar);                        \
                    fh.close();                                 \
                }                                               \
            }                                                   \
        }                                                       \
        s = s->next;                                            \
    }                                                           \
} END_PRIMITIVE


/****************************************************************
 ******************* Copied from dumpmacros.h *******************
 ***************************************************************/

#define WRITE(type,value,fileIO) BEGIN_PRIMITIVE {   \
    type __WRITE_tmp = (type)(value);                \
    fileIO.write( &__WRITE_tmp, 1 );                 \
} END_PRIMITIVE

// Note: strlen does not include the terminating NULL
#define WRITE_STRING(string,fileIO) BEGIN_PRIMITIVE {       \
    int __WRITE_STRING_len = 0;                             \
    if( string!=NULL ) __WRITE_STRING_len = strlen(string); \
    fileIO.write( &__WRITE_STRING_len, 1 );                 \
    if( __WRITE_STRING_len>0 ) {                            \
        fileIO.write( string, __WRITE_STRING_len );         \
    }                                                       \
} END_PRIMITIVE

#define READ(type,value,fileIO) BEGIN_PRIMITIVE {   \
    type __READ_tmp;                                \
    fileIO.read(&__READ_tmp, 1 );                   \
    (value) = __READ_tmp;                           \
} END_PRIMITIVE

