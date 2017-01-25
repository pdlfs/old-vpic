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
 * 	...
 *      define_species(,...);  		// Because I'm lazy, tracer species
 *      ...                            	// must be defined after all other
 *      ...			       	// species !!
 *      tag_tracer(...);		// Call to tag/create a new tracer particle
 *      ...
 *      hijack_tracers(...);   		// MUST BE CALLED AFTER DEFINING TRACERS
 *	...
 * }
 *
 * begin_diagnostics{
 *      ...
 *      tag_tracer(...);               	// Call to tag/create a new tracer particle
 *      ...
 *      dump_tracers(...);             	// Call to dump tracers when required
 *      ...
 *      dump_tracer_restart();      	// Call to dump tracers to a restart file...
 *      ...
 * }
 *
 * begin_particle_injection{
 *      advance_tracers();          	// MUST BE CALLED AT EVERY TIME STEP
 *      ...
 *  }
 *
 */

#include <sys/types.h>
#include <dirent.h> /* Needed for opendir */
#include <unistd.h> /* Needed for getcwd */
#include <assert.h>

// Here, p is the particle to copy and tag and tracer is the species 
// itno which we will inject the tracer. If only one species is 
// defined then tracer = tracers_list. tag shold be a unique 
// 4 byte long identifier for this tracer; it is the only way to 
// distinguish tracers! Also note that inject_particle_raw doesn't 
// check for free storage space, so we're silently assuming all nodes 
// have enough room to hold all tracers. Not a bad assumption since
// typically number of tracers will be small due to output size constaraints.
inline void tag_tracer(particle_t *p, species_t *tracer, long tag) {
  particle_t *tp = tracer->p + (tracer->np++);
  tp->dx = p->dx; tp->dy = p->dy; tp->dz = p->dz; tp->i = p->i;
  tp->ux = p->ux; tp->uy = p->uy; tp->uz = p->uz; tp->q = 0;
  tp->tag = tag;
}

// I'm lazy, so here we are assuming tracers are at the head
// of the species list (defined last), and that we know
// exactly how many tracer species there are.
#define hijack_tracers( num_tracer_species )   BEGIN_PRIMITIVE{ \
  species_t *s         = species_list ;				\
  global->tracers_list = species_list ;				\
								\
  for(int i = 1; i < num_tracer_species ; ++i) s = s->next;	\
  								\
  species_list = s->next;					\
  s->next      = NULL ;						\
								\
} END_PRIMITIVE

// advance_p takes care of all local moves, but to resolve cross-domain
// moves and boundary interactions we need to call boundary_p as well.
#define advance_tracers() BEGIN_PRIMITIVE{			\
  static int a_initted = 0;					\
  static accumulator_t *dummy_a;				\
  if (a_initted == 0){						\
    dummy_a = new_accumulators(grid);				\
    if(step) read_tracer_restart();				\
    a_initted = 1;						\
  }								\
								\
  species_t *s = global->tracers_list ;				\
  while( s ){							\
    s->nm += advance_p(s->p, s->np, s->q_m, s->pm, s->max_nm, 	\
		       dummy_a, interpolator, grid) ;		\
    boundary_p(s, field, dummy_a, grid, rng);			\
    boundary_p(s, field, dummy_a, grid, rng);			\
    boundary_p(s, field, dummy_a, grid, rng);			\
    if (step%s->sort_interval == 0) sort_p(s,grid);		\
    s = s->next ;						\
  }								\
								\
} END_PRIMITIVE



#define tracer_x ( p[j].i%(grid->nx+2)                + (p[j].dx-1)/2.0 ) * grid->dx + grid->x0
#define tracer_y ((p[j].i/(grid->nx+2))%(grid->ny+2)  + (p[j].dy-1)/2.0 ) * grid->dy + grid->y0
#define tracer_z ( p[j].i/((grid->nx+2)*(grid->ny+2)) + (p[j].dz-1)/2.0 ) * grid->dz + grid->z0

// This dump routine probably isn't very fast, but it slightly reduces the
// output size of files and takes care of all the messy calculations we
// would have to do afterwards (local->global coordinate conversion and
// momentum->velocity)
#define dump_tracers(fbase) 		BEGIN_PRIMITIVE{	\
  char dname[256], fname[256] ;					\
  species_t *s = global->tracers_list ;				\
  particle_t *p;						\
  int j;							\
  float pout[13];						\
  double gam;							\
  FileIO f;							\
								\
  sprintf(dname, "%s/T.%d", fbase, step );			\
  dump_mkdir(dname);						\
								\
  while( s ){							\
    if ( s->np > 0 ){						\
      sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank()); \
      f.open(fname,io_write);					\
      WRITE(int,    s->np ,         f ) ;			\
      WRITE(float,  s->q_m ,        f ) ;			\
      WRITE(double, step*grid->dt , f ) ;			\
      WRITE(float,  1.0 ,           f ) ;			\
      WRITE(double, 1.0 ,           f ) ;			\
      p = s->p;							\
      for (j=0 ; j< s->np ; ++j){				\
	gam  = sqrt( 1.0 + p[j].ux * p[j].ux + p[j].uy * p[j].uy\
			 + p[j].uz * p[j].uz );			\
	pout[0] = (float) p[j].q;					\
	pout[1] = (float) tracer_x ;				\
	pout[2] = (float) tracer_y ;				\
	pout[3] = (float) tracer_z ;				\
	pout[4] = (float) p[j].ux;			\
	pout[5] = (float) p[j].uy;			\
	pout[6] = (float) p[j].uz;			\
	pout[7] = (float) field[p->i].ex;			\
	pout[8] = (float) field[p->i].ey;			\
	pout[9] = (float) field[p->i].ez;			\
	pout[10] = (float) field[p->i].cbx;			\
	pout[11] = (float) field[p->i].cby;			\
	pout[12] = (float) field[p->i].cbz;			\
	f.write(pout,13);					\
      }								\
      f.close();						\
    } 								\
    s = s->next;						\
  }								\
} END_PRIMITIVE






#define dump_tracer_restart()		 BEGIN_PRIMITIVE{	\
  species_t *s = global->tracers_list ;				\
  char fname[256] ;						\
  int numspecies = 0;						\
  FileIO f;							\
  sprintf(fname, "restart%d/tracer%d.%d", global->restart_set, global->particle_tracing, 	\
		(int)rank() ) ;					\
  f.open(fname,io_write);					\
								\
  while( s ){							\
    ++numspecies;						\
    s = s->next;						\
  }								\
  WRITE(int, numspecies, f);					\
  s = global->tracers_list;					\
  while( s ){							\
    WRITE_STRING( s->name, f);					\
    WRITE(species_id,   s->id, 			f);		\
    WRITE(float,        s->q_m,			f);		\
    WRITE(int, 		s->max_np,		f);		\
    WRITE(int,		s->max_nm,		f);		\
    WRITE(int,		s->sort_interval, 	f);		\
    WRITE(int,		s->sort_out_of_place,	f);		\
    WRITE(int,		s->np,			f);		\
    if (s->np > 0 ) f.write(s->p, s->np);			\
    s = s->next;						\
  }								\
  f.close();							\
} END_PRIMITIVE

#define read_tracer_restart()  BEGIN_PRIMITIVE{			\
  char fname[256] ;						\
  FileIO f;							\
  int namelen, maxnp, maxnm, sort ;				\
  int sorttype, np, numspecies, i ;				\
  float qm ;							\
  species_id id;						\
  char *name;							\
  species_t * s;						\
								\
  global->tracers_list = NULL ;					\
  sprintf(fname, "restart%d/tracer%d.%d", global->restart_set, global->particle_tracing, \
		(int)rank() ) ;					\
  f.open(fname,io_read);					\
 								\
  READ(int, numspecies, f);					\
  for( i=0 ; i < numspecies ; ++i){				\
    READ(int, namelen, f) ;					\
    name = (char *)malloc(namelen+1) ;				\
    f.read(name, namelen) ;					\
    name[namelen] = '\0' ;					\
    READ(species_id,	id,	  f);				\
    READ(float,		qm,	  f);				\
    READ(int,		maxnp,    f);				\
    READ(int,		maxnm,    f);				\
    READ(int,		sort,     f);				\
    READ(int,		sorttype, f);				\
    READ(int,		np,	  f);				\
    s = new_species(name, qm, maxnp, maxnm, sort, sorttype, 	\
		    &(global->tracers_list));			\
    s->id = id ;						\
    s->np = np ;						\
    if ( np > 0 ) f.read(s->p, np);				\
  }								\
  f.close();							\
} END_PRIMITIVE

// dump tracer by particle trajectory
/*
 * Original 52b output:
 *
 * pout[0]  = step*grid->dt;
 * pout[1]  = (float) tracer_x;
 * pout[2]  = (float) tracer_y;
 * pout[3]  = (float) tracer_z;
 * pout[4]  = (float) p[j].ux;
 * pout[5]  = (float) p[j].uy;
 * pout[6]  = (float) p[j].uz;
 * pout[7]  = (float) field[p->i].ex;
 * pout[8]  = (float) field[p->i].ey;
 * pout[9]  = (float) field[p->i].ez;
 * pout[10] = (float) field[p->i].cbx;
 * pout[11] = (float) field[p->i].cby;
 * pout[12] = (float) field[p->i].cbz;
 */
#define dump_traj(fbase) 		BEGIN_PRIMITIVE{	\
  char dname[256], fname[256];                      \
  char cwd[1024], adname[1024];	                    \
  species_t *s = global->tracers_list ;				\
  particle_t *p;                                    \
  int j;                                            \
  float pout[10];                                   \
  FileIO f;                                         \
  DIR *d;                                           \
                                                    \
  sprintf(dname, "%s", fbase );                     \
  assert(getcwd(cwd, sizeof(cwd)) != NULL);         \
  sprintf(adname, "%s/%s", cwd, dname);             \
  d = opendir(adname);                              \
  dump_mkdir(dname);                                \
                                                    \
  while( s ){                                       \
    if ( s->np > 0 ){                               \
      p = s->p;                                     \
      for (j=0; j<s->np ; ++j){                     \
        int64_t tag = p[j].tag;                     \
        if (tag !=  0 ) {                           \
          sprintf(fname, "%s/%s.%lx",               \
                  dname, s->name, tag);             \
          f.open(fname,io_append);                  \
          pout[0] = step*grid->dt;                  \
          pout[1] = (float) p[j].dx;                \
          pout[2] = (float) p[j].dy;                \
          pout[3] = (float) p[j].dz;                \
          pout[4] = (float) p[j].i;                 \
          pout[5] = (float) p[j].ux;                \
          pout[6] = (float) p[j].uy;                \
          pout[7] = (float) p[j].uz;                \
          memcpy(pout+8, &tag, sizeof(tag));        \
          f.write(pout,10);                         \
          f.close();                                \
        }                                           \
      }                                             \
    }                                               \
    s = s->next;                                    \
  }                                                 \
                                                    \
  closedir(d);                                      \
} END_PRIMITIVE



/****************************************************************
 ******************* Copied from dumpmacros.h *******************
 ***************************************************************/

#define WRITE(type,value,fileIO) BEGIN_PRIMITIVE { \
  type __WRITE_tmp = (type)(value);                \
  fileIO.write( &__WRITE_tmp, 1 );                 \
} END_PRIMITIVE
 
// Note: strlen does not include the terminating NULL
#define WRITE_STRING(string,fileIO) BEGIN_PRIMITIVE {     \
  int __WRITE_STRING_len = 0;                             \
  if( string!=NULL ) __WRITE_STRING_len = strlen(string); \
  fileIO.write( &__WRITE_STRING_len, 1 );                 \
  if( __WRITE_STRING_len>0 )                              \
    fileIO.write( string, __WRITE_STRING_len );           \
} END_PRIMITIVE
 
#define READ(type,value,fileIO) BEGIN_PRIMITIVE { \
  type __READ_tmp;                                \
  fileIO.read(&__READ_tmp, 1 );                   \
  (value) = __READ_tmp;                           \
} END_PRIMITIVE

