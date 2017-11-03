/* Defining this reduces produced output */
#define QUIET_RUN
#define LOG_SYSSTAT

/* Define simulation mode: file per process or particle? */
#define VPIC_FILE_PER_PARTICLE 0

/********************************/
/* Define simulation parameters */
/********************************/

/* Node topology. Total nodes = x * y * z */
#define VPIC_TOPOLOGY_X 256
#define VPIC_TOPOLOGY_Y 1
#define VPIC_TOPOLOGY_Z 2

/* Particle distribution. Total particles = 100 * x * y * z */
#define VPIC_PARTICLE_X 2048
#define VPIC_PARTICLE_Y 1
#define VPIC_PARTICLE_Z 1536
