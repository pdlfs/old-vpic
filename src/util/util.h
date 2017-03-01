#ifndef _util_h_
#define _util_h_

// Expose all public functionality in util.  The below includes bring
// in util_base.h and other low level includes automatically.

#ifdef VPIC_INSTALLED
#include <vpic/mp.h>
#include <vpic/mtrand.h>
#include <vpic/pipelines.h>
#include <vpic/v4.h>
#else
#include "mp/mp.h"
#include "mtrand/mtrand.h"
#include "pipelines/pipelines.h"
#include "v4/v4.h"
#endif

#endif // _util_h_

