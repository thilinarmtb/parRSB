#ifndef _GENMAP_METRICS_H_
#define _GENMAP_METRICS_H_

#include <genmap-gslib.h>

typedef enum {
  RCB,
  WEIGHTEDLAPLACIANSETUP,
  FIEDLER,
  NFIEDLER,
  FIEDLERSORT,
  BISECTANDREPAIR,
  LANCZOS,
  NLANCZOS,
  WEIGHTEDLAPLACIAN,
  TQLI,
  LAPLACIANSETUP,
  LAPLACIANGSSETUP,
  FINDNBRS,
  FIRSTHALF,
  SECONDHALF,
  CSRMATSETUP,
  CSRTOPSETUP,
  PRECONDSETUP,
  RQI,
  NRQI,
  PROJECT,
  NPROJECT,
  GRAMMIAN,
  LAPLACIAN,
  VCYCLE,
  END
} metric;

void metric_init();
void metric_acc(metric m, double count);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);
void metric_finalize();

#endif
