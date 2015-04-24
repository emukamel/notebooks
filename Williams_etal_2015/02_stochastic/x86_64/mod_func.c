#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _h_reg(void);
extern void _kadist_reg(void);
extern void _kaprox_reg(void);
extern void _kdrca1_reg(void);
extern void _km_reg(void);
extern void _na3n_reg(void);
extern void _naxn_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," h.mod");
    fprintf(stderr," kadist.mod");
    fprintf(stderr," kaprox.mod");
    fprintf(stderr," kdrca1.mod");
    fprintf(stderr," km.mod");
    fprintf(stderr," na3n.mod");
    fprintf(stderr," naxn.mod");
    fprintf(stderr," vecevent.mod");
    fprintf(stderr, "\n");
  }
  _h_reg();
  _kadist_reg();
  _kaprox_reg();
  _kdrca1_reg();
  _km_reg();
  _na3n_reg();
  _naxn_reg();
  _vecevent_reg();
}
