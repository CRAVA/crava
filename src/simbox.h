#ifndef SIMBOX_H
#define SIMBOX_H

struct irapgrid;

class Simbox{
public:
  Simbox(void);
  Simbox(double x0, double y0, irapgrid * z0, double lx, double ly, double lz,
         double rot, double dx, double dy, double dz); //Assumes constant thickness.
  Simbox(const Simbox *simbox);
  ~Simbox();

  int        getIndex(double x, double y, double z)   const;
  int        getClosestZIndex(double x, double y, double z);
  void       getIndexes(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void       getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const;
  void       getZInterpolation(double x, double y, double z, 
                              int & index1, int & index2, double & t) const;
  void       getCoord(int xInd, int yInd, int zInd, double &x, double &y, double &z) const;
  int        getnx()                         const { return nx_ ;}
  int        getny()                         const { return ny_ ;}
  int        getnz()                         const { return nz_ ;}
  double     getdx()                         const { return dx_ ;}
  double     getdy()                         const { return dy_ ;}
  double     getdz()                         const { return dz_ ;} // Maximum dz (nz is constant).
  double     getlx()                         const { return lx_ ;}
  double     getly()                         const { return ly_ ;}
  double     getlz()                         const { return lz_ ;} // Maximum thickness
  double     getx0()                         const { return x0_ ;}
  double     gety0()                         const { return y0_ ;}
  double     getAngle()                      const { return rot_ ;}
  int        getIL0()                        const { return inLine0_ ;}
  int        getXL0()                        const { return crossLine0_ ;}
  int        getILStep()                     const { return ilStep_  ;}
  int        getXLStep()                     const { return xlStep_  ;}
  bool       getIsConstantThick()            const { return constThick_;}
  double     getMinRelThick()                const { return minRelThick_;} // Returns minimum relative thickness.
  double     getRelThick(int i, int j)       const;                        // Local relative thickness.
  double     getRelThick(double x, double y) const;                        // Local relative thickness.
  double     getAvgRelThick(void)            const;
  void       getMinMaxZ(double & minZ, double & maxZ) const;
  int        isInside(double x, double y) const;
  int        insideRectangle(double xr, double yr, double rotr, double lxr, double lyr) const;
  double     getTop(double x, double y) const;
  double     getBot(double x, double y) const;
  char     * getStormHeader(int cubetype, int nx, int ny, int nz, bool flat = false, bool ascii = false) const;
  void       writeTopBotGrids(const char * topname, const char * botname);
  int        checkError(double lzLimit, char * errText);
  void       setArea(double x0, double y0, double lx, double ly, double rot, double dx, double dy);
  void       setDepth(irapgrid * zref, double zShift, double lz, double dz);
  void       setDepth(irapgrid * z0, irapgrid * z1, int nz);
  void       setSeisLines(int il0, int cl0, int ilStep, int xlStep);
  int        status() const {return(status_);}
  void       externalFailure() {status_ = EXTERNALERROR;}

  irapgrid * getTopGrid()const {return z0Grid_;};

  enum       simboxstatus{BOXOK, INTERNALERROR, EXTERNALERROR, EMPTY, NOAREA, NODEPTH};

private:
  double     x0_, y0_, lx_, ly_, lz_;  // Simbox is regular, except for top surface.
  irapgrid * z0Grid_;                  // Top surface given as Irapgrid.
  irapgrid * z1Grid_;                  // Bottom surface given as Irapgrid.
  double     rot_;                     // Rotation angle of box.
  double     dx_, dy_, dz_;            // Working resolution.
  int        nx_, ny_, nz_;            // Number of cells in each direction.
  int        status_;                  // Since Simbox may be incomplete or with error
  double     cosrot_, sinrot_;         // Saving time in transformations.
  int        inLine0_, crossLine0_;
  int        ilStep_, xlStep_;
  char     * topName_;
  char     * botName_;
  bool       constThick_;
  double     minRelThick_;
};
#endif
