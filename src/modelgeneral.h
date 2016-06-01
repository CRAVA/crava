/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELGENERAL_H
#define MODELGENERAL_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/cravatrend.h"
#include "src/seismicparametersholder.h"
#include "src/state4d.h"
#include "src/timeevolution.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionssolid.h"

#include "nrlib/grid/grid.hpp"
#include "src/blockedlogscommon.h"
#include "src/commondata.h"

struct irapgrid;
class Corr;
class Simbox;
class FFTGrid;
class RandomGen;
class GridMapping;
class InputFiles;
class TimeLine;
class SeismicParameters;
class CravaTrend;

class ModelGeneral
{
public:
  ModelGeneral(ModelSettings           *& modelSettings,
               const InputFiles         * inputFiles,
               SeismicParametersHolder  & seismicParameters,
               CommonData               * common_data,
               int                        i_interval);

  ~ModelGeneral();

  const Simbox                               * GetSimbox()                const { return simbox_                  ;}
  RandomGen                                  * GetRandomGen()             const { return random_gen_              ;}
  GridMapping                                * GetTimeDepthMapping()      const { return time_depth_mapping_      ;}
  CravaTrend                                 & GetTrendCubes()                  { return trend_cubes_             ;}
  CravaTrend                                   GetTrendCubes()            const { return trend_cubes_             ;}
  bool                                         GetVelocityFromInversion() const { return velocity_from_inversion_ ;}
  State4D                                      GetState4D()               const { return state4d_                 ;}
  TimeLine                                   * GetTimeLine()                    { return time_line_               ;}
  std::map<std::string, BlockedLogsCommon *> & GetBlockedWells()                { return blocked_logs_            ;}
  const std::vector<float>                   & GetPriorFacies()       /*const*/ { return prior_facies_            ;}
  const std::vector<FFTGrid *>               & GetPriorFaciesCubes()  /*const*/ { return prior_facies_prob_cubes_ ;}
  const std::vector<std::string>             & GetFaciesNames(void)       const { return facies_names_            ;}
  std::vector<int>                             GetFaciesLabel()           const { return facies_labels_           ;}
  bool                                         GetIs4DActive()            const { return(do_4D_inversion_)        ;}
  std::string                                  GetIntervalName()                { return interval_name_           ;}

  void AddFaciesLabel(int faciesLabel)                                        { facies_labels_.push_back(faciesLabel) ;}
  void AddFaciesName(const std::string & faciesName)                          { facies_names_.push_back(faciesName)   ;}

  std::map<std::string, DistributionsRock *> GetRockDistributionTime0() const;
  const State4D            & getState4D()               const { return state4d_                ;}
  State4D                  * getState4D()                     { return &state4d_               ;}

  static FFTGrid  * CreateFFTGrid(int nx,
                                  int ny,
                                  int nz,
                                  int nxp,
                                  int nyp,
                                  int nzp,
                                  bool fileGrid);

  void              Complete4DBackground(const int nx,const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad,NRLib::Vector &initial_mean,NRLib::Matrix &initial_cov);

  bool              Do4DRockPhysicsInversion(ModelSettings* modelSettings);

  void              MergeCovariance(std::vector<FFTGrid *> & sigma) {state4d_.mergeCov(sigma);}

  void              AdvanceTime(int time_step, SeismicParametersHolder & seismicParameters,ModelSettings* modelSettings);
  void              LastUpdateOfStaticAndDynamicParts(SeismicParametersHolder &  seismicParameters,ModelSettings* modelSettings);
  void              DumpSeismicParameters(ModelSettings* modelSettings, std::string identifyer, int timestep,SeismicParametersHolder &  current_state);

  void              WriteToFile(const Simbox        * simbox,
                                GridMapping         * time_depth_mapping,
                                const ModelSettings * model_settings,
                                FFTGrid             * grid,
                                const std::string   & file_name,
                                const std::string   & sgri_label,
                                bool                  padding = false);

  void              setTimeDepthMapping(GridMapping * new_timeDepthMapping);
  void              Dump4Dparameters(const ModelSettings* modelSettings, std::string identifyer, int timestep,bool printPadding=true);

  void                      mergeState4D(SeismicParametersHolder &  seismicParameters);
  void                      updateState4D(SeismicParametersHolder &  seismicParameters);

  void                      updateState4DWithSingleParameter(FFTGrid * EPost,
                                                             FFTGrid * CovPost,
                                                             int       parameterNumber);

  void                      updateState4DMu(FFTGrid * mu_vp_static,
                                            FFTGrid * mu_vs_static,
                                            FFTGrid * mu_rho_static,
                                            FFTGrid * mu_vp_dynamic,
                                            FFTGrid * mu_vs_dynamic,
                                            FFTGrid * mu_rho_dynamic);


private:

  void              CopyCorrelationsTo4DState(SeismicParametersHolder                    & seismicParameters,
                                              State4D                                    & state4d);

  void              SetupState4D(SeismicParametersHolder & seismicParameters,
                                 const Simbox            * simbox,
                                 State4D                 & state4d,
                                 NRLib::Vector           & initialMean,
                                 NRLib::Matrix           & initialCov);

  const Simbox                                                * simbox_;                       ///< Information about simulation area.

  RandomGen                                                   * random_gen_;                   ///< Random generator.

  CravaTrend                                                    trend_cubes_;                  ///< Trend cubes used in rock phyiscs prior model
  std::map<std::string, std::vector<DistributionsRock *> >      rock_distributions_;           ///< Rocks used in rock physics model
  std::map<std::string, std::vector<DistributionWithTrend *> >  reservoir_variables_;          ///< Reservoir variables used in the rock physics model

  TimeEvolution                                                 time_evolution_;

  GridMapping                                                 * time_depth_mapping_;           ///< Contains both simbox and mapping used for depth conversion

  bool                                                          velocity_from_inversion_;

  TimeLine                                                    * time_line_;

  std::map<std::string, BlockedLogsCommon *>                    blocked_logs_;

  std::vector<float>                                            prior_facies_;                ///< Prior facies probabilities
  std::vector<FFTGrid *>                                        prior_facies_prob_cubes_;     ///< Cubes for prior facies probabilities

  std::vector<int>                                              facies_labels_;               ///< Facies labels, flyttes til blockedlogs
  std::vector<std::string>                                      facies_names_;                ///< Facies names   (nFacies = faciesNames.size()). Use for ordering of facies

  bool                                                          do_4D_inversion_;
  bool                                                          do_4D_rock_physics_vnversion_;
  State4D                                                       state4d_;                     ///< State4D holds the 27 grdis needed for 4D inversion.

  std::string                                                   interval_name_;

};

#endif
