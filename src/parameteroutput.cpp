/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/parameteroutput.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
#include "src/gridmapping.h"
#include "src/io.h"
#include "lib/utils.h"
#include "fft/include/fftw.h"

void
ParameterOutput::WriteParameters(const Simbox        * simbox,
                                 GridMapping         * time_depth_mapping,
                                 const ModelSettings * model_settings,
                                 StormContGrid       * vp,
                                 StormContGrid       * vs,
                                 StormContGrid       * rho,
                                 int                   output_flag,
                                 int                   sim_num,
                                 bool                  kriged)
{
  std::string prefix;
  std::string suffix;
  std::string file_name;

  if(sim_num >= 0) {
    prefix = IO::PrefixSimulations();
    suffix = "_"+NRLib::ToString(sim_num+1);
  }
  else {
    prefix = IO::PrefixPredictions();
    suffix = "";
  }

  if(kriged)
    suffix = "_Kriged"+suffix;

  if((output_flag & IO::MURHO) > 0) {
    file_name = prefix+"MuRho"+suffix;
    ComputeMuRho(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMBDARHO) > 0) {
    file_name = prefix+"LambdaRho"+suffix;
    ComputeLambdaRho(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMELAMBDA) > 0) {
    file_name = prefix+"LameLambda"+suffix;
    ComputeLameLambda(simbox, time_depth_mapping, model_settings, vp, vs, rho, file_name);
  }
  if((output_flag & IO::LAMEMU) > 0) {
    file_name = prefix+"LameMu"+suffix;
    ComputeLameMu(simbox, time_depth_mapping,  model_settings, vs, rho, file_name);
  }
  if((output_flag & IO::POISSONRATIO) > 0) {
    file_name = prefix+"PoissonRatio"+suffix;
    ComputePoissonRatio(simbox, time_depth_mapping, model_settings, vp, vs, file_name);
  }
  if((output_flag & IO::AI) > 0) {
    file_name = prefix+"AI"+suffix;
    ComputeAcousticImpedance(simbox, time_depth_mapping, model_settings, vp, rho, file_name);
  }
  if((output_flag & IO::SI) > 0) {
    file_name = prefix+"SI"+suffix;
    ComputeShearImpedance(simbox, time_depth_mapping, model_settings, vs, rho, file_name);
  }
  if((output_flag & IO::VPVSRATIO) > 0) {
    file_name = prefix+"VpVsRatio"+suffix;
    ComputeVpVsRatio(simbox, time_depth_mapping, model_settings, vp, vs, file_name);
  }
  if((output_flag & IO::VP) > 0) {
    file_name = prefix+"Vp"+suffix;

    ExpTransf(vp);
    WriteToFile(simbox, time_depth_mapping, model_settings, vp, file_name, "Inverted Vp");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  vp->logTransf();

  }
  if((output_flag & IO::VS) > 0) {
    file_name = prefix+"Vs"+suffix;

    ExpTransf(vs);
    WriteToFile(simbox, time_depth_mapping, model_settings, vs, file_name, "Inverted Vs");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  vs->logTransf();

  }
  if((output_flag & IO::RHO) > 0) {
    file_name = prefix+"Rho"+suffix;

    ExpTransf(rho);
    WriteToFile(simbox, time_depth_mapping, model_settings, rho, file_name, "Inverted density");
    //if (sim_num < 0) //prediction, need grid unharmed.
    //  rho->logTransf();

  }
}

void
ParameterOutput::ComputeAcousticImpedance(const Simbox        * simbox,
                                          GridMapping         * time_depth_mapping,
                                          const ModelSettings * model_settings,
                                          StormContGrid       * vp,
                                          StormContGrid       * rho,
                                          const std::string   & file_name)
{
  StormContGrid * pr_impedance = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = exp(ijk_a + ijk_r);

        pr_impedance->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, pr_impedance, file_name, "Acoustic Impedance");

  delete pr_impedance;
}


void
ParameterOutput::ComputeShearImpedance(const Simbox        * simbox,
                                       GridMapping         * time_depth_mapping,
                                       const ModelSettings * model_settings,
                                       StormContGrid       * vs,
                                       StormContGrid       * rho,
                                       const std::string   & file_name)
{
  StormContGrid * sh_impedance = new StormContGrid(*vs);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vs->GetNI(); i++) {
    for (size_t j = 0; j < vs->GetNJ(); j++) {
      for (size_t k = 0; k < vs->GetNK(); k++) {

        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = exp(ijk_b + ijk_r);

        sh_impedance->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, sh_impedance, file_name, "Shear impedance");

  delete sh_impedance;
}


void
ParameterOutput::ComputeVpVsRatio(const Simbox        * simbox,
                                  GridMapping         * time_depth_mapping,
                                  const ModelSettings * model_settings,
                                  StormContGrid       * vp,
                                  StormContGrid       * vs,
                                  const std::string   & file_name)
{
  StormContGrid * ratio_vp_vs = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        comp_val = exp(ijk_a - ijk_b);

        ratio_vp_vs->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, ratio_vp_vs, file_name, "Vp-Vs ratio");

  delete ratio_vp_vs;
}

void
ParameterOutput::ComputePoissonRatio(const Simbox        * simbox,
                                     GridMapping         * time_depth_mapping,
                                     const ModelSettings * model_settings,
                                     StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     const std::string   & file_name)
{
  StormContGrid * poi_rat = new StormContGrid(*vp);

  float ijk_a      = 0.0f;
  float ijk_b      = 0.0f;
  float comp_val   = 0.0f;
  float v_ratio_sq = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a      = vp->GetValue(i, j, k);
        ijk_b      = vs->GetValue(i, j, k);
        v_ratio_sq = exp(2*(ijk_a-ijk_b));
        comp_val   = static_cast<float>(0.5*(v_ratio_sq - 2)/(v_ratio_sq - 1));

        poi_rat->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, poi_rat, file_name, "Poisson ratio");

  delete poi_rat;
}

void
ParameterOutput::ComputeLameMu(const Simbox        * simbox,
                               GridMapping         * time_depth_mapping,
                               const ModelSettings * model_settings,
                               StormContGrid       * vs,
                               StormContGrid       * rho,
                               const std::string   & file_name)
{
  StormContGrid * mu = new StormContGrid(*vs);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vs->GetNI(); i++) {
    for (size_t j = 0; j < vs->GetNJ(); j++) {
      for (size_t k = 0; k < vs->GetNK(); k++) {

        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(ijk_r+2*ijk_b-13.81551)); // -13.81551 in the exponent divides by 1 000 000

        mu->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, mu, file_name, "Lame mu");

  delete mu;
}

void
ParameterOutput::ComputeLameLambda(const Simbox        * simbox,
                                   GridMapping         * time_depth_mapping,
                                   const ModelSettings * model_settings,
                                   StormContGrid       * vp,
                                   StormContGrid       * vs,
                                   StormContGrid       * rho,
                                   const std::string   & file_name)
{
  StormContGrid * lambda = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(ijk_r)*(exp(2*ijk_a-13.81551)-2*exp(2*ijk_b-13.81551))); // -13.81551 in the exponent divides by 1 000 000

        lambda->SetValue(i, j, k, comp_val);
      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, lambda, file_name, "Lame lambda");

  delete lambda;
}

void
ParameterOutput::ComputeLambdaRho(const Simbox        * simbox,
                                  GridMapping         * time_depth_mapping,
                                  const ModelSettings * model_settings,
                                  StormContGrid       * vp,
                                  StormContGrid       * vs,
                                  StormContGrid       * rho,
                                  const std::string   & file_name)
{
  StormContGrid * lambda_rho = new StormContGrid(*vp);

  float ijk_a    = 0.0f;
  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;

  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_a    = vp->GetValue(i, j, k);
        ijk_b    = vs->GetValue(i, j, k);
        ijk_r    = rho->GetValue(i, j, k);
        comp_val = static_cast<float>(exp(2.0*(ijk_a +ijk_r)-13.81551)-2.0*exp(2.0*(ijk_b +ijk_r)-13.81551)); // -13.81551 in the exponent divides by 1e6=(1 000 000)

        lambda_rho->SetValue(i, j, k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, lambda_rho, file_name, "Lambda rho");

  delete lambda_rho;
}

void
ParameterOutput::ComputeMuRho(const Simbox        * simbox,
                              GridMapping         * time_depth_mapping,
                              const ModelSettings * model_settings,
                              StormContGrid       * vp,
                              StormContGrid       * vs,
                              StormContGrid       * rho,
                              const std::string   & file_name)
{
  StormContGrid * mu_rho;
  mu_rho = new StormContGrid(*vp);

  float ijk_b    = 0.0f;
  float ijk_r    = 0.0f;
  float comp_val = 0.0f;
  for (size_t i = 0; i < vp->GetNI(); i++) {
    for (size_t j = 0; j < vp->GetNJ(); j++) {
      for (size_t k = 0; k < vp->GetNK(); k++) {

        ijk_b    = vs->GetValue(i,j,k);
        ijk_r    = rho->GetValue(i,j,k);
        comp_val = static_cast<float>(exp(2.0*(ijk_b +ijk_r)-13.81551)); // -13.81551 in the exponent divides by 1e6=(1 000 000)

        mu_rho->SetValue(i,j,k, comp_val);

      }
    }
  }

  WriteToFile(simbox, time_depth_mapping, model_settings, mu_rho, file_name, "Mu rho");

  delete mu_rho;
}

//FFTGrid*
//ParameterOutput::createFFTGrid(FFTGrid * referenceGrid, bool file_grid)
//{
//  int nx  = referenceGrid->getNx();
//  int ny  = referenceGrid->getNy();
//  int nz  = referenceGrid->getNz();
//  int nxp = referenceGrid->getNxp();
//  int nyp = referenceGrid->getNyp();
//  int nzp = referenceGrid->getNzp();
//
//  FFTGrid * fftGrid;
//
//  if(file_grid)
//    fftGrid = new FFTFileGrid(nx,ny,nz,nxp,nyp,nzp);
//  else
//    fftGrid = new FFTGrid(nx,ny,nz,nxp,nyp,nzp);
//
//  return(fftGrid);
//}

//void
//ParameterOutput::WriteToFile(const Simbox        * simbox,
//                             GridMapping         * time_depth_mapping,
//                             const ModelSettings * model_settings,
//                             FFTGrid             * grid,
//                             const std::string   & file_name,
//                             const std::string   & sgri_label,
//                             bool                  padding)
//{
//  //GridMapping * timeDepthMapping = modelGeneral->GetTimeDepthMapping();
//  //GridMapping * timeCutMapping;//   = modelGeneral->getTimeCutMapping(); //Included in the new simbox format.
//  float seismic_start_time  = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
//  TraceHeaderFormat *format = model_settings->getTraceHeaderFormatOutput();
//
//  grid->writeFile(file_name,
//                  IO::PathToInversionResults(),
//                  simbox,
//                  sgri_label,
//                  seismic_start_time,
//                  time_depth_mapping,
//                  *format,
//                  padding);
//}

void
ParameterOutput::WriteToFile(const Simbox        * simbox,
                             GridMapping         * time_depth_mapping,
                             const ModelSettings * model_settings,
                             StormContGrid       * grid,
                             const std::string   & file_name,
                             const std::string   & sgri_label,
                             bool                  padding)
{

  float seismic_start_time   = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
  TraceHeaderFormat * format = model_settings->getTraceHeaderFormatOutput();

  WriteFile(model_settings,
            grid,
            file_name,
            IO::PathToInversionResults(),
            simbox,
            false,
            sgri_label,
            seismic_start_time,
            time_depth_mapping,
            *format,
            padding);
}

void
ParameterOutput::ExpTransf(StormContGrid * grid)
{
  float value = 0.0f;
  for (size_t i = 0; i < grid->GetNI(); i++) {
    for (size_t j = 0; j < grid->GetNJ(); j++) {
      for (size_t k = 0; k < grid->GetNK(); k++) {

        value = grid->GetValue(i, j, k);

        if (value != RMISSING) {
          value = exp(value);
          grid->SetValue(i, j, k, value);
        }

      }
    }
  }
}

void
ParameterOutput::WriteFile(const ModelSettings     * model_settings,
                           StormContGrid           * storm_grid,
                           const std::string       & f_name,
                           const std::string       & sub_dir,
                           const Simbox            * simbox,
                           bool                      is_seismic,
                           const std::string         label,
                           const float               z0,
                           const GridMapping       * depth_map,
                           const TraceHeaderFormat & thf,
                           bool                      padding)
{
  //All crava files are written out directly in CravaResult
  (void) padding;
  std::string file_name = IO::makeFullFileName(sub_dir, f_name);
  int format_flag       = model_settings->getOutputGridFormat();
  int domain_flag       = model_settings->getOutputGridDomain();

  if (format_flag > 0) {//Output format specified.
    StormContGrid * output = NULL;
    if(is_seismic == true) {
      output = new StormContGrid(*storm_grid);
      SeismicShift(output);
    }
    else
      output = storm_grid;

    NRLib::Volume * output_vol = output;
    *output_vol = *simbox;

    if ((domain_flag & IO::TIMEDOMAIN) > 0) {
      if ((format_flag & IO::STORM) > 0) {
        const std::string header = simbox->getStormHeader(1, simbox->getnx(), simbox->getny(), simbox->getnz(), false, false);
        output->SetFormat(NRLib::StormContGrid::STORM_BINARY);
        std::string file_name_storm = file_name + IO::SuffixStormBinary();
        LogKit::LogFormatted(LogKit::Low," Writing STORM file "+file_name_storm+"...");
        output->WriteToFile(file_name_storm, header, false);
        LogKit::LogFormatted(LogKit::Low,"done\n");
      }

      if ((format_flag & IO::ASCII) > 0) {
        output->SetFormat(NRLib::StormContGrid::STORM_ASCII);
        const std::string header = simbox->getStormHeader(1, simbox->getnx(), simbox->getny(), simbox->getnz(), false, true);
        std::string file_name_ascii = file_name + IO::SuffixGeneralData();
        LogKit::LogFormatted(LogKit::Low," Writing ASCII file "+file_name_ascii+"...");
        output->WriteToFile(file_name_ascii, header, true);
        LogKit::LogFormatted(LogKit::Low,"done\n");
      }

      //SEGY, SGRI CRAVA are never resampled in time.
      if ((format_flag & IO::SEGY) > 0) {
        std::string file_name_segy = file_name + IO::SuffixSegy();
        LogKit::LogFormatted(LogKit::Low," Writing SEGY file "+file_name_segy+"...");

        SegyGeometry geometry(simbox->getx0(), simbox->gety0(), simbox->getdx(), simbox->getdy(),
                              simbox->getnx(), simbox->getny(),simbox->getIL0(), simbox->getXL0(),
                              simbox->getILStepX(), simbox->getILStepY(),
                              simbox->getXLStepX(), simbox->getXLStepY(),
                              simbox->getAngle());
        SegY * segy = new SegY(output,
                               &geometry,
                               z0,
                               file_name_segy,
                               true, //Write to file
                               thf,
                               is_seismic);

        LogKit::LogFormatted(LogKit::Low,"done\n");

        delete segy;
      }
      if ((format_flag & IO::SGRI) >0) {
        std::string file_name_sgri   = file_name + IO::SuffixSgri();
        std::string file_name_header = file_name + IO::SuffixSgriHeader();

        LogKit::LogFormatted(LogKit::Low," Writing SGRI header file "+ file_name_header + "...");
        output->WriteToSgriFile(file_name_sgri, file_name_header, label, simbox->getdz());
        LogKit::LogFormatted(LogKit::Low,"done\n");

      }
    }

    if (depth_map != NULL && (domain_flag & IO::DEPTHDOMAIN) > 0) { //Writing in depth. Currently, only stormfiles are written in depth.
      std::string depth_name = file_name+"_Depth";

      if (depth_map->getMapping() == NULL) {
        if (depth_map->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::Warning,
            "WARNING: Depth interval lacking when trying to write %s. Write cancelled.\n",depth_name.c_str());
          return;
        }
        if ((format_flag & IO::STORM) > 0) {
          output->SetFormat(NRLib::StormContGrid::STORM_BINARY);
          std::string file_name_storm = depth_name + IO::SuffixStormBinary();
          int nx = static_cast<int>(output->GetNI());
          int ny = static_cast<int>(output->GetNJ());
          int nz = static_cast<int>(output->GetNK());
          std::string header = depth_map->getSimbox()->getStormHeader(FFTGrid::PARAMETER, nx, ny, nz, false, false);
          LogKit::LogFormatted(LogKit::Low," Writing STORM file "+file_name_storm+"...");
          output->WriteToFile(file_name_storm, header, false);
          LogKit::LogFormatted(LogKit::Low,"done\n");
        }
        if ((format_flag & IO::ASCII) > 0) {
          output->SetFormat(NRLib::StormContGrid::STORM_ASCII);
          std::string file_name_ascii = depth_name + IO::SuffixGeneralData();
          int nx = static_cast<int>(output->GetNI());
          int ny = static_cast<int>(output->GetNJ());
          int nz = static_cast<int>(output->GetNK());
          LogKit::LogFormatted(LogKit::Low," Writing ASCII file "+file_name_ascii+"...");
          std::string header = depth_map->getSimbox()->getStormHeader(FFTGrid::PARAMETER, nx, ny, nz, false, true);
          output->WriteToFile(file_name_ascii, header, true);
          LogKit::LogFormatted(LogKit::Low,"done\n");
        }
        if ((format_flag & IO::SEGY) >0) {
          //makeDepthCubeForSegy(depth_map->getSimbox(),depth_name);
          StormContGrid * storm_cube_depth = new StormContGrid(*depth_map->getSimbox(), output->GetNI(), output->GetNJ(), output->GetNK());

          for (size_t i = 0; i < output->GetNI(); i++) {
            for (size_t j = 0; j < output->GetNJ(); j++) {
              for (size_t k = 0; k < output->GetNK(); k++) {
                (*storm_cube_depth)(i, j, k) = output->GetValue(i, j, k);
              }
            }
          }

          SegyGeometry geometry(simbox->getx0(), simbox->gety0(), simbox->getdx(), simbox->getdy(),
                                simbox->getnx(), simbox->getny(),simbox->getIL0(), simbox->getXL0(),
                                simbox->getILStepX(), simbox->getILStepY(),
                                simbox->getXLStepX(), simbox->getXLStepY(),
                                simbox->getAngle());

          std::string file_name_segy = file_name + IO::SuffixSegy();

          LogKit::LogFormatted(LogKit::Low," Writing SEGY file "+file_name_segy+"...");
          SegY * segy = new SegY(storm_cube_depth, &geometry, z0, file_name_segy, true);
          delete segy;
          delete storm_cube_depth;
          LogKit::LogFormatted(LogKit::Low,"done\n");
        }
      }
      else {
        if (depth_map->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::Warning,
            "WARNING: Depth mapping incomplete when trying to write %s. Write cancelled.\n",depth_name.c_str());
          return;
        }
        // Writes also segy in depth if required
        WriteResampledStormCube(output, depth_map, depth_name, simbox, format_flag, z0);
      }
    }

    //If seismic, the grid is copied
    if (is_seismic == true) {
      if (output != NULL)
        delete output;
    }

  }
}

void
ParameterOutput::WriteResampledStormCube(const StormContGrid * storm_grid,
                                         const GridMapping   * gridmapping,
                                         const std::string   & file_name,
                                         const Simbox        * simbox,
                                         const int             format,
                                         float                 z0)
{
  // simbox is related to the cube we resample from. gridmapping contains simbox for the cube we resample to.

  float time, kindex;
  StormContGrid *mapping = gridmapping->getMapping();
  StormContGrid *outgrid = new StormContGrid(*mapping);

  double x,y;
  int nz = static_cast<int>(mapping->GetNK());
  for (int i = 0; i < static_cast<int>(storm_grid->GetNI()); i++) {
    for (int j = 0; j < static_cast<int>(storm_grid->GetNJ()); j++) {
      simbox->getXYCoord(i,j,x,y);
      for (int k = 0; k < nz; k++) {
        time   = (*mapping)(i,j,k);
        kindex = float((time - static_cast<float>(simbox->getTop(x,y)))/simbox->getdz());

        float value = storm_grid->GetValueInterpolated(i, j, kindex);

        (*outgrid)(i,j,k) = value;
      }
    }
  }

  std::string gf_name;
  std::string header;
  if ((format & IO::ASCII) > 0) { // ASCII
    gf_name = file_name + IO::SuffixGeneralData();
    int nx = static_cast<int>(storm_grid->GetNI());
    int ny = static_cast<int>(storm_grid->GetNJ());
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER, nx, ny, nz, 0, 1);
    outgrid->SetFormat(StormContGrid::STORM_ASCII);
    LogKit::LogFormatted(LogKit::Low," Writing ASCII file "+gf_name+"...");
    outgrid->WriteToFile(gf_name, header);
    LogKit::LogFormatted(LogKit::Low,"done\n");
  }

  if ((format & IO::STORM) > 0) {
    gf_name =  file_name + IO::SuffixStormBinary();
    int nx = static_cast<int>(storm_grid->GetNI());
    int ny = static_cast<int>(storm_grid->GetNJ());
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER, nx, ny, nz, 0, 0);
    outgrid->SetFormat(StormContGrid::STORM_BINARY);
    LogKit::LogFormatted(LogKit::Low," Writing STORM file "+gf_name+"...");
    outgrid->WriteToFile(gf_name,header);
    LogKit::LogFormatted(LogKit::Low,"done\n");
  }
  if((format & IO::SEGY) > 0) {
    gf_name =  file_name + IO::SuffixSegy();

    SegyGeometry geometry(simbox->getx0(), simbox->gety0(), simbox->getdx(), simbox->getdy(),
                          simbox->getnx(), simbox->getny(),simbox->getIL0(), simbox->getXL0(),
                          simbox->getILStepX(), simbox->getILStepY(),
                          simbox->getXLStepX(), simbox->getXLStepY(),
                          simbox->getAngle());
    LogKit::LogFormatted(LogKit::Low," Writing SEGY file "+gf_name+"...");
    SegY * segy = new SegY(outgrid, &geometry, z0, gf_name, true);
    delete segy;
    LogKit::LogFormatted(LogKit::Low,"done\n");

    //writeSegyFromStorm(outgrid, gf_name);
  }
  delete outgrid;
}


void
ParameterOutput::SeismicShift(NRLib::Grid<float> * grid)
{
  for(size_t i=0;i<grid->GetNI();i++) {
    for(size_t j=0;j<grid->GetNJ();j++) {
      std::vector<fftw_real> trace(grid->GetNK());
      for(size_t k=0;k<grid->GetNK();k++)
        trace[k] = (*grid)(i,j,k);
      Utils::ShiftTrace(trace);
      for(size_t k=0;k<grid->GetNK();k++)
        (*grid)(i,j,k) = trace[k];
    }
  }
}

