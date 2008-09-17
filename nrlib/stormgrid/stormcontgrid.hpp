// $Id$

#ifndef NRLIB_STORMCONTGRID_HPP
#define NRLIB_STORMCONTGRID_HPP

#include <string>

#include "../volume/volume.hpp"
#include "../grid/grid.hpp"

namespace NRLib2 {
  class StormContGrid : public Grid<float>, public Volume {
  public:
    enum FileFormat {STORM_BINARY = 0, STORM_ASCII};

    explicit StormContGrid(const std::string& filename);
    explicit StormContGrid(int nx = 0, int ny = 0, int nz = 0);
    explicit StormContGrid(const Volume &vol, int nx = 0, int ny = 0, int nz = 0);

    void SetMissingCode(double missing_code)
    { missing_code_ = missing_code; }

    double GetMissingCode() const
    { return missing_code_; }

    bool IsMissing(double val) const
    { return val == missing_code_; }

    void SetFormat(FileFormat format) 
    { file_format_ = format; }

    FileFormat GetFormat() const 
    { return file_format_; }

    void SetModelFileName(const std::string& filename) 
    { model_file_name_ = filename; }

    std::string GetModelFileName() const 
    { return model_file_name_; }

    void SetVariableName(const std::string& name) 
    { variable_name_ = name; }

    std::string GetVariableName() const 
    { return variable_name_; }

    void SetZoneNumber(int zone_number) 
    { zone_number_ = zone_number; }

    int GetZoneNumber() const 
    { return zone_number_; }
  
    /// Write to file. If predefinedHeader is not empty, this header is written instead
    ///                of standard, and surfaces are not written.
    void WriteToFile(const std::string& filename, const std::string& predefinedHeader = "") const;
    void ReadFromFile(const std::string& filename, bool commonPath = true);

    double GetDX() const {return GetLX()/GetNI();}
    double GetDY() const {return GetLY()/GetNJ();}
    int FindIndex(double x, double y, double z) const;
    double getValueZInterpolated(double x, double y, double z)const;
    double getValueClosestInZ(double x, double y, double z)const;
  private:
    double RecalculateLZ();

    FileFormat file_format_;
    double missing_code_;
    int zone_number_;
    std::string model_file_name_;
    std::string variable_name_;
};

}

#endif // NRLIB_STORMCONTGRID_HPP
