// $Id: volume.hpp 90 2008-06-19 19:22:52Z perroe $

#ifndef NRLIB_VOLUME_HPP
#define NRLIB_VOLUME_HPP

#include <fstream>

namespace NRLib2 {
  class Surface;

  class Volume {
  public:
    Volume();
    Volume(const Volume& volume);
    ~Volume();
    Volume& operator=(const Volume& rhs);

    void SetDimensions(double x_min, double y_min,
                       double lx, double ly);
    void SetAngle(double angle);

    double GetXMin() const {return x_min_;}
    double GetYMin() const {return y_min_;}
    double GetLX() const {return lx_;}
    double GetLY() const {return ly_;}
    double GetAngle() const {return angle_;}
    
    /// \brief Maximum height of grid.
    double GetLZ() const {return lz_;}

    /// \brief Set surfaces.
    /// \note  The StormVolume object takes control of the pointers, and 
    ///        deletes them when appropriate.
    /// \todo  Is it better to make copies of the surfaces instead???
    void SetSurfaces(Surface* top_surf, 
                     Surface* bot_surf,
                     Surface* erosion_top = 0, 
                     Surface* erosion_bot = 0);

    const Surface& GetTopSurface() const {return *z_top_;}
    const Surface& GetBotSurface() const {return *z_bot_;}
    const Surface& GetErosionTop() const {return *erosion_top_;}
    const Surface& GetErosionBot() const {return *erosion_bot_;}

  protected:
    /// \brief Reader and writer on storm-format.
    /// \todo  Maybe move to storm-specific files.
    void WriteVolumeToFile(std::ofstream& file, 
                           const std::string& filename) const;
    void ReadVolumeFromFile(std::ifstream& file, int line);

    /// \brief The local coorinates are (0,0) in (x_min, y_min), and
    ///        have the same orientation as the volume.
    void GlobalToLocalCoord(double global_x, double global_y,
                            double& local_x, double& local_y) const;

    void LocalToGlobalCoord(double local_x, double local_y,
                            double& global_x, double& global_y) const;

  private:
    virtual double RecalculateLZ();
    /// \brief Checks if surface covers the whole volume.
    bool CheckSurface(const Surface& surface) const;

    double x_min_;
    double y_min_;
    double lx_;
    double ly_;
    double lz_;
    Surface* z_top_;
    Surface* z_bot_;
    Surface* erosion_top_;
    Surface* erosion_bot_;
    double angle_;
  };
} // namespace NRLib2

#endif // NRLIB_VOLUME_HPP
