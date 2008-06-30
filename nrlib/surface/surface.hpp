// $Id: surface.hpp 69 2008-05-23 13:56:27Z perroe $

#ifndef NRLIB_SURFACE_HPP
#define NRLIB_SURFACE_HPP

namespace NRLib2 {
  class Surface {
  public:
    virtual ~Surface();

    /// \brief Generate a copy of the underlying object.
    virtual Surface* Clone() const = 0;

    virtual double GetZ(double x, double y) const = 0;
    
    virtual bool EnclosesRectangle(double x_min, double x_max, 
                                 double y_min, double y_max) const = 0;

    virtual void Add(double c) = 0;

    virtual bool IsMissing(double) const {return(false);}
    virtual bool IsMissing(float) const {return(false);}

    virtual double Min() const = 0;
    virtual double Max() const = 0;
  };

  class ConstantSurface : public Surface {
  public:
    ConstantSurface(double z) : z_(z) {}

    Surface* Clone() const
    { return new ConstantSurface(*this); }

    double GetZ() const {
      return z_;
    }
    
    double GetZ(double /*x*/, double /*y*/) const
    { return z_; }
        
    bool EnclosesRectangle(double /*x_min*/, double /*x_max*/, 
                           double /*y_min*/, double /*y_max*/) const
    { return true; }

    void Add(double c) {
      z_ += c;
    }

    double Min() const {return(z_);}
    double Max() const {return(z_);}

  private:
    double z_;
  };
} // namespace NRLib2

#endif // NRLIB_SURFACE_HPP
