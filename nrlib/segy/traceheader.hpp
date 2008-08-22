// $Id$
#ifndef TRACEHEADER_HPP
#define TRACEHEADER_HPP

#include <cstdio>

#include <sstream>
#include <iostream>
/**
  The format of the trace header. Specifies the location 
  of the fields of interest. All locations are 1-based, i.e.
  the location of the first field in the header is 1.
  The locations are set to -1 if the fields are not given in the header.
 */
class TraceHeaderFormat {
public:
  /**
     Possible coordinate systems.
   */
  enum coordSys_t {
    UTM  = 0,   ///< UTM coordinates.
    ILXL = 1    ///< Inline/Crossline coordinates.
  };
  enum headers {
    SEISWORKS = 1,
    IESX = 2
  };
  /**
    Default constructor, the standard format.
   */
  TraceHeaderFormat(int headerformat);

  /**
    Constructor. Locations should be set to -1 if not set.
    @param[in] utmxLoc      Location scaling coefficient for UTM X and Y.
    @param[in] utmxLoc      Location of the UTM-X field.
    @param[in] utmyLoc      Location of the UTM-Y field.
    @param[in] inlineLoc    Location of the inline coordinate field.
    @param[in] crosslineLoc Location of crossline coordinate field.
   */
  TraceHeaderFormat(int scaleCoLoc,
                    int utmxLoc,
                    int utmyLoc,
                    int inlineLoc,
                    int crosslineLoc,
                    coordSys_t coordSys);

 
  void bypassCoordinateScaling();

  /// Get location of the UTM-X field. (-1 if non-existant)  
  int getUtmxLoc() const {return utmxLoc_;}

  /// Get location of the UTM-Y field. (-1 if non-existant)  
  int getUtmyLoc() const {return utmyLoc_;}

  /// Get location of the inline field. (-1 if non-existant)  
  int getInlineLoc() const {return inlineLoc_;}

  /// Get location of the crossline field. (-1 if non-existant)  
  int getCrosslineLoc() const {return crosslineLoc_;}

  /// Get location of the scaling cooefficient field. (-1 if non-existant)
  int getScalCoLoc() const {return scalCoLoc_;}

  /// Get coordinate system.
  coordSys_t getCoordSys() const {return coordSys_;}

  /// String representation.
  const char* toString() const;
private:
  /// Location of scaling coefficient. (-1 if non-existant)
  int scalCoLoc_;
  /// Location of the UTM-X field. (-1 if non-existant)
  int utmxLoc_;
  /// Location of the UTM-Y field. (-1 if non-existant)
  int utmyLoc_;
  /// Location of inline coordinate field. (-1 if non-existant)
  int inlineLoc_;
  /// Location of crossline coordinate field. (-1 if non-existant)
  int crosslineLoc_;
  /// Coordinate system to use.
  coordSys_t coordSys_;
};


/**
  The trace header for a SEGY file.
 */
class TraceHeader {
public:
  /**
    Constructor generating an empty header.
    @param[in] format  header format.
   */
  TraceHeader(const TraceHeaderFormat& format = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));

  /**
    Read in a new header.
    @param[in] inFile  input file.
    @param[in] lineNo  line number. (from binary header.) -1 if not used.
   */
  void read(std::istream& inFile,
            int lineNo = -1);

  /**
    Write header to file.
    @param[in]  outFile output file.
   */
  int write(std::ostream& outFile);

  /**
    Get UTM X coordinate.
    returns #RMISSING if UTM X location is not set in the format.
   */
  float getUtmx() const;

  /**
    Set UTM X coordinate.
    Does nothing if UTM X location is not set in the format.
   */
  void setUtmx(float utmx);

  /**
    Get UTM Y coordinate.
    returns #RMISSING if UTM Y location is not set in the format.
   */
  float getUtmy() const;

  /**
    Set UTM Y coordinate.
    Does nothing if UTM Y location is not set in the format.
  */
  void setUtmy(float utmy);

  /**
    Get inline coordinate.
    returns #IMISSING if inline location is not set in the format.
   */
  int getInline() const;

  /**
    Set inline coordinate.
    Does nothing if inline location is not set in the format.
   */ 
  void setInline(int inLine);

  /**
    Get crossline coordinate.
    returns #IMISSING if crossline location is not set in the format.
   */
  int getCrossline() const;

  /**
    Set crossline coordinate.
    Does nothing if crossline location is not set in the format.
   */ 
  void setCrossline(int crossLine);

  /**
    Set number of samples.
   */
  void setNSamples(int ns);

  /**
    Set sample interval in mikroseconds.
  */
  void setDt(int dt);
  short getDt() {return(dt_);}

  /**
    Get status code.
    0  - everything went OK.
    -1 - not a trace header, but EDBDIC header.
    -2 - error reading from file.
   */
  int getStatus() const {return status_;}
private:
  /// Header buffer in machine-specific byte order. 
  char buffer_[240];
  /// string with the conetent of buffer_
 // std::istringstream ist_;
  /// Header format.
  TraceHeaderFormat format_;

  /// Status code.
  int status_;

  /// Scaling coefficient for UTM X and UTM Y.
  float scalCo_;

  short scalcoinitial_;
  float utmx_;
  float utmy_;
  int inline_;
  int crossline_;
  short ns_;
  short dt_;
  int imissing_;
  float rmissing_;

 bool useBinaryInline;
  /**
    Get scaling coefficient for SX and SY from buffer.
   */
  short getScalCo() const;

  /**
    Get 32-bit value from buffer.
    @param[in] loc  Location in buffer (1-based).
    @return value from buffer.
   */
  int getInt32(int loc) const;

  /**
    Set 32-bit value in buffer.
    @param[in] loc  Location in buffer (1-based).
    @param[in] val  Value.
   */
  void setInt32(int loc, int val);

  /**
    Get 16-bit value from buffer.
    @param[in] loc  Location in buffer (1-based).
    @return value from buffer.
   */
  short getInt16(int loc) const;

  /**
    Set 16-bit value in buffer.
    @param[in] loc  Location in buffer (1-based).
    @param[in] val  Value.
   */
  void setInt16(int loc, short val);

 
};

#endif
