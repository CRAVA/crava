#include "src/inputfiles.h"

InputFiles::InputFiles(void)
  : seedFile_(""),          
    logNames_(5),         
    wellFiles_(0),           
    seismicFiles_(0),        
    waveletFiles_(0),        
    waveletEstIntFile_(2),  
    faciesEstIntFile_(2),   
    timeSurfFiles_(0),       
    depthSurfFiles_(2),      
    velocityField_(""),     
    backFile_(3),           
    backVelFile_(""),       
    corrDirFile_(""),       
    reflMatrFile_(""),      
    paramCorrFile_("")
{
}

InputFiles::~InputFiles(void)
{
}

