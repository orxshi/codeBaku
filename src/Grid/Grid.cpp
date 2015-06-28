#include "Grid.h"

// Constructor
Grid::Grid (string mainDir, int id)
{
    meshFile  = "undefined";    
    bcVerbose[0] = "EMPTY";
    bcVerbose[1] = "SLIP_WALL";
    bcVerbose[2] = "DIRICHLET";
    this->id = id;    
    
    // mainDir
    this->mainDir = mainDir;
    
    // outputDir
    string ids = to_string (id);    
    outputDir = mainDir;
    outputDir.append ("/Grid_");
    outputDir.append (ids);
    outputDir.append ("/");
    
    mkdir (outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);    
    
    // logDir
    logDir = outputDir;
    string temps = "log.dat";
    string slash = "/";
    logDir.append (slash);
    logDir.append (temps);
}

