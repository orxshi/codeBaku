#include "Grid.h"

Grid::~Grid ()
{
    cout << "destroying grid" << endl;
    cellADT.destroy_tree();
    cout << "done destroying grid" << endl;
}

// Constructor
Grid::Grid (string mainDir, int id)
{
    cout << "constructing grid" << endl;

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
    
    nHoles = 0;    
    
    cout << "done constructing grid" << endl;
}

Grid::Grid (const Grid& other)
{
    cout << "copy constructing grid" << endl;

    n_bou_elm = other.n_bou_elm;
    n_in_elm = other.n_in_elm;
    totalNElms = other.totalNElms;
    id = other.id;
    phys_count = other.phys_count;
    wallDistance = other.wallDistance;
    //in = other.in;
    meshFile = other.meshFile;
    nHoles = other.nHoles;
    outputDir = other.outputDir;
    mainDir = other.mainDir;
    logDir = other.logDir;
    phys = other.phys;
    bc = other.bc;
    pt = other.pt;
    face = other.face;
    cell = other.cell;
    bt = other.bt;    
    bcVerbose = other.bcVerbose;
    cellADT = other.cellADT;
    holes = other.holes;    
    
    cout << "done copy constructing grid" << endl;
}