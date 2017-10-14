
#ifndef incl_Definitions_h
#define incl_Definitions_h


#define DOMAIN_TYPE_NAMES {"Domain",             \
	                   "Geometry",           \
	                   "Mesh",               \
	                   "Discrete",           \
	                   "FiniteElementBVP",   \
	                   "FiniteElementBVPWI", \
	                   "FiniteElementBVPWNI",\
	                   "Fluid",              \
	                   "Solid",              \
	                   "Interface",          \
	                   "MicroCellWulf",      \
	                   "MicroCell",          \
                           "Fem1D",              \
	                   "BeamSection",        \
	                   "BeamBending",        \
                           "RigidBody",          \
                           "Aeroplane",          \
                           "InterfaceMatch",     \
                           "FreeSurface",        \
                           "InterfaceN",         \
                           "InterfaceGS",        \
                           "InterfacePB",        \
                           "VortexSheet",        \
                           "LiftingLine",        \
                           "FlexibleWing",       \
                           "SolidALE",           \
                           "IsogeometricFEM",    \
                           "HBSplineFEM",        \
                           "HBSplineCutFEM",     \
                           "StandardFEM",        \
                           "HBScutFEMElasticity",        NULL}

#define DOMAIN_KEY  {"DOMAIN",             \
	             "GEOMETRY",           \
	             "MESH",               \
	             "DISCRETE",           \
	             "FINITEELEMENTBVP",   \
	             "FINITEELEMENTBVPWI", \
	             "FINITEELEMENTBVPWNI",\
		     "FLUID",              \
		     "SOLID",              \
		     "INTERFACE",          \
	             "MICROCELLWULF",      \
	             "MICROCELL",          \
	             "FEM1D",              \
	             "BEAMSECTION",        \
	             "BEAMBENDING",        \
	             "RIGIDBODY",          \
                     "AEROPLANE",          \
                     "INTERFACEMATCH",     \
                     "FREESURFACE",        \
                     "INTERFACEN",         \
                     "INTERFACEGS",        \
                     "INTERFACEPB",        \
                     "VORTEXSHEET",        \
                     "LIFTINGLINE",        \
                     "FLEXIBLEWING",       \
                     "SOLIDALE",           \
                     "ISOGEOMETRICFEM",    \
                     "HBSPLINEFEM",        \
                     "HBSPLINECUTFEM",     \
                     "STANDARDFEM",        \
                     "HBSCUTFEMELASTICITY",        NULL}


#define DOMAIN_TYPE_ENUM { ROOTDOMAIN,         \
	                   GEOMETRY,           \
	                   MESH,               \
		           DISCRETE,           \
	                   FINITEELEMENTBVP,   \
	                   FINITEELEMENTBVPWI, \
	                   FINITEELEMENTBVPWNI,\
		           FLUID,              \
		           SOLID,              \
		           INTERFACE,          \
                           MICROCELLWULF,      \
                           MICROCELL,          \
	                   FEM1D,              \
                           BEAMSECTION,        \
                           BEAMBENDING,        \
                           RIGIDBODY,          \
                           AEROPLANE,          \
                           INTERFACEMATCH,     \
                           FREESURFACE,        \
                           INTERFACEN,         \
                           INTERFACEGS,        \
                           INTERFACEPB,        \
                           VORTEXSHEET,        \
                           LIFTINGLINE,        \
                           FLEXIBLEWING,       \
                           SOLIDALE,           \
                           ISOGEOMETRICFEM,    \
                           HBSPLINEFEM,        \
                           HBSPLINECUTFEM,     \
                           STANDARDFEM,        \
                           HBSCUTFEMELASTICITY }


#define ELEMENT_GEOMETRY_2D_ENUM { TRIANGLE3, \
	                           TRIANGLE6, \
	                           QUAD4,     \
                                   QUAD8,     \
                                   QUAD9   }


#define GEOMETRY_TYPE_ENUM { POINT,   \
	                     SPLINE,  \
	                     SURFACE }


#define PROPERTY_TYPE_ENUM { ELEMENTTYPE, \
	                     MATERIAL,    \
	                     ALETYPE }


#define ELEMENT_TYPE_NAMES {"2D3nodedTriangle",              \
	                    "2D6nodedTriangle",              \
	                    "2D4nodedQuadrilateral",         \
	                    "2D8nodedQuadrilateral",         \
	                    "2D9nodedQuadrilateral",         \
                            "3D4nodedTetrahedron",           \
                            "3D10nodedTetrahedron",          \
                            "3D8nodedBrick",                 \
                            "3D20nodedBrick",                \
	                    "2D3nodedStabIncompFluid",       \
	                    "2D3nodedStabIncompHighReFluid", \
	                    "2D3nodedStabCompFluid",         \
	                    "2D3nodedLinearSolid",           \
	                    "2D6nodedQuadraticSolid",        \
	                    "2D4nodedLinearSolid",           \
	                    "2D4nodedFBarSolid",             \
	                    "2D8nodedQuadraticSolid",        \
	                    "2D9nodedQuadraticSolid",        \
	                    "3D4nodedStabIncompFluid",       \
	                    "3D4nodedStabIncompHighReFluid", \
	                    "3D4nodedStabCompFluid",         \
	                    "3D4nodedLinearSolid",           \
                            "3D10nodedQuadraticSolid",       \
	                    "3D8nodedLinearSolid",           \
	                    "3D8nodedFBarSolid",             \
                            "3D20nodedQuadraticSolid",       \
                            "2D2nodedLine",                  \
                            "2D3nodedLine",                  \
	                    "2D2nodedGeomExSmallStrainBeam", \
	                    "2D2nodedKirchhoffBeam",         \
	                    "2D2nodedTruss",                 \
	                    "1D2nodedLine",                  \
	                    "1D2nodedAdvectionDiffusion",    \
	                    "1D2nodedSchroedingerCN",        \
	                    "1D2nodedSchroedingerST",        \
                            "2D2nodedPressureLoad",          \
                            "2D3nodedPressureLoad",          \
                            "1D2nodedPipeFlowST",            \
                            "1D2nodedFlexibleWing",          \
                            "2D3nodedLinearPoisson",         \
                            "1D2nodedLinearSolidALE",        \
                            "2D2nodedFreeSurface",   NULL}



#define ISOGEOM_ELEMENT_TYPE_NAMES {"NurbsElem1DAdvectionDiffusion",  \
                                    "NurbsElem1DElasticBar",     \
                                    "NurbsElem1DEulerBeam",     \
                                    "NurbsElem1DElasticBarLSFEM",     \
                                    "NurbsElem2DStructSolid", \
                                    "NurbsElem2DStructFbarSolid", \
                                    "NurbsElem2DStructBbarSolid", \
                                    "NurbsElem2DStructMixed2field", \
                                    "NurbsElem2DStructMixed3field", \
                                    "NurbsElemKirchhoffPlate", \
                                    "NurbsElemMindlinPlate", \
                                    "NurbsElem3DStructSolid", \
                                    "NurbsElem3DStructMixed2field", \
                                    "NurbsElem2DAdvectionDiffusion", \
                                    "NurbsElem2DStokes", \
                                    "NurbsElem2DNavierStokes3dof", \
                                    "NurbsElem2DNavierStokes4dof", \
                                    "NurbsElem2DStructSolidLSFEM2dof",\
                                    "NurbsElem2DStructSolidLSFEM3dof",\
                                    "NurbsElem2DHeatTransfer", \
                                    "NurbsElem2DTempCoupled4dof", \
                                    "NurbsElem2DStructMixed2fieldStabilised", \
                                    "NurbsElem3DStructMixed2fieldStabilised", NULL}




#define MATERIAL_TYPE_NAMES {"smallStrainElasticity",          \
	                           "NeoHookeElasticity",             \
     	                       "OgdenElasticity",                \
	                           "vonMisesElastoPlasticity",       \
	                           "smallStrainVonMisesEP",          \
	                           "vonMisesEPWithKinematicHardening", \
			                       "multiscaleMaterial2D",           \
			                       "multiscaleMaterial3D",           \
	                           "dummy",                          \
	                           "NeoHookeElasticity1D",           \
	                           "sfVonMisesIsotropicElpl",        \
   	 	                       "sHillAnisotropicElpl",           \
                             "testTangent",                    \
                             "EVEP",                           \
                             "HenckyElasticity",               \
                             "StVenantKirchhoffElasticity",    \
                             "genericHyperelasticity",         \
                             "smallStrainAnisotropicElasticity",\
                             "smallStrainElasticLamina2D",     \
 	                           "HyplasELASTC",                   \
	                           "HyplasTRESCA",                   \
	                           "HyplasVMISES",                   \
	                           "HyplasMOHCOU",                   \
	                           "HyplasDRUPRA",                   \
	                           "HyplasCAPDP",                    \
	                           "HyplasLEMDAM",                   \
	                           "HyplasDAMELA",                   \
	                           "HyplasOGDEN",                    \
	                           "HyplasPDSCRY",                   \
	                           "HyplasVMMIXD", NULL }


#define MATERIAL_DIMENSIONS { 3,  \
	                      3,  \
			      3,  \
			      3,  \
	                      3,  \
			      3,  \
			      2,  \
			      3,  \
			      3,  \
			      1,  \
                              3,  \
                              3,  \
                              3,  \
                              3,  \
                              3,  \
                              3,  \
                              3,  \
                              3,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2,  \
                              2   }


#define ALE_TYPE_NAMES {"fixed",                     \
	                "cellCentroid",              \
                        "aspectRatio",               \
	                "linearPseudoElasticity",    \
	                "nonlinearPseudoElasticity", \
                        "dummy",                   NULL}


#define WRND_TYPE_NAMES { "u", "du", "ddu", "reac", "x", "x0", "int(uu+vv)", "outp", NULL }


#define COMMAND_LINE_ARGUMENTS { "-noGUI",           \
                                 "-batch",           \
                                 "-debug",           \
                                 "-wulf",            \
                                 "-sony",            \
                                 "-deniz",           \
	                         "-keep",            \
                                 "-aspectRatioCorr", \
	                         "-readIFileInfo",   \
	                         "-lastProj",        \
                                 "-test",   NULL }


#define LAST_PROJECT "mpap2.name"


#define RUN_CONTROL_COMMANDS { "BATCH", "INTER", NULL } 



#define STATUS_ENUM { UNDEF,            \
	              NOPROJECT,        \
	              LOADING,          \
	              INTERACTIVE,      \
	              INTERACTIVEBUSY,  \
	              BATCHMODE,        \
                      SELECTZOOMBOX,    \
                      SELECTSEARCHBOX,  \
                      EXECTESTMACRO,    \
                      PRESSMOUSE,       \
                      CHANGEPERSPECTIVE,\
                      SELECTNODE        }  

#define STATUS_NAMES { "status undefined",                \
	               "no project loaded",               \
	               "loading ...",                     \
	               "interactive",                     \
	               "interactive, but busy ...",       \
	               "batch mode",                      \
	               "select zoom box",                 \
                       "select search box",               \
                       "test macro",                      \
	               "press mouse button",              \
                       "change perspective",              \
                       "select node",  NULL }

#define RUN_MODE_PROMPTS { "mpap", "batch", "pre", "noproj", NULL}

#define RUN_MODE_ENUM { INTER, BATCH, PRE, NOPROJ }



#define MAX_MACRO 230

#define TD_DIM    100

#define MAX_PROCESSORS 64


#ifndef BIT64

  #define VOID_PTR int

#else

  #define VOID_PTR long

#endif



#define COLOURS_ENUM { RED, BLUE, GREEN, YELLOW, CYAN, MAGENTA, LIGHTBLUE, WHITE, BLACK }

#define COLOUR_NAMES \
	"red","blue","green","yellow","cyan","magenta","lightblue","white","black" 



#ifndef WINDOWS 

  #define SLASH '/'
  #define DELETE_FILE_COMMAND "rm"

  #define COLOURS           \
	"red","blue","green","yellow","cyan","magenta","lightblue","white","black" 
  #define COLOURS_RED       \
	"*red","blue","green","yellow","cyan","magenta","lightblue","white","black"
  #define COLOURS_BLUE      \
	"red","*blue","green","yellow","cyan","magenta","lightblue","white","black"
  #define COLOURS_GREEN     \
	"red","blue","*green","yellow","cyan","magenta","lightblue","white","black"
  #define COLOURS_YELLOW    \
	"red","blue","green","*yellow","cyan","magenta","lightblue","white","black"
  #define COLOURS_CYAN      \
	"red","blue","green","yellow","*cyan","magenta","lightblue","white","black"
  #define COLOURS_MAGENTA   \
	"red","blue","green","yellow","cyan","*magenta","lightblue","white","black"
  #define COLOURS_LIGHTBLUE \
        "red","blue","green","yellow","cyan","magenta","*lightblue","white","black"
  #define COLOURS_WHITE     \
	"red","blue","green","yellow","cyan","magenta","lightblue","*white","black"

#else

  #define SLASH '\\'
  #define DELETE_FILE_COMMAND "del"

  #define COLOURS           \
	"RED","BLUE","GREEN","YELLOW","CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_RED       \
	"*RED","BLUE","GREEN","YELLOW","CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_BLUE      \
	"RED","*BLUE","GREEN","YELLOW","CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_GREEN     \
	"RED","BLUE","*GREEN","YELLOW","CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_YELLOW    \
	"RED","BLUE","GREEN","*YELLOW","CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_CYAN      \
	"RED","BLUE","GREEN","YELLOW","*CYAN","MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_MAGENTA   \
	"RED","BLUE","GREEN","YELLOW","CYAN","*MAGENTA","LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_LIGHTBLUE \
        "RED","BLUE","GREEN","YELLOW","CYAN","MAGENTA","*LIGHT BLUE","WHITE","BLACK"
  #define COLOURS_WHITE     \
	"RED","BLUE","GREEN","YELLOW","CYAN","MAGENTA","LIGHT BLUE","*WHITE","BLACK"
#endif


#define NOO     {"n", "N", "no", "No", NULL}
#define YESS    {"y", "Y", "yes", "Yes", NULL}
#define QUIT   {"q", "Q", "x", "X", "quit", "Quit", "exit", "Exit", NULL}
#define CANCEL {"c", "C", "cancel", "Cancel", NULL}

#define COUT std::cout << "          "

#define hello1 std::cout << "hello 1\n";
#define hello2 std::cout << "hello 2\n";
#define hello3 std::cout << "hello 3\n";
#define hello4 std::cout << "hello 4\n";
#define hello5 std::cout << "hello 5\n";
#define hello6 std::cout << "hello 6\n";
#define hello7 std::cout << "hello 7\n";
#define hello8 std::cout << "hello 8\n";
#define hello9 std::cout << "hello 9\n";
#define hello10 std::cout << "hello 10\n";


#endif

