#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Map_Tpetra.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Map test

 @brief  Map test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
using namespace Teuchos;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
using namespace FEDD;
int main(int argc, char *argv[]) {

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > commWorld = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    GO numGlobalElements = 10;
    myCLP.setOption("nge",&numGlobalElements,"numGlobalElements.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    
    typedef Map_Tpetra<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;


    Array<GO> indices(numGlobalElements);
    for (UN i=0; i<indices.size(); i++) {
        indices[i] = i;
    }

    MapPtr_Type map = rcp( new Map_Type(commWorld->getSize()*numGlobalElements, indices(), 0, commWorld) );

    map->print();

    return(EXIT_SUCCESS);
}
