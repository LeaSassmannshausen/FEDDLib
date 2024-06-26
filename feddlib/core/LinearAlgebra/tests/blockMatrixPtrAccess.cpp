#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Matrix test

 @brief  Matrix test
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
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    GO numGlobalElements = 4;
    myCLP.setOption("nge",&numGlobalElements,"numGlobalElements.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

     typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef RCP<TpetraMap_Type> TpetraMapPtr_Type;
    typedef RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef RCP<Matrix_Type> MatrixPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef Map_Tpetra<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;

    TpetraMapConstPtr_Type tmap = RCP(new TpetraMap_Type(numGlobalElements, 0, commWorld));
    TpetraMatrixPtr_Type tmatrix = RCP( new TpetraMatrix_Type(tmap, 1));
   
    MatrixPtr_Type matrix = rcp( new Matrix_Type( tmatrix ) );

    for (int i=0; i<matrix->getMap()->getNodeNumElements(); i++) {
        Array<SC> values( 1 , 1.);
        Array<GO> indicesCol( 1 , matrix->getMap()->getGlobalElement( i ) );
        matrix->insertGlobalValues( indicesCol[0], indicesCol(), values() );
    }
    matrix->fillComplete();
    matrix->print();

    BlockMatrixPtr_Type bmatrix = Teuchos::rcp(new BlockMatrix_Type( 1 ) );
    bmatrix->addBlock(matrix,0,0);
    
    bmatrix->print();
    
    matrix->resumeFill();
    matrix->scale(2.);
    matrix->fillComplete();
    
    matrix->print();
    bmatrix->print();
    
    return(EXIT_SUCCESS);
}
