#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 blockView test

 @brief  blockView
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

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    GO numGlobalElements = 2;
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

    TpetraMapConstPtr_Type tmap = RCP(new TpetraMap_Type(numGlobalElements, 0, comm));
    TpetraMatrixPtr_Type tmatrix = RCP( new TpetraMatrix_Type(tmap, 1));
    TpetraMatrixPtr_Type tmatrix2 = RCP( new TpetraMatrix_Type(tmap, 1));


    MatrixPtr_Type matrix1 = rcp( new Matrix_Type( tmatrix ) );
    MatrixPtr_Type matrix2 = rcp( new Matrix_Type( tmatrix2 ) );
    MapConstPtr_Type map = matrix1->getMap();
    for (UN i=0 ; i<matrix1->getNodeNumRows(); i++) {
        Array<GO> indices(1);
        Array<SC> values1(1,1.);
        Array<SC> values2(1,-1.);
        indices[0] = map->getGlobalElement( i );
        matrix1->insertGlobalValues( indices[0], indices(), values1() );
        matrix2->insertGlobalValues( indices[0], indices(), values2() );
    }

    matrix1->fillComplete();
    matrix2->fillComplete();

    BlockMatrixPtr_Type system = rcp(new BlockMatrix_Type(2));

    system->addBlock( matrix1, 0, 0 );
    system->addBlock( matrix2, 1, 1 );

    system->print(VERB_EXTREME);

    system->merge();

    system->printMerge(VERB_EXTREME);

    matrix1->resumeFill();
    MapConstPtr_Type colMap = matrix1->getMap("col");
    for (UN i=0 ; i<matrix1->getNodeNumRows(); i++) {
        Array<LO> indices(1);
        Array<SC> values1(1,2.);
        indices[0] = colMap->getLocalElement( map->getGlobalElement( i ) );
        matrix1->insertLocalValues( i, indices(), values1() );
    }

    matrix1->fillComplete();
    system->print(VERB_EXTREME);
    system->printMerge(VERB_EXTREME);
    cout << " Done " << endl;

    return(EXIT_SUCCESS);
}
