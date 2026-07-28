#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/FE/Domain.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Mesh Partitioner test

 @brief  MeshUnstructured test
 @author Lea Saßmannshausen
 @version 1.0
 @copyright LS
 */



using namespace std;
using namespace Teuchos;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
int main(int argc, char *argv[]) {

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string filename = "plaque_solid_length_0_5.mesh";
    int dim = 3;
    string delimiter = " ";
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    // Mesh
    
    int numProcsCoarseSolve = 0;
  
    DomainPtr_Type domain1; // unweighted graph
    DomainPtr_Type domain2; // weighted graph

    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    pListPartitioner->set( "Weight ID", 21 );

    // Unweighted graph
    domain1.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainArray1(1);
    domainArray1[0] = domain1;

    MeshPartitioner<SC,LO,GO,NO> partitioner1 ( domainArray1, pListPartitioner, "P1", dim );
    
    partitioner1.readAndPartition();

    domain1->exportElementFlags("metis_elements");
    domain1->exportDistribution("metis_unweighted_graph_partition");

    // Weighted graph
    domain2.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainArray2(1);
    domainArray2[0] = domain2;

    MeshPartitioner<SC,LO,GO,NO> partitioner2 ( domainArray2, pListPartitioner, "P1", dim );
    partitioner2.readAndPartitionCustom();

    domain2->exportDistribution("metis_weighted_graph_partition");

    return(EXIT_SUCCESS);
}
