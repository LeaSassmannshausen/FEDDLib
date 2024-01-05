#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/BCBuilder.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 MeshUnstructured test

 @brief  MeshUnstructured test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */



using namespace std;
using namespace Teuchos;

typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
using namespace FEDD;
int main(int argc, char *argv[]) {

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra"; //this does nothing atm
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string filename = "fluidBenchmark2.mesh";
    myCLP.setOption("file",&filename,"Mesh filename");
    string exportfilename = "export.mesh";
    myCLP.setOption("fileExport",&exportfilename,"Export Mesh filename");
    int dim = 2;
    myCLP.setOption("dim",&dim,"Dimension");
    string delimiter = " ";
    myCLP.setOption("delimiter",&delimiter,"Delimiter in mesh-file");
    string FEType="P1";
    myCLP.setOption("FEType",&FEType,"FEType");
   /* bool exportEdges= true;
    myCLP.setOption("exportEdges",&exportEdges,"Exporting Edges");
    bool exportSurfaces= true;
    myCLP.setOption("exportSurfaces",&exportSurfaces,"Exporting Surfaces");
    */
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    // Mesh
    
    int numProcsCoarseSolve = 0;
    bool boolExportMesh = true;
    bool boolExportSubdomains = false;
    int volumeID = 10;
    if (filename=="some_tetrahedron.mesh")
        volumeID = 12;

    DomainPtr_Type domainP1;
    DomainPtr_Type domainP2;
    DomainPtr_Type domain;
    
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    
    domainP1.reset( new Domain_Type( comm, dim ) );
    domainP2.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domainP1;

    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
    
    partitionerP1.readAndPartition(volumeID);

    if (FEType == "P2") {
        domainP2->buildP2ofP1Domain( domainP1 );
        domain = domainP2;
    }
    else
        domain = domainP1;

    domain->exportMesh(true,true,exportfilename);
    
 	FE<SC,LO,GO,NO> fe_1;
    fe_1.addFE(domain);
    fe_1.checkMeshOrientation(dim,FEType);
    
    comm->barrier();
    
    if(FEType =="P1"){
		ParameterListPtr_Type pListPartitionerTest = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
		pListPartitionerTest->set( "Mesh 1 Name", exportfilename );
		
		DomainPtr_Type domainP1Test;
		DomainPtr_Type domainTest;

		domainP1Test.reset( new Domain_Type( comm, dim ) );
		MeshPartitioner_Type::DomainPtrArray_Type domainP1ArrayTest(1);
		domainP1ArrayTest[0] = domainP1Test;

		MeshPartitioner<SC,LO,GO,NO> partitionerP1Test ( domainP1ArrayTest, pListPartitionerTest, "P1", dim );
		comm->barrier();

		partitionerP1Test.readAndPartition(volumeID);

		domainTest = domainP1Test;
	   
		FE<SC,LO,GO,NO> fe_2;
		fe_2.addFE(domainTest);
		fe_2.checkMeshOrientation(dim,FEType);
   	}
    

    


    return(EXIT_SUCCESS);
}
