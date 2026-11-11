#ifndef SCI_decl_hpp
#define SCI_decl_hpp
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/problems/abstract/TimeProblem.hpp"
#include "feddlib/problems/specific/DiffusionReaction.hpp"
#include "feddlib/core/General/ExporterTxt.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/Solver/TimeSteppingTools.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"


namespace FEDD{

/*!
    Solid–chemistry/pharmaco-mechanical interaction problem.
    Displacement and concentration form a coupled 2 x 2 system, or chemistry
    can be solved separately with the supplied explicit coupling option.

*/
template <class SC , class LO , class GO , class NO >
class TimeProblem;
template <class SC , class LO , class GO , class NO >
class DiffusionReaction;
template <class SC , class LO , class GO , class NO >
class NonLinElasticity;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class SCI : public NonLinearProblem<SC,LO,GO,NO>  {

public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::Domain_Type Domain_Type;
    typedef Teuchos::RCP<Domain_Type > DomainPtr_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;

    typedef typename Problem_Type::Domain_Type::Mesh_Type Mesh_Type;
    typedef typename Problem_Type::Domain_Type::MeshPtr_Type MeshPtr_Type;

    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;

    typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type BlockMultiVectorPtrArray_Type;

    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef Teuchos::RCP<TimeProblem_Type> TimeProblemPtr_Type;

    typedef DiffusionReaction<SC,LO,GO,NO> ChemProblem_Type;
    typedef NonLinElasticity<SC,LO,GO,NO> StructureNonLinProblem_Type;

    typedef Teuchos::RCP<ChemProblem_Type> ChemProblemPtr_Type;
    typedef Teuchos::RCP<StructureNonLinProblem_Type> StructureNonLinProblemPtr_Type;

    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename Problem_Type::BC_Type BC_Type;
    typedef typename Teuchos::RCP<BC_Type> BCPtr_Type;

    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;

    typedef ExporterParaView<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;

    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;

    // Coupled AceGen assembly requires matching 3D P2 solid and chemistry meshes.
    SCI(const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
                    const DomainConstPtr_Type &domainChem, std::string FETypeChem,vec2D_dbl_Type diffusionTensor, RhsFunc_Type reactionFunc,
                    ParameterListPtr_Type parameterListStructure, ParameterListPtr_Type parameterListChem,
                    ParameterListPtr_Type parameterListSCI, Teuchos::RCP<SmallMatrix<int> > &defTS);

    ~SCI();

    virtual void info();

    virtual void assemble( std::string type = "" ) const;

    // Share the coupled vectors with the individual solid and chemistry problems.
    void setFromPartialVectorsInit() const;

    // Setze die aktuelle Loesung als vergangene Loesung
    void updateMeshDisplacement() const;

    // Berechnet die Massematrix und die daraus resultierende rechte Seite nach BDF2
    // void getFluidMassmatrixAndRHSInTime(BlockMatrixPtr_Type massmatrix, BlockMultiVectorPtr_Type rhs) const;

    // Berechne die Massematrix fuer das FluidProblem.
    void setChemMassmatrix(MatrixPtr_Type& massmatrix) const;

    // Berechnet die Massematrix und die daraus resultierende rechte Seite nach Newmark
    // und macht direkt ein Update. Dies koennen wir bei Struktur machen, da Massematrix
    // innerhalb einer Zeitschleife konstant ist.
    void setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const;

    // Advance Newmark history only at the start of a step; add the follower-load
    // tangent only after assembling a fresh material Jacobian.
    void computeSolidRHSInTime(bool updateHistory = true, bool updateJacobian = false) const;

    void computeChemRHSInTime() const;

    // Hier wird timeSteppingTool_->t_ inkrementiert
    void updateTime() const;

    // Hier wird im Prinzip updateSolution() fuer problemTimeChem_ aufgerufen
    void updateChemInTime() const;

    // Solving chemistry component of problem in case of explicit chemistry
    void solveChemistryProblem() const;

    // Verschiebt die notwendigen Gitter
    void moveMesh() const;

    void initializeCE();
    // Set up the solid and chemistry time problems.
    void setupSubTimeProblems(ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListStructure) const;

    void setBoundariesSubProblems() const;
    ChemProblemPtr_Type getChemProblem(){
        return problemChem_;
    }

    StructureNonLinProblemPtr_Type getStructureProblem(){
        return problemStructureNonLin_;
    }

    // Berechnet von einer dofID, d.h. dim*nodeID+(0,1,2), die entsprechende nodeID.
    // IN localDofNumber steht dann, ob es die x- (=0), y- (=1) oder z-Komponente (=2) ist.
    void toNodeID(UN dim, GO dofID, GO& nodeID, LO& localDofNumber ) const
    {
        nodeID = (GO) (dofID/dim);
        localDofNumber = (LO) (dofID%dim);
    }

    // Diese Funktion berechnet genau das umgekehrte. Also von einer nodeID die entsprechende dofID
    void toDofID(UN dim, GO nodeID, LO localDofNumber, GO& dofID)  const
    {
        dofID = (GO) ( dim * nodeID + localDofNumber);
    }

    // void getValuesOfInterest2DBenchmark( vec_dbl_Type& values );

    // void getValuesOfInterest3DBenchmark( vec_dbl_Type& values );

    virtual void getValuesOfInterest( vec_dbl_Type& values ) {}  ;

    virtual void getValuesOfInterest( BlockMultiVectorPtr_Type& values );

    virtual void exportValuesOfInterest(double time) override;
    virtual void importValuesOfInterest(double time) override;

    virtual void computeValuesOfInterestAndExport() {} ;

    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const{};

    virtual void reAssemble(std::string type) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions) {};

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse

    // Element output contains only the first four fields (stresses).
    BlockMultiVectorPtr_Type getPostProcessingData(bool elementAverage = false) const;

    vec_string_Type getPostprocessingNames() const;

    /*####################*/

    MultiVectorPtr_Type meshDisplacementOld_rep_;
    MultiVectorPtr_Type meshDisplacementNew_rep_;
    MultiVectorPtr_Type c_rep_;
    MultiVectorPtr_Type d_rep_;

    // stationaere Systeme
    ChemProblemPtr_Type problemChem_;
    StructureNonLinProblemPtr_Type problemStructureNonLin_;

    // zeitabhaengige Systeme
    mutable TimeProblemPtr_Type problemTimeChem_;
    mutable TimeProblemPtr_Type problemTimeStructure_;

    Teuchos::RCP<SmallMatrix<int>> defTS_;
    mutable Teuchos::RCP<TimeSteppingTools> timeSteppingTool_;

private:
    std::string materialModel_;
    vec_dbl_Type valuesForExport_;
    vec_string_Type postProcessingnames_;
    bool geometryExplicit_;
    mutable BlockMatrixPtr_Type systemC_;
    ExporterTxtPtr_Type exporterIterationsChem_;
    mutable ExporterPtr_Type exporterEMod_;
    mutable ExporterPtr_Type exporterChem_;
    mutable bool exportedEMod_ ;
    mutable bool setUpTimeStep_;
    mutable MultiVectorPtr_Type eModVec_;
    bool loadStepping_;
    bool chemistryExplicit_;
    bool externalForce_;
    bool nonlinearExternalForce_;
    /*####################*/

public:


private:

};
}
#endif
