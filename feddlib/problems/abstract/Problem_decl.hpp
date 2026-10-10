#ifndef PROBLEM_DECL_hpp
#define PROBLEM_DECL_hpp

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_PreconditionerBase.hpp>

#include "git_version.h"

#ifdef FEDD_HAVE_TEKO
#include <Teko_StratimikosFactory.hpp>
#endif

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/Checkpointing/RecoveryCheckpoint.hpp"

/*!
 Declaration of Problem

 @brief  Problem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template<class SC_, class LO_, class GO_, class NO_>
class BlockMultiVector;
template<class SC_, class LO_, class GO_, class NO_>
class BCBuilder;
template<class SC_, class LO_, class GO_, class NO_>
class BlockMatrix;
template<class SC_, class LO_, class GO_, class NO_>
class Domain;
template<class SC_, class LO_, class GO_, class NO_>
class FE;
template<class SC_, class LO_, class GO_, class NO_>
class Preconditioner;

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Problem {

public:

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<Domain_Type> DomainPtr_Type;
    typedef Teuchos::RCP<const Domain_Type> DomainConstPtr_Type;
    typedef std::vector<DomainConstPtr_Type> DomainConstPtr_vec_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef typename Matrix_Type::Map_Type Map_Type;
    typedef typename Matrix_Type::MapPtr_Type MapPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
    typedef Teuchos::RCP<const BlockMatrix_Type> BlockMatrixConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef Preconditioner<SC,LO,GO,NO> Preconditioner_Type;
    typedef Teuchos::RCP<Preconditioner_Type> PreconditionerPtr_Type;
    typedef Teuchos::RCP<const Preconditioner_Type> PreconditionerConstPtr_Type;

    typedef BCBuilder<SC,LO,GO,NO> BC_Type;
    typedef Teuchos::RCP<BC_Type> BCPtr_Type;
    typedef Teuchos::RCP<const BC_Type> BCConstPtr_Type;

    typedef FE<SC,LO,GO,NO> FEFac_Type;
    typedef Teuchos::RCP<FEFac_Type> FEFacPtr_Type;
    typedef Teuchos::RCP<const FEFac_Type> FEFacConstPtr_Type;

    typedef Teuchos::ParameterList ParameterList_Type;
    typedef Teuchos::RCP<ParameterList_Type> ParameterListPtr_Type;

    typedef Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder >  LinSolverBuilderPtr_Type;
    typedef Teuchos::RCP<const Stratimikos::DefaultLinearSolverBuilder >  LinSolverBuilderConstPtr_Type;

    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    typedef Teuchos::ArrayRCP<GO> GOVecPtr;

    typedef std::vector<std::string> string_vec_Type;

    typedef Teuchos::RCP<Thyra::PreconditionerBase<SC> > ThyraPrecPtr_Type;
    using ThyraTypes = ThyraTypedefs<SC>;
    using ThyraLinOpPtr_Type = Teuchos::RCP<typename ThyraTypes::ThyraOp_Type>;

    Problem(CommConstPtr_Type comm);

    Problem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    virtual ~Problem() = default;

    virtual void info() = 0;

    /*! If we have a Nonlinear Problem, we need the opportunity to know about the current Newton Step -> override in NonLinearProblem
        Default implementation just returns 0
    */
    virtual int getNonlinearIterationStep() const{return 0; } 

    void infoProblem();

    void infoParameter(bool full = true, std::string ="empty");

    void addVariable(const DomainConstPtr_Type &domain, std::string FEType, std::string name, int dofsPerNode);

    /*! Add right hand side function for each block, if you want to skip a block add a dummy function
        -> Error Warning. In case initializeProblem() is called beforehand the rhsVec is already initialized with a certain size.
        --> This leads to a problem, if addRhsFunction is called after, as it is a push_back operation.
    
    */
    void addRhsFunction(RhsFunc_Type func);

    /*! Add right hand side function for block i   */
    void addRhsFunction(RhsFunc_Type func,int i);

    /*! Adds rhs function and specifies it so certain flag */
   // void addRhsFunctionAndFlag(RhsFunc_Type func, int i, int flag);    

    RhsFunc_Type& getRhsFunction( int i );

    virtual void assemble( std::string type ) const = 0;
    
//    virtual void assembleExternal( std::string type ) = 0;

    //void reAssemble();
    
    void assembleSourceTerm( double time = 0. ) const;
    
    void assembleVolumeTerm( double time ) const;
    
    void assembleSurfaceTerm( double time ) const;
    
    bool hasSourceTerm() const;

    int solve( BlockMultiVectorPtr_Type rhs = Teuchos::null );

    void setupPreconditioner( std::string type="Monolithic" ) const;

    void initializePreconditioner( std::string type="Monolithic" ) const;

    void addBoundaries(const BCConstPtr_Type &bcFactory);

    void setBoundaries(double time=.0) const;

    void setBoundariesRHS(double time=.0) const;
    
    void setAllDirichletZero( BlockMultiVectorPtr_Type rhs) const;
    
    void setBoundariesSystem() const;

    void initializeProblem(int nmbVectors=1);
    
    /// Allocate solution, RHS and source vectors without reading checkpoint files.
    void initializeVectors(int nmbVectors=1);

    /**
     * @brief Restore primary solution fields from a checkpoint.
     *
     * Reads each field from `Solution<variable>.h5` in the configured restart
     * directory, using std::to_string(restartTime) as the checkpoint key and
     * the allocated solution block's map for the distributed import.
     * The ALE geometry field d_f is skipped; FSI restores it separately.
     * Time integration history, the simulation clock and operators are not
     * restored here.
     *
     * @param[in] restartTime Physical time used to identify the checkpoint.
     * @pre initializeVectors() has allocated all solution blocks and maps.
     * @pre Call before assembly or linking solution blocks to subproblems,
     *      since importing replaces the solution block pointers.
     * @note Call collectively on the problem communicator. The caller decides
     *       whether restart is enabled; this function always reads the files.
     * @throws std::logic_error If the solution vector has not been allocated.
     * @see initializeProblem()
     */
    void restoreSolutionFromCheckpoint(double restartTime);

    /** @brief Load Navier–Stokes velocity and pressure as a new initial condition.
     * Validates mesh/discretization metadata and stages all fields before changing
     * solution blocks. Reads "Initial solution directory" at sourceTime; does not
     * restore history, alter the simulation clock, assemble, or solve a timestep.
     * @param sourceTime Time key in the source checkpoint, independent of the new dt.
     * @pre initializeVectors() has allocated solution blocks; call before assembly.
     * @note Collective operation. Requires standalone multistep Navier–Stokes.
     */
    void initializeSolutionFromCheckpoint(double sourceTime);

    /** @brief Cache the reference mesh/discretization and validate restart data.
     * Collects field descriptions; core checkpoint helpers construct the schema
     * and compute mesh identity independently of the Problem class.
     * Call collectively before the first restore or initial solution load,
     * including FSI constructor reads. Initial solution mode requires the field
     * schema but validates primary fields separately from restart history.
     * A component role selects FSI fluid mass history or structure Newmark history.
     * Cached mesh identity is retained when ALE subsequently moves the coordinates.
     */
    void prepareCheckpointMetadata(const std::string& role = "");

    /// Validate an explicitly selected checkpoint, once per time, before loading values.
    void validateRestartCheckpoint(double restartTime);

    /// Write the versioned manifest alongside an already requested checkpoint.
    void writeCheckpointMetadata(double time);


    BlockMultiVectorPtr_Type getRhs();

    BlockMultiVectorPtr_Type getRhs() const;

    BlockMultiVectorPtr_Type getSolution();

    BlockMatrixPtr_Type getSystem() const;

    PreconditionerPtr_Type getPreconditioner();

    PreconditionerConstPtr_Type getPreconditionerConst() const;
 
    void setPreconditionerThyraFromLinOp( ThyraLinOpPtr_Type precLinOp );
    
    void initializeSolverBuilder() const;

    bool getVerbose() const;

    FEFacConstPtr_Type getFEFactory();

    BCConstPtr_Type getBCFactory();

    DomainConstPtr_Type getDomain(int i) const;

    DomainConstPtr_vec_Type getDomainVector() const{
        return domainPtr_vec_;
    }
    
    std::string getFEType(int i) const;

    std::string getVariableName(int i) const;

    int getDofsPerNode(int i) const;

    ParameterListPtr_Type getParameterList() const;

    void addToRhs(BlockMultiVectorPtr_Type x) const;
    
    BlockMultiVectorPtr_Type getSourceTerm();

    void initSolutionWithVector(MultiVector_Type& mv);

    LinSolverBuilderPtr_Type getLinearSolverBuilder() const{return linearSolverBuilder_;}

    CommConstPtr_Type getComm() const{return comm_;}

    virtual void getValuesOfInterest( vec_dbl_Type& values ) = 0 ;

    virtual void computeValuesOfInterestAndExport() = 0;
    
    /// Write additional problem-specific state at a checkpoint. Nothing by default.
    virtual void exportValuesOfInterest(double time) {}

    /// Restore additional problem-specific state at the restart time. Nothing by default.
    virtual void importValuesOfInterest(double time) {}

    using RecoverySnapshot = checkpoint::RecoveryCheckpointSnapshot<SC,LO,GO,NO>;
    using RecoverySnapshotPtr = Teuchos::RCP<RecoverySnapshot>;
    /// Attach the current independent recovery capture; FSI also attaches its components.
    virtual void setRecoverySnapshot(RecoverySnapshotPtr snapshot) { recoverySnapshot_ = snapshot; }
    RecoverySnapshotPtr getRecoverySnapshot() const { return recoverySnapshot_; }
    /// Capture model-specific state without exporting or modifying it.
    virtual void captureAdditionalRecoveryState(double time) {}
    /// Add this problem's compatibility manifest to the attached recovery capture.
    void captureRecoveryMetadata() {
        prepareCheckpointMetadata();
        if (!recoverySnapshot_.is_null()) recoverySnapshot_->addManifest(checkpointSchema_,
            checkpoint::clockState(parameterList_, recoverySnapshot_->time()));
    }

    void addParemeterRhs(double para){ parasSourceFunc_.push_back( para ); }
    
    void changeAssFELinearization(std::string linearization); // Function in order to be able to change e.g. from FixedPoint to Newton linearization on element level

	double calculateH1Norm(MultiVectorConstPtr_Type mv, int blockId1=0, int blockId2=0, int domainInd=0); // Function that calculates H1 Error in the 'mv * K * mv' sense, with K beeing the Stiffness Matrix

	double calculateL2Norm(MultiVectorConstPtr_Type mv, int domainInd=0); // Function that calculates L2 Error in the 'mv * M * mv' sense, with M beeing the Mass Matrix


    int dim_;
    mutable CommConstPtr_Type comm_;
    mutable BlockMatrixPtr_Type system_;
    mutable BlockMultiVectorPtr_Type rhs_;
    mutable BlockMultiVectorPtr_Type solution_;
    PreconditionerPtr_Type preconditioner_;
    LinSolverBuilderPtr_Type linearSolverBuilder_;

    bool verbose_;

    std::vector<RhsFunc_Type>   rhsFuncVec_; // RHS functions of different blocks
    vec_dbl_Type parasSourceFunc_; //
    
protected:

    /// Stage and load primary fields after the caller has validated their compatibility.
    void loadSolutionFields(const std::string& directory, double sourceTime);

    Teuchos::ParameterList checkpointSchema_;
    bool checkpointSchemaPrepared_ = false;
    bool checkpointValidated_ = false;
    double checkpointValidatedTime_ = 0.;

    mutable ParameterListPtr_Type	parameterList_;
    RecoverySnapshotPtr recoverySnapshot_;
    mutable DomainConstPtr_vec_Type domainPtr_vec_;
    string_vec_Type                 domain_FEType_vec_;
    string_vec_Type                 variableName_vec_;
    mutable BCConstPtr_Type         bcFactory_;

    FEFacPtr_Type feFactory_;
    std::vector<int> dofsPerNode_vec_;
    
    /*!  sourceTerm_: Is a source term or a surface integral. Fill parasSourceFunc_ for additional parameters */
    BlockMultiVectorPtr_Type    sourceTerm_; // BlockMV of all assembled RHS functions
    
#ifdef FEDD_TIMER
    TimePtr_Type solveProblemTimer_;
    TimePtr_Type bcMatrixTimer_;
    TimePtr_Type bcRHSTimer_;
#endif



};
}
#endif
