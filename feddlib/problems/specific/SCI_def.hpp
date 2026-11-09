#ifndef SCI_def_hpp
#define SCI_def_hpp
#include "SCI_decl.hpp"
#include "feddlib/core/General/BCBuilder.hpp"


namespace FEDD {
// Funktionen fuer die rechte Seite der Struktur/ Chem/ Geometrie sind im jeweiligen Problem

template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::SCI(const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
                    const DomainConstPtr_Type &domainChem, std::string FETypeChem, vec2D_dbl_Type diffusionTensor, RhsFunc_Type reactionFunc,
                    ParameterListPtr_Type parameterListStructure, ParameterListPtr_Type parameterListChem,
                    ParameterListPtr_Type parameterListSCI, Teuchos::RCP<SmallMatrix<int> > &defTS):
NonLinearProblem<SC,LO,GO,NO>( parameterListSCI, domainChem->getComm() ),
problemChem_(),
problemStructureNonLin_(),
meshDisplacementOld_rep_(),
meshDisplacementNew_rep_(),
c_rep_(),
defTS_(defTS),
timeSteppingTool_(),
exporterEMod_(),
materialModel_( parameterListSCI->sublist("Parameter").get("Structure Model","SCI_NH") )
{
    //this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);

    this->initNOXParameters();

    timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));

    //std::string linearization = parameterListSCI->sublist("General").get("Linearization","FixedPoint");

    //TEUCHOS_TEST_FOR_EXCEPTION( !(linearization == "Newton"|| linearization == "NOX")  && materialModel_ != "linear", std::runtime_error, "Nonlinear material models can only be used with Newton's method or FixedPoint (nonlinear material Jacobian will still be used).");
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); // Structure
    this->addVariable( domainChem, FETypeChem, "c", 1); // Chemistry scalar valued problem

    this->dim_ = this->getDomain(0)->getDimension();

    TEUCHOS_TEST_FOR_EXCEPTION(materialModel_ != "SCI_NH" && materialModel_ != "SCI_SMC_MLCK" &&
        materialModel_ != "SCI_SMC_Active_Growth_Reorientation", std::logic_error,
        "SCI requires a nonlinear SCI material; the linear option has been removed.");
    auto& stepping = this->parameterList_->sublist("Timestepping Parameter");
    TEUCHOS_TEST_FOR_EXCEPTION(stepping.get("Restart", false) || stepping.get("Checkpointing", false) ||
        stepping.get("Import history", false) || stepping.get("Export history", false), std::logic_error,
        "SCI material-history file I/O and restart are not ported; disable Restart, Checkpointing, Import history and Export history.");
    problemStructureNonLin_ = Teuchos::rcp(new StructureNonLinProblem_Type(domainStructure, FETypeStructure, parameterListStructure));
    problemStructureNonLin_->initializeProblem();
    problemChem_ = Teuchos::rcp( new ChemProblem_Type( domainChem, FETypeChem, parameterListChem, diffusionTensor, reactionFunc ) );
    problemChem_->initializeProblem();

    // assemble() shares the coupled solution vectors with the subproblems.
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

    d_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    c_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );

    exportedEMod_ = false;
    setUpTimeStep_=false;
    eModVec_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getElementMap() ) );
    chemistryExplicit_ =    parameterListSCI->sublist("Parameter").get("Chemistry Explicit",false);
    loadStepping_ =    parameterListSCI->sublist("Parameter").get("Load Stepping",false);
    externalForce_ =   parameterListSCI->sublist("Parameter").get("External Force",false);
    nonlinearExternalForce_ = parameterListSCI->sublist("Parameter").get("Nonlinear External Force",false);
    this->info();

    if ( chemistryExplicit_){
        exporterIterationsChem_ = Teuchos::rcp(new ExporterTxt());
        exporterIterationsChem_->setup( "linearIterations_chem", this->comm_ );
    }

    if (materialModel_ == "SCI_SMC_Active_Growth_Reorientation")
        postProcessingnames_ = {"vonMisesStress", "SCirc", "SAxial", "SRadial", "W",
            "a11", "a12", "a13", "a21", "a22", "a23", "Stretch1", "Stretch2",
            "nC1", "nC2", "nD1", "nD2"};

}


template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::~SCI()
{
    if (!exporterEMod_.is_null()) {
       exporterEMod_->closeExporter();
    }
    if (!exporterChem_.is_null()) {
       exporterChem_->closeExporter();
    }
    if ( chemistryExplicit_){
        exporterIterationsChem_->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::info()
{
    this->infoProblem();
    if (this->verbose_)
    {
        std::cout << "\t SCI:: info() ... " << std::endl;
        std::cout << "\t SCI:: loadStepping_ " << loadStepping_  << " ... " << std::endl;
        std::cout << "\t SCI:: chemistryExplicit_ " << chemistryExplicit_  << " ... " << std::endl;
        std::cout << "\t SCI:: externalForce_ " << externalForce_  << " ... " << std::endl;
        std::cout << "\t SCI:: nonlinearExternalForce_ " << nonlinearExternalForce_  << " ... " << std::endl;
    }

    //this->infoNonlinProblem();
}

/*!
    Generally in assemble() all linear equations are assembled. Suppose we have a linear elasticity and reaction-diffusion equation,
    we get a system with linear Diagonal entries.
    The Offdiagonal parts of the system Matrix can then either be treated as nonliarities or eliminated by treating them explicitly
    in the time stepping scheme at hand.

    We start with the latter approach.

    type:
    " " -> assembly of constant matrices
    " EMod " -> assembly E module with concentration of previous timestep
    " MoveMesh " -> the move mesh operation moves the mesh, then the reaction diffusion includes the diplacement

*/
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::assemble( std::string type ) const
{

    if (type == "") {
        if (this->verbose_)
        {
            std::cout << "SCI_def: assemble() ... " << std::endl;
        }

        // Maybe nothing should happen here as there are no constant matrices
        this->system_.reset(new BlockMatrix_Type(2));
        if ( chemistryExplicit_){
            this->system_.reset(new BlockMatrix_Type(1));
            this->systemC_.reset(new BlockMatrix_Type(1));
        }


        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

        if (this->residualVec_.is_null())
            this->residualVec_.reset(new BlockMultiVector_Type(2));

        MatrixPtr_Type B(new Matrix_Type(this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type C(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type D(new Matrix_Type( this->getDomain(1)->getMapUnique(),this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ));

        // For implicit the system is ordered differently with solid block in 0,0 and diffusion in 1,1

        // As the implementation is based on a 2x2 System we use a tmp system for the actual assembly routine.
        BlockMatrixPtr_Type systemTmp(new BlockMatrix_Type(2));
        systemTmp->addBlock(A,0,0);
        systemTmp->addBlock(B,0,1);
        systemTmp->addBlock(C,1,0);
        systemTmp->addBlock(D,1,1);


        this->setupSubTimeProblems(this->problemChem_->getParameterList(), this->problemStructureNonLin_->getParameterList());

        // External pressure uses SCI's dedicated surface assembly at each step.
        // Its parameters do not follow the generic surface quadrature convention.
        if (!externalForce_) {
            this->problemTimeStructure_->assembleSourceTerm( 0. );
            this->problemTimeStructure_->addToRhs( this->problemTimeStructure_->getSourceTerm() );
        }
        this->problemTimeStructure_->setBoundariesRHS();

        this->setFromPartialVectorsInit();

        MultiVectorConstPtr_Type c;
        if(chemistryExplicit_)
            c= this->problemTimeChem_->getSolution()->getBlock(0);
        else
            c = this->solution_->getBlock(1);

        c_rep_->importFromVector(c, true);

        MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
        d_rep_->importFromVector(d, true);

        this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,systemTmp,this->residualVec_, this->parameterList_, "Jacobian", true/*call fillComplete*/);

        if(chemistryExplicit_){
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->systemC_->addBlock(systemTmp->getBlock(1,1),0,0);
        }
        else{
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->system_->addBlock(systemTmp->getBlock(0,1),0,1);
            this->system_->addBlock(systemTmp->getBlock(1,0),1,0);
            this->system_->addBlock(systemTmp->getBlock(1,1),1,1);
        }

        if ( chemistryExplicit_ && this->parameterList_->sublist("General").get("ParaViewExport",false)){
            exporterChem_ = Teuchos::rcp(new Exporter_Type());

            DomainConstPtr_Type dom = this->getDomain(1);

            int exportEveryXTimesteps = this->parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string varName = "c";

            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
            exporterChem_->setup(varName, meshNonConst, dom->getFEType(), exportEveryXTimesteps, this->parameterList_);

            MultiVectorConstPtr_Type exportVector = this->problemTimeChem_->getSolution()->getBlock(0);

            exporterChem_->addVariable( exportVector, varName, "Scalar", 1, dom->getMapUnique() );
        }

        if (this->verbose_)
        {
            std::cout << "done -- " << std::endl;
        }
    }
    else
        reAssemble(type);

}
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::solveChemistryProblem() const
{
    this->updateChemInTime();

    // Move Mesh operation was already called for chemistry domain
    SmallMatrix<double> massCoeffChem(1);
    SmallMatrix<double> problemCoeffChem(1);
    SmallMatrix<int> defChem(1);
    double dt = timeSteppingTool_->get_dt();

    // Always bdf 1
    defChem[0][0] = 1;
    massCoeffChem[0][0] =  1. / dt;
    problemCoeffChem[0][0] =1.0; //
    double  coeffSourceTermChem = 1.0;

    this->problemTimeChem_->setTimeDef(defChem);
    this->problemTimeChem_->setTimeParameters(massCoeffChem,problemCoeffChem);
    // 1. Assemble Mass System
    MatrixPtr_Type massmatrix;
    this->setChemMassmatrix(massmatrix);
    // 1. Assemble Chemisty Problem
    this->problemTimeChem_->assemble();
    // 2. Rhs
    this->computeChemRHSInTime();

    //this->residualVec_->getBlockNonConst(0)->update(1.,*this->rhs_->getBlockNonConst(0),-1.);

    //this->problemTimeChem_->getRhs()->getBlockNonConst(0)->scale(-1.0);
    this->problemTimeChem_->combineSystems();


    this->problemTimeChem_->setBoundaries(timeSteppingTool_->currentTime());

    double iter= 0.;

    iter = this->problemTimeChem_->solve(this->problemTimeChem_->getRhs());

    exporterIterationsChem_->exportData(  timeSteppingTool_->currentTime() , iter);

    if (!exporterChem_.is_null())
            this->exporterChem_->save( this->timeSteppingTool_->currentTime() );

}
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::reAssemble(std::string type) const
{

    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);

    if (type == "MoveMesh"){
        moveMesh();
        return;
    }
    else if (type == "UpdateMeshDisplacement"){
        updateMeshDisplacement();
        return;
    }

    else if (type == "ComputeSolidRHSInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (ComputeSolidRHSInTime)" << '\n';

        computeSolidRHSInTime();
        return;
    }
    else if (type == "ComputeChemRHSInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (ComputeChemRHSInTime)" << '\n';

        computeChemRHSInTime();
        return;
    }
    else if (type == "UpdateChemInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (UpdateChemInTime)" << '\n';
        updateChemInTime();
        return;
    }

    else if(type == "UpdateTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateTime)" << '\n';

        updateTime();
        return;
    }
     else if(type == "SolveChemistryProblem")
    {
        if(this->verbose_)
            std::cout << "-- SolveChemistryProblem" << '\n';

        solveChemistryProblem();
        return;
    }
    else if(type == "SetBoundaries")
    {
        if(this->verbose_)
            std::cout << "-- set Boundaries" << '\n';

        setBoundariesSubProblems();
        return;
    }

    if(type == "FixedPoint")
    {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "should always be called by the residual and never here.");
    }
    else if(type == "Newton")
    {

        if(this->verbose_)
            std::cout << "-- reassemble Newton " << '\n';

        // Maybe nothing should happen here as there are no constant matrices
        this->system_.reset(new BlockMatrix_Type(2));
        if (chemistryExplicit_)
            this->system_.reset(new BlockMatrix_Type(1));


        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type B(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type BT(new Matrix_Type(this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type C(new Matrix_Type( this->getDomain(1)->getMapUnique(),this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ));

        // For implicit the system is ordered differently with solid block in 0,0 and diffusion in 1,1

        // As the implementation is based on a 2x2 System we use a tmp system for the actual assembly routine.
        BlockMatrixPtr_Type systemTmp(new BlockMatrix_Type(2));
        systemTmp->addBlock(A,0,0);
        systemTmp->addBlock(BT,0,1);
        systemTmp->addBlock(B,1,0);
        systemTmp->addBlock(C,1,1);

        MultiVectorConstPtr_Type c;
        if(chemistryExplicit_)
            c= this->problemTimeChem_->getSolution()->getBlock(0);
        else
            c = this->solution_->getBlock(1);

        c_rep_->importFromVector(c, true);

        MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
        d_rep_->importFromVector(d, true);

        this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,systemTmp,this->residualVec_, this->parameterList_, "Jacobian", true/*call fillComplete*/);
        //this->feFactory_->globalAssembly(materialModel_, this->dim_, 2, blockSol, this->system_, this->residualVec_,this->parameterList_,"Jacobian",true);

        if(chemistryExplicit_){
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->systemC_->addBlock(systemTmp->getBlock(1,1),0,0);
        }
        else{
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->system_->addBlock(systemTmp->getBlock(0,1),0,1);
            this->system_->addBlock(systemTmp->getBlock(1,0),1,0);
            this->system_->addBlock(systemTmp->getBlock(1,1),1,1);
        }
        if(nonlinearExternalForce_)
            computeSolidRHSInTime(false, true);

        if(this->verbose_)
            std::cout << " done -- " << '\n';


    }

}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{

    MultiVectorConstPtr_Type c;
    if(chemistryExplicit_)
        c= this->problemTimeChem_->getSolution()->getBlock(0);
    else
        c = this->solution_->getBlock(1);

    c_rep_->importFromVector(c, true);

    MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
    d_rep_->importFromVector(d, true);

    BlockMultiVectorPtr_Type resTmp(new BlockMultiVector_Type(2));
    MultiVectorPtr_Type res_d = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    MultiVectorPtr_Type res_c = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );
    resTmp->addBlock(res_d,0);
    resTmp->addBlock(res_c,1);

    this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,this->system_,resTmp, this->parameterList_, "Rhs", true/*call fillComplete*/);

    this->residualVec_->addBlock(Teuchos::rcp_const_cast<MultiVector_Type>(resTmp->getBlock(0)),0);

    if(!chemistryExplicit_){
        this->residualVec_->addBlock(Teuchos::rcp_const_cast<MultiVector_Type>(resTmp->getBlock(1)),1);
    }

    if(nonlinearExternalForce_)
        computeSolidRHSInTime(false);

    if (!type.compare("standard")){
        //this->rhs_->getBlockNonConst(0)->scale(-1.0);
        if(this->verbose_)
            std::cout << " Residual Type : " << type  << std::endl;
        this->residualVec_->getBlockNonConst(0)->update(-1.,*this->rhs_->getBlockNonConst(0),1.);
        //if ( !this->problemTimeStructure_->getSourceTerm()->getBlock(0).is_null() )
        //   this->residualVec_->getBlockNonConst(0)->update(-1.,*this->problemTimeStructure_->getSourceTerm()->getBlockNonConst(0),1.);

    }
    else if(!type.compare("reverse")){
        if(!chemistryExplicit_)
            this->residualVec_->getBlockNonConst(1)->scale(-1.0);

        if(this->verbose_)
            std::cout << " Residual Type : " << type  << std::endl;

        this->residualVec_->getBlockNonConst(0)->update(1.,*this->rhs_->getBlockNonConst(0),-1.);
        //if ( !this->problemTimeStructure_->getSourceTerm()->getBlock(0).is_null() )
        //     this->residualVec_->getBlockNonConst(0)->update(1.,*this->problemTimeStructure_->getSourceTerm()->getBlockNonConst(0),1.);

    }

    /*this->problemChem_->assemble();

    MultiVectorPtr_Type residualChemSCI =  Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(1) );
    this->problemChem_->getSystem()->getBlock(0,0)->apply( *this->problemChem_->getSolution()->getBlock(0), *residualChemSCI, Teuchos::NO_TRANS, -1. ); // y= -Ax + 0*y
    MultiVectorPtr_Type resChemNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(1) );
    resChemNonConst->update(1., *this->problemChem_->getRhs()->getBlock(0), 1.);
    */
    if(!chemistryExplicit_){
        MultiVectorPtr_Type resChemNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(1) );
        resChemNonConst->update(1., *this->problemChem_->getRhs()->getBlock(0), 1.);
    }



   // might also be called in the sub calculateNonLinResidualVec() methods which were used above
    if (type == "standard")
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    else if (type == "reverse"){
        //this->residualVec_->scale(-1.);
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
    }

    this->setBoundariesRHS(this->timeSteppingTool_->currentTime());


    Teuchos::Array<SC> norm_d(1);
    this->residualVec_->getBlock(0)->norm2(norm_d);

    if(this->verbose_)
        std::cout << "2-Norm of residual of displacement: " << norm_d[0] << std::endl;

    if(!chemistryExplicit_){
        Teuchos::Array<SC> norm_c(1);
        this->residualVec_->getBlock(1)->norm2(norm_c);

        if(this->verbose_)
            std::cout << "2-Norm of residual of concentration: " << norm_c[0] << std::endl;

    }

}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateMeshDisplacement() const
{

     *meshDisplacementOld_rep_ = *meshDisplacementNew_rep_;

}

// Share the coupled vectors with the solid and chemistry subproblems.
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setFromPartialVectorsInit() const
{

    //Chem
    if(!chemistryExplicit_){
        this->solution_->addBlock( this->problemChem_->getSolution()->getBlockNonConst(0), 1);
        //this->residualVec_->addBlock( this->problemChem_->getResidualVector()->getBlockNonConst(0), 1 );
        this->rhs_->addBlock( this->problemChem_->getRhs()->getBlockNonConst(0), 1 );
        this->sourceTerm_->addBlock( this->problemChem_->getSourceTerm()->getBlockNonConst(0), 1 );
    }


        this->solution_->addBlock( this->problemStructureNonLin_->getSolution()->getBlockNonConst(0), 0);
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0), 0 );
        this->rhs_->addBlock( this->problemStructureNonLin_->getRhs()->getBlockNonConst(0), 0 );
        this->previousSolution_->addBlock( this->problemStructureNonLin_->getPreviousSolution()->getBlockNonConst(0), 0 );
        this->sourceTerm_->addBlock( this->problemStructureNonLin_->getSourceTerm()->getBlockNonConst(0), 0 );


}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListStructure) const
{
    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblems " << std::endl;

    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();

    int sizeChem = this->problemChem_->getSystem()->size();
    const int sizeStructure = this->problemStructureNonLin_->getSystem()->size();

    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Chem ";

    problemTimeChem_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemChem_, this->comm_));

    if(this->verbose_)
        std::cout << "-- done " << std::endl;

    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Elasticity " ;

    problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructureNonLin_, this->comm_));

    if(this->verbose_)
        std::cout << "-- done " << std::endl;

    // ######################
    // Chem: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffChem(sizeChem);
    SmallMatrix<double> problemCoeffChem(sizeChem);
    SmallMatrix<int> defChem(sizeChem);

    double coeffSourceTermChem = 0.0;

    if ( this->getParameterList()->sublist("Timestepping Parameter").get("Class","Multistep") == "Multistep" ) {
        defChem[0][0] = 1;
        massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt;
        problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1);
        coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);

        this->problemTimeChem_->setTimeDef(defChem);
        this->problemTimeChem_->setTimeParameters(massCoeffChem,problemCoeffChem);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Implement other SCI Chem time stepping than BDF.");
    }
    // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems

    SmallMatrix<double> massCoeffStructure(sizeStructure);
    SmallMatrix<double> problemCoeffStructure(sizeStructure);
    SmallMatrix<int> defStructure(sizeStructure);
    double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    // Koeffizient vor der Massematrix
    for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
            // i == j, da vektorwertige Massematrix blockdiagonal ist
            if((*defTS_)[i][j] == 1 && i == j) // Weil: (d_s,c) und timeStepDef_ von SCI
            {
                defStructure[i][j] = 1;
            // Vorfaktor der Massematrix in der LHS
                massCoeffStructure[i][j] = 1.0/(dt*dt*beta);
            }
            else
            {
                massCoeffStructure[i][j] = 0.;
            }
        }
    }

    // Die anderen beiden Koeffizienten
    for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            if((*defTS_)[i ][j] == 1)
            {
                problemCoeffStructure[i][j] =  1.0;
                // Der Source Term ist schon nach der Assemblierung mit der Dichte \rho skaliert worden
                coeffSourceTermStructure = 1.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!
            }
            else // Die steady-Systemmatrix ist nicht zwingend blockdiagonal
            {
                problemCoeffStructure[i][j] = 1.0;
            }
        }
    }
    this->problemTimeStructure_->setTimeDef(defStructure);
    this->problemTimeStructure_->setTimeParameters(massCoeffStructure,problemCoeffStructure);

    this->problemTimeChem_->assemble( "MassSystem" );
    this->problemTimeStructure_->assemble( "MassSystem" );

    if(this->verbose_)
        std::cout << "done -- \n" << std::endl;

}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer SCI combineSystems(), ggf nichtlinear. Mass Matrix same as for Chem.
    //######################
    int size = this->problemTimeChem_->getSystem()->size();

    this->problemTimeChem_->systemMass_.reset(new BlockMatrix_Type(size));
    {
        BlockMatrixPtr_Type massSystem = Teuchos::rcp(new BlockMatrix_Type(1));

        massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        // 1 = Chem
        this->feFactory_->assemblyMass( this->dim_, this->problemTimeChem_->getFEType(0), "Scalar",  massmatrix, 1, true );
        //massSystem->addBlock(massmatrix,0,0);
       //this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2,1,this->dim_,c_rep_,d_rep_,massSystem,this->residualVec_, this->parameterList_, "MassMatrix", true/*call fillComplete*/);

        /*if(chemistryExplicit_){
            massmatrix->resumeFill();
            massmatrix->scale(0.);
            massmatrix->fillComplete( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->problemTimeChem_->getDomain(0)->getMapUnique() );
        }*/
        this->problemTimeChem_->systemMass_->addBlock(massmatrix, 0, 0);
    }
}

// The supplied staggered chemistry scheme uses backward Euler.
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeChemRHSInTime( ) const
{
    //######################
    // RHS nach BDF1
    //######################

    double dt = timeSteppingTool_->get_dt();
    double dt_prev = timeSteppingTool_->get_dt_prev();

    vec_dbl_Type coeffPrevSteps(1);
    coeffPrevSteps.at(0) = 1. / dt;
    SmallMatrix<double> massCoeff(1);
    SmallMatrix<double> problemCoeff(1);
    massCoeff[0][0] = 1. / dt;
    problemCoeff[0][0] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat
    this->problemTimeChem_->setTimeParameters(massCoeff, problemCoeff);

    this->problemTimeChem_->updateMultistepRhsFSI(coeffPrevSteps,1);/*apply (mass matrix_t / dt) to u_t and more*/


    // TODO
    if (this->problemTimeChem_->hasSourceTerm()) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Check sourceterm.");
    }

    // Wieder zu den eigentlichen Parametern zuruecksetzen, nachdem die temporaeren
    // genommen wurden.

}


// This is equivalent to the FSI Structure part.
// We add one additionale feature where we can distinguish between load or timestep
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeSolidRHSInTime(bool updateHistory, bool updateJacobian) const {
    //######################
    // RHS nach Newmark
    //######################
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();

    //double density = this->problemTimeStructure_->getParameterList()->sublist("Parameter Solid").get("Density",1.e-0);

    // Temporaerer Koeffizienten fuer die Skalierung der Massematrix in der rechten Seite des Systems in UpdateNewmarkRhs()
    vec_dbl_Type coeffTemp(1);
    coeffTemp.at(0) = 1.0;

    if( loadStepping_ == false){
        // Update u und berechne u' und u'' mit Hilfe der Newmark-Vorschrift
        if (updateHistory)
            this->problemTimeStructure_->updateSolutionNewmarkPreviousStep(timeSteppingTool_->get_dt_prev(), beta, gamma);

        // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
        // Bei Newmark lautet dies:
        // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
        // wobei u' = v (velocity) und u'' = w (acceleration).
        this->problemTimeStructure_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
    }
    else{
        // The if condition resets the rhs. If we skip it when we attemp loadstepping, the rhs will be updated continously and wrongly increase with each timestep
        this->problemTimeStructure_->getRhs()->scale(0.0);
    }

    //can we get rid of this?

    // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
    // if(time == 0){nur dann konstanten SourceTerm berechnen}

    if (this->problemTimeStructure_->hasSourceTerm())
    {
        if(externalForce_){

            MultiVectorPtr_Type FERhs = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));
            /*vec_dbl_Type funcParameter(6,0.);
            funcParameter[0] = timeSteppingTool_->t_;
            // how can we use different parameters for different blocks here?
            funcParameter[1] =this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Volume force",0.00211);
            funcParameter[2]= this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Load Step Size",1.);
            funcParameter[3] = this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Load Ramp End",1.);

            funcParameter[4] =this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Heart Beat Start",70.);
            funcParameter[5] = 0.;*/
            vec_dbl_Type funcParameter(1,0.);
            funcParameter[0] = timeSteppingTool_->t_;
            // how can we use different parameters for different blocks here?
            const auto& rhsParameters = this->problemTimeStructure_->getUnderlyingProblem()->parasSourceFunc_;
            funcParameter.insert(funcParameter.end(), rhsParameters.begin(), rhsParameters.end());
            funcParameter.push_back(0.);


            if(nonlinearExternalForce_){

                MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
                d_rep_->importFromVector(d, true);
                MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow()*2 ) );

                this->feFactory_->assemblyNonlinearSurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, d_rep_,Kext, funcParameter, this->problemTimeStructure_->getUnderlyingProblem()->rhsFuncVec_[0],this->parameterList_);

                if (updateJacobian) {
                    // Residual evaluations must not accumulate load tangents.
                    MatrixPtr_Type AKext(new Matrix_Type(this->getDomain(0)->getMapVecFieldUnique(),
                        this->dim_ * this->getDomain(0)->getApproxEntriesPerRow()));
                    this->system_->getBlock(0,0)->addMatrix(1., AKext, 0.);
                    Kext->addMatrix(-1., AKext, 1.);
                    AKext->fillComplete(this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
                    this->system_->addBlock(AKext, 0, 0);
                }

            }
            else
                this->feFactory_->assemblySurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, d_rep_, funcParameter, this->problemTimeStructure_->getUnderlyingProblem()->rhsFuncVec_[0],this->parameterList_);


            this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs, false, "Add" );



        }
        else
            this->problemTimeStructure_->assembleSourceTerm( timeSteppingTool_->t_ );
        //this->problemTimeStructure_->getSourceTerm()->scale(density);
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!

        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpSourceterm = Teuchos::rcp(new BlockMultiVector_Type(1)) ;
        tmpSourceterm->addBlock(this->sourceTerm_->getBlockNonConst(0),0);

        BlockMultiVectorPtr_Type tmpPtr ;
        if(externalForce_)
            tmpPtr =  tmpSourceterm;
        else
            tmpPtr = this->problemTimeStructure_->getSourceTerm();



        this->problemTimeStructure_->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);
        this->rhs_->addBlock( this->problemTimeStructure_->getRhs()->getBlockNonConst(0), 0 );


    }

}

// Can stay the same
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix
    //######################
    double density = this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Density",1.0);
    int size = this->problemTimeStructure_->getSystem()->size();

    // bool restart = this->parameterList_->sublist("Timestepping Parameter").get("Restart", false);
    // double timeStepRestart = this->parameterList_->sublist("Timestepping Parameter").get("Time step", 0.0);

    // if(timeSteppingTool_->currentTime() == 0.0 || (restart &&  timeSteppingTool_->currentTime() -1.e-5 < timeStepRestart ))
    // {
        this->problemTimeStructure_->systemMass_.reset(new BlockMatrix_Type(size));
        {

            massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // 0 = Struktur
            this->feFactory_->assemblyMass(this->dim_, this->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 0, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
            if(loadStepping_ == true)
                massmatrix->scale(0.0);

            massmatrix->fillComplete( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique());

            this->problemTimeStructure_->systemMass_->addBlock( massmatrix, 0, 0 );
        }
    // }
}

// Can stay the same
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setBoundariesSubProblems( ) const
{
    //######################
    // Boundaries
    //######################

    this->problemTimeStructure_->setBoundaries();
    this->problemTimeChem_->setBoundaries();


}


// We update the time in multiple component. Also the elements need to be updated with respect to time
// since internal processes like 'active response' are triggered after certain time
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->advanceTime();
    timeSteppingTool_->prepareStep();

    // std::cout << " ###### Timestep in SCI dt_prev" << timeSteppingTool_->dt_prev_ << " dt= " << timeSteppingTool_->dt_ <<" time= " << timeSteppingTool_->t_ << " ####### " << std::endl;

    MultiVectorConstPtr_Type c;
    if(chemistryExplicit_)
        c= this->problemTimeChem_->getSolution()->getBlock(0);
    else
        c = this->solution_->getBlock(1);

    MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
    d_rep_->importFromVector(d, true);

    c_rep_->importFromVector(c, true);

    this->feFactory_->advanceInTimeAssemblyFEElements(timeSteppingTool_->dt_, d_rep_, c_rep_ );

    this->problemTimeChem_->updateTime(timeSteppingTool_->t_);
    this->problemTimeStructure_->updateTime(timeSteppingTool_->t_);

    this->newtonStep_ = 0; // reset newton step counter at the beginning of each time step
}


// This is needed for explicit coupling
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;

    displacementUniqueConst = this->solution_->getBlock(0);
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Chemistry
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSI-Domain, da selbes Objekt bei problemChem_ und problemTimeChem_.
   // ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
   // ( Teuchos::rcp_const_cast<Domain_Type>(this->problemChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);

   // ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
   // ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeStructure_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);

}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::initializeCE(){
    // initializeProblem() initially creates both fields. In explicit chemistry
    // mode only displacement belongs to the coupled nonlinear solve.
    if (chemistryExplicit_) {
        this->solution_->resize( 1 );
        this->rhs_->resize( 1 );
        this->sourceTerm_->resize( 1 );
        this->rhsFuncVec_.resize( 1 );
        this->previousSolution_->resize( 1 );
        this->residualVec_->resize( 1 );
        this->initVectorSpaces();  //reinitialize NOX vector spaces
    }
}


// Am Anfang der Zeititeration erst updateSolutionMultiPreviousStep() aufrufen und dann erst updateMultistepRhs(),
// damit die previousSolution_ initialisiert sind. Genauso fuer SystemMass
// TODO: updateSystemMassMultiPreviousStep() fertig programmieren
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateChemInTime() const
{
    // int nmbBDF = 1;
    // FEDDLIB_NOTIFICATION("updateChemInTime", this->getComm()->getRank() == 0, " This was changed to accomodate also BDF 2 -- Keep in mind, AceGEN does BDF-1");

    this->problemTimeChem_->updateSolutionMultiPreviousStep(1);
    this->problemTimeChem_->updateSystemMassMultiPreviousStep(1);

}
template<class SC,class LO,class GO,class NO>
typename SCI<SC,LO,GO,NO>::BlockMultiVectorPtr_Type SCI<SC,LO,GO,NO>::getPostProcessingData() const
{
    if (postProcessingnames_.empty())
        return Teuchos::rcp(new BlockMultiVector_Type(0));
    d_rep_->importFromVector(this->solution_->getBlock(0), true);
    c_rep_->importFromVector(chemistryExplicit_ ? this->problemTimeChem_->getSolution()->getBlock(0)
                                               : this->solution_->getBlock(1), true);
    this->feFactory_->updateSolutionAssemblyFEElements(d_rep_, c_rep_);
    const int columns[] = {10, 11, 12, 13, 23, 39, 40, 41, 42, 43, 44, 27, 28, 45, 46, 47, 48};
    auto result = Teuchos::rcp(new BlockMultiVector_Type(postProcessingnames_.size()));
    for (UN i = 0; i < postProcessingnames_.size(); ++i) {
        auto field = Teuchos::rcp(new MultiVector_Type(this->getDomain(0)->getMapUnique(), 1));
        this->feFactory_->postProcessing(columns[i], field);
        result->addBlock(field, i);
    }
    return result;
}

template<class SC,class LO,class GO,class NO>
vec_string_Type SCI<SC,LO,GO,NO>::getPostprocessingNames() const
{
    return postProcessingnames_;
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::getValuesOfInterest(BlockMultiVectorPtr_Type&)
{
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Legacy SCI history file I/O is unavailable in this branch.");
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::exportValuesOfInterest(double)
{
    TEUCHOS_TEST_FOR_EXCEPTION(this->parameterList_->sublist("Timestepping Parameter").get("Export history", false),
        std::logic_error, "Legacy SCI history file export is unavailable in this branch.");
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::importValuesOfInterest(double)
{
    TEUCHOS_TEST_FOR_EXCEPTION(this->parameterList_->sublist("Timestepping Parameter").get("Import history", false),
        std::logic_error, "Legacy SCI history file import is unavailable in this branch.");
}

}
#endif
