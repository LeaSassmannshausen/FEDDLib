#ifndef DAESOLVERINTIME_DEF_hpp
#define DAESOLVERINTIME_DEF_hpp
#include "DAESolverInTime_decl.hpp"


/*!
 Definition of DAESolverInTime

 @brief  DAESolverInTime
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {


template<class SC,class LO,class GO,class NO>
DAESolverInTime<SC,LO,GO,NO>::DAESolverInTime( CommConstPtr_Type comm ):
comm_(comm),
verbose_(comm->getRank()==0),
parameterList_(),
isTimeSteppingDefined_(false),
problem_(),
problemTime_(),
timeStepDef_(),
timeSteppingTool_(),
exporter_vector_(),
export_solution_vector_(),
boolExporterSetup_(false),
boolExporterSetupPostprocess_(false)
#ifdef FEDD_TIMER
,solveProblemTimer_ (Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime - Solve Problem"))
#endif
#ifdef FEDD_DETAIL_TIMER
,reassmbleAddInterfaceRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Interface RHS"))
,reassmbleUpdateMeshDisplacementTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Mesh Displ."))
,reassmbleSolveGeometryTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Solve Geometry"))
,reassmbleMoveMeshTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Move Mesh"))
,reassmbleSolidMassAndRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Solid Massmatrix and RHS"))
,reassmbleForTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Fluid Problem (steady Navier-Stokes)"))
,reassmbleUpdateFluidInTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Fluid Solution"))
#endif
{


}

template<class SC,class LO,class GO,class NO>
DAESolverInTime<SC,LO,GO,NO>::DAESolverInTime( ParameterListPtr_Type &parameterList, CommConstPtr_Type comm):
comm_(comm),
verbose_(comm->getRank()==0),
parameterList_(parameterList),
isTimeSteppingDefined_(false),
problem_(),
problemTime_(),
timeStepDef_(),
timeSteppingTool_(),
exporter_vector_(),
export_solution_vector_(),
boolExporterSetup_(false),
boolExporterSetupPostprocess_(false)
#ifdef FEDD_TIMER
,solveProblemTimer_ (Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime - Solve Problem"))
#endif
#ifdef FEDD_DETAIL_TIMER
,reassmbleAddInterfaceRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Interface RHS"))
,reassmbleUpdateMeshDisplacementTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Mesh Displ."))
,reassmbleSolveGeometryTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Solve Geometry"))
,reassmbleMoveMeshTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Move Mesh"))
,reassmbleSolidMassAndRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Solid Massmatrix and RHS"))
,reassmbleForTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Fluid Problem (steady Navier-Stokes)"))
,reassmbleUpdateFluidInTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Fluid Solution"))
#endif
    
{

}

template<class SC,class LO,class GO,class NO>
DAESolverInTime<SC,LO,GO,NO>::DAESolverInTime( SmallMatrix<int> &timeStepDef, ParameterListPtr_Type &parameterList, CommConstPtr_Type comm):
comm_(comm),
verbose_(comm->getRank()==0),
parameterList_(parameterList),
isTimeSteppingDefined_(false),
problem_(),
problemTime_(),
timeStepDef_(),
timeSteppingTool_(),
exporter_vector_(),
export_solution_vector_(),
boolExporterSetup_(false),
boolExporterSetupPostprocess_(false)
#ifdef FEDD_TIMER
,solveProblemTimer_ (Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime - Solve Problem"))
#endif
#ifdef FEDD_DETAIL_TIMER
,reassmbleAddInterfaceRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Interface RHS"))
,reassmbleUpdateMeshDisplacementTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Mesh Displ."))
,reassmbleSolveGeometryTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Solve Geometry"))
,reassmbleMoveMeshTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Move Mesh"))
,reassmbleSolidMassAndRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Solid Massmatrix and RHS"))
,reassmbleForTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Assemble Fluid Problem (steady Navier-Stokes)"))
,reassmbleUpdateFluidInTimeTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - DAETime Detail - FSI Update Fluid Solution"))
#endif
    
{
    this->defineTimeStepping(timeStepDef);
}

template<class SC,class LO,class GO,class NO>
DAESolverInTime<SC,LO,GO,NO>::~DAESolverInTime(){

}


template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setProblem(Problem_Type& problem){

    problem_ = Teuchos::rcpFromRef(problem); /*now a NON-OWNING TEUCHOS::RCP to the object which was probably constructed in main function */

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::defineTimeStepping(SmallMatrix<int> &timeStepDef){

    timeStepDef_ = timeStepDef;
    timeSteppingTool_.reset(new TimeSteppingTools(sublist(parameterList_,"Timestepping Parameter") , comm_));
    isTimeSteppingDefined_ = true;


}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTime(){

    TEUCHOS_TEST_FOR_EXCEPTION(problemTime_.is_null(), std::logic_error, "Mass system is null.");

    checkTimeSteppingDef();

    if(this->parameterList_->sublist("Parameter").get("FSI",false))
    {
        advanceInTimeFSI();
    }
    else if(this->parameterList_->sublist("Parameter").get("SCI",false) && !parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("Multistep"))
    {
        advanceInTimeSCI();
        
    }
    else if(this->parameterList_->sublist("Parameter").get("FSCI",false))
    {
        advanceInTimeFSCI();
        
    }
    else{
        if (!parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("Loadstepping")) {
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
            if (nonLinProb.is_null()){ 
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Loadstepping only available to nonlinear Problems. (It is a tool to increase Nonlinear Solver convergence.)");
            }
            else{
                advanceWithLoadStepping();
            }
        }
        else if (!parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("Singlestep")) {
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
            if (nonLinProb.is_null()) {
                advanceInTimeLinear();
            }
            else{
                advanceInTimeNonLinear();
            }
        }
        else if(!parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("Newmark"))
        {
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
            if (nonLinProb.is_null()) {
                advanceInTimeLinearNewmark();
            }
            else
            {
                 advanceInTimeNonLinearNewmark();
            }
        }
        else if( !parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("Multistep")  ){
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
            if (nonLinProb.is_null()) {
                advanceInTimeLinearMultistep();
            }
            else{
                advanceInTimeNonLinearMultistep();
            }
        }
        else if( !parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep").compare("External")  ){
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
            if (nonLinProb.is_null()) {
                advanceInTimeLinearExternal();
            }
            else{
                advanceInTimeNonLinearExternal();
            }
        }
        else{
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown time discretization type.");
        }
    }
    
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeLinear(){

    
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print)
        exportTimestep();
    bool fullImplicitPressure = false;
    int size = timeStepDef_.size();
    double dt = 0.0;
    while (timeSteppingTool_->continueTimeStepping()) {
        problemTime_->updateTime ( timeSteppingTool_->currentTime() );
        dt = timeSteppingTool_->get_dt();
        problemTime_->updateSolutionPreviousStep();
        
        Teuchos::Array<BlockMatrixPtr_Type> bMatNonLin_vec_allstages;
        BlockMultiVectorPtrArray_Type	solutionRK_stages;
        BlockMultiVectorPtrArray_Type	sourceTermRK_stages;
        
        for (int s=0; s<timeSteppingTool_->getNmbStages(); s++) {
            double time = timeSteppingTool_->currentTime() + dt * timeSteppingTool_->getButcherTableC(s);
            if (verbose_)
                cout << "Currently in stage " << s+1 << " of "<< timeSteppingTool_->getNmbStages() << endl;
            
            problemTime_->updateRhs();/*apply (mass matrix / dt) to u_t*/
            if ( problemTime_->hasSourceTerm() ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Fix source term.");
                problemTime_->assembleSourceTerm( time );
            }
            
            if (s==0 && timeSteppingTool_->getButcherTableCoefficient(s , s) == 0.0) {
                /* solution of last time step is later added to rk solution vector */
            }
            else if (s==0 && timeSteppingTool_->getButcherTableCoefficient(s , s) != 0.0) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented butcher table!");
            }
            else{
                for (int s_prior=0; s_prior<s; s_prior++) {
                    SmallMatrix<double> coeff(size);
                    for (int i=0; i<size; i++) {
                        for (int j=0; j<size; j++) {
                            if (timeStepDef_[i][j]>0) {
                                if (i==0 && j==1 && fullImplicitPressure) {
                                    coeff[i][j] = 0;
                                }
                                else{
                                    coeff[i][j] = - timeSteppingTool_->getButcherTableCoefficient(s , s_prior);
                                }
                            }
                        }
                    }
                    if (problemTime_->hasSourceTerm())
                        this->addRhsDAE(coeff, sourceTermRK_stages[s_prior] );
                    
                    this->addRhsDAE( coeff, problemTime_->getSystem(), solutionRK_stages[s_prior] );
                    
                }
                
                SmallMatrix<double> massCoeff(size);
                SmallMatrix<double> problemCoeff(size);
                double coeffSourceTerm = 0.;
                for (int i=0; i<size; i++) {
                    for (int j=0; j<size; j++) {
                        if (timeStepDef_[i][j]>0 && i==j) {
                            massCoeff[i][j] = 1. / dt;
                        }
                        else if (timeStepDef_[i][j]==2) /*force off-diagnonal mass matrix*/
                            massCoeff[i][j] = 1. / dt;
                        else{
                            massCoeff[i][j] = 0.;
                        }
                    }
                }
                for (int i=0; i<size; i++) {
                    for (int j=0; j<size; j++){
                        if (timeStepDef_[i][j]>0){
                            if (i==0 && j==1 && fullImplicitPressure) {
                                problemCoeff[i][j] = 1.;
                                coeffSourceTerm = timeSteppingTool_->getButcherTableCoefficient(s , s); // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT
                            }
                            else{
                                problemCoeff[i][j] = timeSteppingTool_->getButcherTableCoefficient(s , s);
                                coeffSourceTerm = timeSteppingTool_->getButcherTableCoefficient(s , s); // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!
                            }
                        }
                        else{
                            problemCoeff[i][j] = 1.;
                        }
                    }
                }
                
                if (problemTime_->hasSourceTerm())
                    addSourceTermToRHS(coeffSourceTerm);
                
                problemTime_->setTimeParameters(massCoeff, problemCoeff);
                problemTime_->combineSystems();
                problemTime_->setBoundaries(time);
                
                problemTime_->solve();

            }
            
            if (s+1 == timeSteppingTool_->getNmbStages()) {
                BlockMultiVectorPtr_Type tmpSolution = Teuchos::rcp( new BlockMultiVector_Type ( problemTime_->getSolution() ) );
                solutionRK_stages.push_back(tmpSolution);
                BlockMultiVectorPtr_Type tmpSolutionPtr = problemTime_->getSolution();
                BlockMatrixPtr_Type tmpMassSystem = problemTime_->getMassSystem();
                timeSteppingTool_->calculateSolution( tmpSolutionPtr, solutionRK_stages, tmpMassSystem);
            }
            else{
                solutionRK_stages.push_back( problemTime_->getSolution() );
                if ( problemTime_->hasSourceTerm() )
                    sourceTermRK_stages.push_back( problemTime_->getSourceTerm() );
            }
        }
        
        timeSteppingTool_->advanceTime(true/*output info*/);
        
        if (print)
            exportTimestep();
    }
    
    if (print) {
        closeExporter();
    }
    if (parameterList_->sublist("NSParameter").get("Calculate Coefficients",false)) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "close txt exporters here.");
    }

}
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::getMassCoefficients(SmallMatrix<double> &massCoeff){
    int size = timeStepDef_.size();
    massCoeff.resize(size);
    
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (timeStepDef_[i][j]>0 && i==j)
                massCoeff[i][j] = 1.;
            else
                massCoeff[i][j] = 0.;
        }
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::getMultiStageCoefficients(SmallMatrix<double> &problemCoeff, int stage, int stagePrior, bool forRhs){
    int size = timeStepDef_.size();
    problemCoeff.resize(size);
    
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (timeStepDef_[i][j]>0)
                problemCoeff[i][j] = timeSteppingTool_->get_dt()* timeSteppingTool_->getButcherTableCoefficient(stage , stagePrior);
            else{
                if (forRhs)
                    problemCoeff[i][j] = 0.;
                else
                    problemCoeff[i][j] = 1.;
            }
        }
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::buildMultiStageRhs( int stage, Teuchos::Array<BlockMatrixPtr_Type>& matrixPrevStages, BlockMultiVectorPtrArray_Type& solutionPrevStages ){
    
    int size = timeStepDef_.size();
    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    
    getMassCoefficients(massCoeff);
    
    problemTime_->setTimeParameters(massCoeff, problemCoeff);//problemCoeff not needed here
    problemTime_->updateRhs();/*apply (mass matrix / dt) to u_t*/
    
    bool fullImplicitPressure = parameterList_->sublist("Timestepping Parameter").get("Full implicit pressure",false);
    for (int prevStage=0; prevStage<stage; prevStage++) {
        
        getMultiStageCoefficients(problemCoeff, stage, prevStage, true);
        if (fullImplicitPressure && size>1 )
            problemCoeff[0][1] = 0.;
        
        problemCoeff.scale(-1.);

        TEUCHOS_TEST_FOR_EXCEPTION(problemTime_->hasSourceTerm(), std::logic_error, "Using source term must be implemented for single-step methods.");
//                        addRhsDAE(coeff,beNL_vec_allstages.at(s_prior), solutionRK_stages.at(s_prior), sourceTermRK_stages.at(s_prior));
        addRhsDAE( problemCoeff, matrixPrevStages[prevStage], solutionPrevStages[prevStage] );

    }
}



template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeNonLinear(){

    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    BlockMultiVectorPtr_Type solShort;
    if (print){
        if (parameterList_->sublist("Timestepping Parameter").get("Print Solution Short",false)) {
            solShort = Teuchos::rcp(new BlockMultiVector_Type (problemTime_->getSolution()) );
            exportTimestep(solShort);
        }
        else{
            exportTimestep();
        }
    }
    
    // Navier-Stokes treatment of pressure
    bool fullImplicitPressure = parameterList_->sublist("Timestepping Parameter").get("Full implicit pressure",false);
    bool semiImplicitPressure = parameterList_->sublist("Timestepping Parameter").get("Semi implicit pressure",false);
    bool correctPressure = parameterList_->sublist("Timestepping Parameter").get("Correct pressure",false);
    int size = timeStepDef_.size();
    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    double dt = 0.0;
    int timeit = 0;
    while (timeSteppingTool_->continueTimeStepping()) {
        
        dt = timeSteppingTool_->get_dt();
        // Which values should we use for extrapolation? u_m and u_m-1 (current implementation) or should we use most recent U_i-1 of multi stages
        if(!parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation")) {
            if (timeSteppingTool_->currentTime()!=0.)
                problemTime_->updateSolutionMultiPreviousStep(2); //2nd order
            else
                problemTime_->updateSolutionMultiPreviousStep(1);
        }
        else{
            problemTime_->updateSolutionPreviousStep();
        }
        
        
//        BlockMultiVectorPtrArray_Type	sourceTermRK_stages;

        TEUCHOS_TEST_FOR_EXCEPTION(timeSteppingTool_->getButcherTableCoefficient(0 , 0) != 0.0, std::logic_error, "Not implemented butchertable! First stage should have 0 diagonal value");
        if (verbose_)
            cout << "Currently in stage " << 1 << " of "<< timeSteppingTool_->getNmbStages() <<" (dummy stage)"<< endl;
        // Multistage stepping, in general we use at least 2 stages (implicit Euler and Crank-Nicolson)
        Teuchos::Array<BlockMatrixPtr_Type> matrixPrevStages;
        BlockMultiVectorPtrArray_Type       solutionPrevStages;
        
        BlockMatrixPtr_Type blockMatrix = Teuchos::rcp( new BlockMatrix_Type( problemTime_->getSystem()->size() ) );
        problemTime_->reAssembleAndFill( blockMatrix ); // Reassemble FixedPoint
        
        matrixPrevStages.push_back( blockMatrix );
        
        BlockMultiVectorPtr_Type sol =
            Teuchos::rcp( new BlockMultiVector_Type( problemTime_->getSolution() ) );
        solutionPrevStages.push_back( sol );
        
        for (int stage=1; stage<timeSteppingTool_->getNmbStages(); stage++) {
            double time = timeSteppingTool_->currentTime() + dt * timeSteppingTool_->getButcherTableC(stage);
            problemTime_->updateTime( time );
            if (verbose_)
                cout << "Currently in stage " << stage+1 << " of "<< timeSteppingTool_->getNmbStages() << endl;
            
            buildMultiStageRhs( stage, matrixPrevStages, solutionPrevStages );
                                    
            TEUCHOS_TEST_FOR_EXCEPTION(problemTime_->hasSourceTerm(), std::logic_error, "Using source term must be implemented for single-step methods.");
//                problemTime_->AssembleSourceTerm(time);
            

            SmallMatrix<double> massCoeff(size);
            SmallMatrix<double> problemCoeff(size);
                        
            getMassCoefficients(massCoeff);
            getMultiStageCoefficients(problemCoeff, stage, stage, false);
            if (fullImplicitPressure && size>1 )
                problemCoeff[0][1] = timeSteppingTool_->get_dt();
                                                       
            problemTime_->setTimeParameters(massCoeff, problemCoeff);

            NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));
            nlSolver.solve(*problemTime_,time);
            
            if (correctPressure) {
                TEUCHOS_TEST_FOR_EXCEPTION( !fullImplicitPressure && !semiImplicitPressure, std::logic_error,"There is no pressure that can be corrected." );
                
                MultiVectorPtr_Type sol = problemTime_->getSolution()->getBlockNonConst(1);
                MultiVectorConstPtr_Type solLast = problemTime_->getSolutionPreviousTimestep()->getBlock(1);
                timeSteppingTool_->correctPressure( sol, solLast );
                if (verbose_)
                    std::cout << "Pressure corrected." << std::endl;
            }
            
            if (stage+1 == timeSteppingTool_->getNmbStages()) {
                // save last solution
                BlockMultiVectorPtr_Type tmpSolution = Teuchos::rcp( new BlockMultiVector_Type ( problemTime_->getSolution() ) );
                
                solutionPrevStages.push_back( tmpSolution );
                BlockMultiVectorPtr_Type tmpSolutionPtr = problemTime_->getSolution();
                BlockMatrixPtr_Type tmpMassSystem = problemTime_->getMassSystem();
                timeSteppingTool_->calculateSolution( tmpSolutionPtr, solutionPrevStages, tmpMassSystem, solShort);
            }
            else{
                BlockMatrixPtr_Type blockMatrix = Teuchos::rcp( new BlockMatrix_Type( problemTime_->getSystem()->size() ) );
                problemTime_->reAssembleAndFill( blockMatrix );
                matrixPrevStages.push_back( blockMatrix );
                
                BlockMultiVectorPtr_Type sol =
                    Teuchos::rcp( new BlockMultiVector_Type( problemTime_->getSolution() ) );
                solutionPrevStages.push_back( sol );
                
                TEUCHOS_TEST_FOR_EXCEPTION(problemTime_->hasSourceTerm(), std::logic_error, "Using source term must be implemented for single-step methods.");
//                    sourceTermRK_stages.push_back(*problemTime_->GetSourceTerm());
            }
        }
        timeSteppingTool_->advanceTime(true/*output info*/);
        timeit++;
        if (print) {
            exportTimestep();
        }
        if (parameterList_->sublist("NSParameter").get("Calculate Coefficients",false)) {
            vec_dbl_ptr_Type values(new vec_dbl_Type(4));
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Drag and Lift are not implemented.");
//            nonLinearProblem_->ComputeDragLift(values);
        }
    }

    if (print) {
        closeExporter();
    }
    if (parameterList_->sublist("NSParameter").get("Calculate Coefficients",false)) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "close txt exporters here.");
    }

}



template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeLinearNewmark()
{
    // problemCoeff vor A (= kompletes steady-System)
    // massCoeff vor M (= Massematrix)
    // coeffSourceTerm vor f (= rechte Seite der DGL)

    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print)
    {
        exportTimestep();
    }

    // Groesse des Problems, Zeitschrittweite und Newmark-Parameter
    int size = timeStepDef_.size();
    TEUCHOS_TEST_FOR_EXCEPTION( size>1, std::runtime_error, "Newmark only implemented for systems of size 1x1.");
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();

    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    double coeffSourceTerm = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    // Koeffizient vor der Massematrix
    massCoeff[0][0] = 1.0/(dt*dt*beta);
    problemCoeff[0][0] =  1.0;
    // Der Source Term ist schon nach der Assemblierung mit der Dichte \rho skaliert worden
    coeffSourceTerm = 1.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!

    // Temporaerer Koeffizienten fuer die Skalierung der Massematrix in der rechten Seite des Systems in UpdateNewmarkRhs()
    vec_dbl_Type coeffTemp(1);
    coeffTemp.at(0) = 1.0;

    // Uebergebe die Parameter fuer Masse- und steady-Problemmatrix an TimeProblem
    // Wegen moeglicher Zeitschrittweitensteuerung, rufe CombineSystems()
    // in jedem Zeitschritt auf, um LHS neu aufzustellen.
    // Bei AdvanceInTimeNonLinear... wird das in ReAssemble() gemacht!!!
    problemTime_->setTimeParameters(massCoeff, problemCoeff);

    // ######################
    // Time loop
    // ######################
    while(timeSteppingTool_->continueTimeStepping())
    {
        // Stelle (massCoeff*M + problemCoeff*A) auf
        problemTime_->combineSystems();

        // Update u und berechne u' und u'' mit Hilfe der Newmark-Vorschrift
        problemTime_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);

        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateTime ( timeSteppingTool_->currentTime() );
        // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
        // Bei Newmark lautet dies:
        // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
        // wobei u' = v (velocity) und u'' = w (acceleration).
        problemTime_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);

        // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
        // if(time == 0){nur dann konstanten SourceTerm berechnen}
        if ( problemTime_->hasSourceTerm() )
        {
//            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Fix source term.");
            problemTime_->assembleSourceTerm(time);

            // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
            // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
            addSourceTermToRHS(coeffSourceTerm);
        }

        // Uebergabeparameter fuer BC noch hinzu nehmen!
        problemTime_->setBoundaries(time);
        problemTime_->solve();

        timeSteppingTool_->advanceTime(true/*output info*/);

        if (print) {
            exportTimestep();
        }

    }

    comm_->barrier();
    if (print)
    {
        closeExporter();
    }
}


// TODO: Irgendwann einmal fuer St. Venant-Kirchoff oder so programmieren.
// Erstmal nicht wichtig
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeNonLinearNewmark()
{
    
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    bool printExtraData = parameterList_->sublist("General").get("Export Extra Data",false);
    bool printData = parameterList_->sublist("General").get("Export Data",false);
    if (print)
    {
        exportTimestep();
    }
    
    vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations
    ExporterTxtPtr_Type exporterTimeTxt;
    ExporterTxtPtr_Type exporterIterations;
    ExporterTxtPtr_Type exporterNewtonIterations;
    if (printData) {
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterTimeTxt = Teuchos::rcp(new ExporterTxt());
        exporterTimeTxt->setup( "time" + suffix, this->comm_ );
        
        exporterNewtonIterations = Teuchos::rcp(new ExporterTxt());
        exporterNewtonIterations->setup( "newtonIterations" + suffix, this->comm_ );
        
        exporterIterations = Teuchos::rcp(new ExporterTxt());
        exporterIterations->setup( "linearIterations" + suffix, this->comm_ );
    }
    // Groesse des Problems, Zeitschrittweite und Newmark-Parameter
    int size = timeStepDef_.size();
    TEUCHOS_TEST_FOR_EXCEPTION( size>1, std::runtime_error, "Newmark only implemented for systems of size 1x1.");
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    double coeffSourceTerm = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren
    
    // Koeffizient vor der Massematrix
    massCoeff[0][0] = 1.0/(dt*dt*beta);
    problemCoeff[0][0] =  1.0;
    // Der Source Term ist schon nach der Assemblierung mit der Dichte \rho skaliert worden
    coeffSourceTerm = 1.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!
    
    // Temporaerer Koeffizienten fuer die Skalierung der Massematrix in der rechten Seite des Systems in UpdateNewmarkRhs()
    vec_dbl_Type coeffTemp(1);
    coeffTemp.at(0) = 1.0;
    
    // Uebergebe die Parameter fuer Masse- und steady-Problemmatrix an TimeProblem
    // Wegen moeglicher Zeitschrittweitensteuerung, rufe CombineSystems()
    // in jedem Zeitschritt auf, um LHS neu aufzustellen.
    problemTime_->setTimeParameters(massCoeff, problemCoeff);
    
    // ######################
    // Time loop
    // ######################
    while(timeSteppingTool_->continueTimeStepping())
    {
        cout << "  ############## Timeloop Newmark ##########" << endl;
        // Stelle (massCoeff*M + problemCoeff*A) auf
        problemTime_->combineSystems();
        
        // Update u und berechne u' und u'' mit Hilfe der Newmark-Vorschrift
        problemTime_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);
        
        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateTime ( time );
        // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
        // Bei Newmark lautet dies:
        // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
        // wobei u' = v (velocity) und u'' = w (acceleration).
        cout << " Beta " << beta << " gamma " << gamma << endl;
        problemTime_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
        
        // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
        // if(time == 0){nur dann konstanten SourceTerm berechnen}
        if (problemTime_->hasSourceTerm())
        {
//            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Fix source term.");
            problemTime_->assembleSourceTerm(time);
            
            // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
            // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
            // This is done in  calculateNonLinResidualVec for this nonlinear problem
            addSourceTermToRHS(coeffSourceTerm);
        }
        
        // Uebergabeparameter fuer BC noch hinzu nehmen!
//        problemTime_->setBoundaries(time);
        
        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","Newton"));
        nlSolver.solve( *problemTime_, time, its );
        
        timeSteppingTool_->advanceTime(true/*output info*/);
        if (printData) {
            exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
            exporterIterations->exportData( (*its)[0] );
            exporterNewtonIterations->exportData( (*its)[1] );
        }
        if (print) {
            exportTimestep();
        }
        
    }
    
    comm_->barrier();
    if (print)
    {
        closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeSCI()
{
    // problemCoeff vor A (= komplettes steady-System)
    // massCoeff vor M (= Massematrix)
    // coeffSourceTerm vor f (= rechte Seite der DGL)
    
    SCIProblemPtr_Type sci = Teuchos::rcp_dynamic_cast<SCIProblem_Type>( this->problemTime_->getUnderlyingProblem() );
    
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    bool printStress = parameterList_->sublist("General").get("Export Stress",false);
    bool printData = parameterList_->sublist("General").get("Export Data",true);
    bool printExtraData = parameterList_->sublist("General").get("Export Extra Data",false);

    BlockMultiVectorPtr_Type stressVec;
    if (print)
    {
        exportTimestep();
    }

    vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations
    ExporterTxtPtr_Type exporterTimeTxt;
    ExporterTxtPtr_Type exporterDisplXTxt;
    ExporterTxtPtr_Type exporterDisplYTxt;
    ExporterTxtPtr_Type exporterIterations;
    ExporterTxtPtr_Type exporterNewtonIterations;
    ExporterTxtPtr_Type exporterCornerValue;

    GO idExport = parameterList_->sublist("General").get("ExportID",0);
    LO valueCorner =-1;

    
    if (printData) {
        exporterTimeTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplXTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplYTxt = Teuchos::rcp(new ExporterTxt());
        exporterTimeTxt->setup( "time", this->comm_ );

        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterNewtonIterations = Teuchos::rcp(new ExporterTxt());
        exporterNewtonIterations->setup( "newtonIterations" + suffix, this->comm_ );
        
        exporterIterations = Teuchos::rcp(new ExporterTxt());
        exporterIterations->setup( "linearIterations" + suffix, this->comm_ );

        exporterCornerValue = Teuchos::rcp(new ExporterTxt());
        exporterCornerValue->setup( "cornerValue" + suffix, this->comm_ );

        MapConstPtr_Type map = problemTime_->getDomain(0)->getMapUnique();
        vec2D_dbl_ptr_Type points = problemTime_->getDomain(0)->getPointsUnique();

        valueCorner = map->getLocalElement(idExport);
        if(valueCorner != -1)
            cout <<" Value corner local ID " << valueCorner  << " with node values " << points->at(valueCorner).at(0) << " " << points->at(valueCorner).at(1) << " " << points->at(valueCorner).at(2)  << endl;


    }
    if (printExtraData) {

        vec_dbl_Type v(3,-9999.);
        this->problemTime_->getValuesOfInterest(v);
        vec_dbl_Type vGathered(this->comm_->getSize());
        Teuchos::gatherAll<int,double>( *this->comm_, 1, &v[0], vGathered.size(), &vGathered[0] );
        int targetRank=0;
        while (vGathered[targetRank] < 0){
            targetRank++;
            TEUCHOS_TEST_FOR_EXCEPTION( targetRank == vGathered.size(), std::runtime_error, "No targetRank for export of displacements was found!" );
        }
        
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterDisplXTxt->setup( "displ_x" + suffix, this->comm_ , targetRank);
        exporterDisplYTxt->setup( "displ_y" + suffix, this->comm_ , targetRank);
        
    }

    vec2D_dbl_Type timeParametersVec(0,vec_dbl_Type(2));
    
    int numSegments = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments; i++){

        double startTime = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec.push_back(segment);
    }
    double loadStepSize = parameterList_->sublist("Parameter").get("Load Step Size",1.);

    if(numSegments > 0 ){
        TEUCHOS_TEST_FOR_EXCEPTION( loadStepSize != timeParametersVec[0][1], std::runtime_error, "Load Step Size and First Time Interval Size appear different" );
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( loadStepSize != timeSteppingTool_->dt_, std::runtime_error, "Load Step Size and dt appear different" );
    }
    double dt;
    for(int i=0; i<numSegments ; i++){
        if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0]){
            dt=timeParametersVec[i][1];
            timeSteppingTool_->dt_ = dt;
        }

    }
    
    // Notwendige Parameter
    int sizeSCI = timeStepDef_.size();

    // ACHTUNG
    int sizeChem = 1; //  c
    int sizeStructure = 1; // d_s

    dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    // ######################
    // Chem: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffChem(sizeChem);
    SmallMatrix<double> problemCoeffChem(sizeChem);
    double coeffSourceTermChem = 0.0;

    /*for (int i=0; i<sizeChem; i++) {
        for (int j=0; j<sizeChem; j++) {
            if (timeStepDef_[i][j]>0 && i==j) {
                massCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
            }
            else{
                massCoeffChem[i][j] = 0.0;
            }
        }
    }
    for (int i=0; i<sizeChem; i++) {
        for (int j=0; j<sizeChem; j++){
            if (timeStepDef_[i][j]>0){
                problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);
            }
            else{
                problemCoeffChem[i][j] = 1.;
            }
        }
    }*/
    
    massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt; // 3/(2\Delta t)
    problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1); // 1
    coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1); // 1

    // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    SmallMatrix<double> massCoeffStructure(sizeStructure);
    SmallMatrix<double> problemCoeffStructure(sizeStructure);
    double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    // Koeffizient vor der Massematrix
    /*for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
            // i == j, da vektorwertige Massematrix blockdiagonal ist
            if(timeStepDef_[i + sizeChem][j + sizeChem] > 0  && i == j) // Weil: (c, d_s,...) und timeStepDef_ von FSI
            {
               // Vorfaktor der Massematrix in der LHS
                massCoeffStructure[i][j] = 1.0/(dt*dt*beta);
            }
            else
            {
                massCoeffStructure[i][j] = 0.;
            }
        }
    }*/
    
    // Die anderen beiden Koeffizienten
    /*for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            if(timeStepDef_[i + sizeChem][j + sizeChem] > 0 )
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
    }*/

    massCoeffStructure[0][0] = 1./(dt*dt*beta);
    problemCoeffStructure[0][0] =  1.0;
    coeffSourceTermStructure = 1.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!

    bool chemistryExplicit_ =    parameterList_->sublist("Parameter").get("Chemistry Explicit",false);

    // ######################
    // SCI: Mass-, Problem-Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffSCI(sizeSCI);
    SmallMatrix<double> problemCoeffSCI(sizeSCI);
    for (int i = 0; i < sizeStructure; i++)
    {
        for (int j = 0; j < sizeStructure; j++)
        {
            massCoeffSCI[i][j] = massCoeffStructure[i][j];
            problemCoeffSCI[i][j] = problemCoeffStructure[i][j];
        }
    }

    if(!chemistryExplicit_){
        for (int i = 0; i < sizeChem; i++)
        {
            for (int j = 0; j < sizeChem; j++)
            {
                massCoeffSCI[i + sizeStructure][j + sizeStructure] = massCoeffChem[i][j];
                problemCoeffSCI[i + sizeStructure][j + sizeStructure] = problemCoeffChem[i][j];
            }
        }
    }
   
    this->problemTime_->setTimeParameters(massCoeffSCI, problemCoeffSCI);
    
    if (printExtraData) {
        exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
        vec_dbl_Type v(3,0.);
        this->problemTime_->getValuesOfInterest( v );

        exporterDisplXTxt->exportData( v[0] );
        exporterDisplYTxt->exportData( v[1] );
    }

    // ######################
    // Time loop
    // ######################


    double inflowRamp = parameterList_->sublist("Parameter").get("Inflow Ramp",0.01);
    std::string structureModel = parameterList_->sublist("Parameter").get("Structure Model","SCI_NH");
    std::string couplingType = parameterList_->sublist("Parameter").get("Coupling Type","explicit");

    double timeStep = 0;

    while(timeSteppingTool_->continueTimeStepping())
    {
        for(int i=0; i<numSegments ; i++){
            if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0])
                dt=timeParametersVec[i][1];
        }
        timeSteppingTool_->dt_= dt;
        sci->timeSteppingTool_->dt_ = dt;
        if(timeSteppingTool_->currentTime() <= 0. + 1e-12){
            timeSteppingTool_->dt_prev_= dt;        
            sci->timeSteppingTool_->dt_prev_= dt;        
        }
        else{
            timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_;

            this->problemTime_->assemble("UpdateTime"); // Updates to next timestep

            sci->timeSteppingTool_->dt_prev_ = timeSteppingTool_->dt_;

        }

        

        timeSteppingTool_->printInfo();

        massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt; // 3/(2\Delta t)
        problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1); // 1
        coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1); // 1

        if(!chemistryExplicit_){
            if(timeSteppingTool_->currentTime() > 0. + 1.e-12) 
            {
                for (int i = 0; i < sizeChem; i++)
                {
                    for (int j = 0; j < sizeChem; j++){
                        massCoeffSCI[i+sizeStructure][j+sizeStructure] = massCoeffChem[i][j];
                    }
                }
                this->problemTime_->setTimeParameters(massCoeffSCI, problemCoeffSCI);
            }
        }
        problemTime_->updateTime ( timeSteppingTool_->currentTime() );

        //string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");

        // Ist noetig, falls wir extrapolieren, damit wir
        // immer die korrekten previousSolution_ haben.
        // TODO: Vermutlich reicht lediglich (da erstmal nur BDF2):
        // this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        if(nmbBDF<2  && !parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation"))
        {// we need the last two solution for a second order extrapolation.
            if (timeSteppingTool_->currentTime() != 0.0)
            {
                this->problemTime_->updateSolutionMultiPreviousStep(2);
            }
            else
            {
                this->problemTime_->updateSolutionMultiPreviousStep(1);
            }
        }
        else
        {
          this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        }
        

        // Alte Gitterbewegung mit der Geometrieloesung ueberschreiben.
        // -- we can keep this as expicit update for the reaction-diffusion displacement
        this->problemTime_->assemble("UpdateMeshDisplacement");   
       
        if(chemistryExplicit_)
        {
            this->problemTime_->assemble("UpdateChemInTime");

            this->problemTime_->assemble("MoveMesh");

            this->problemTime_->assemble("SolveChemistryProblem");
            
        }

        // ######################
        // Struktur Zeitsystem
        // ######################
        // In jedem Zeitschritt die RHS der Struktur holen.
        // Die Massematrix wird in FSI jedoch nur fuer t = 0 berechnet, da Referenzkonfiguration
        // in der Struktur.
        {
            
            // Hier wird auch direkt ein Update der Loesung bei der Struktur gemacht.
            // Aehnlich zu "UpdateFluidInTime".
            
            if(timeSteppingTool_->currentTime() == 0.0)
            {
                // We extract the underlying FSI problem
                MatrixPtr_Type massmatrix;
                sci->setSolidMassmatrix( massmatrix );
                this->problemTime_->systemMass_->addBlock( massmatrix, 0, 0 );
            }
            // this should be done automatically rhs will not be used here
            //  this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 2 );
            this->problemTime_->assemble("ComputeSolidRHSInTime");
        }

        // ######################
        // Chem Zeitsystem
        // ######################
        // Fluid-Loesung aktualisieren fuer die naechste(n) BDF2-Zeitintegration(en)
        // in diesem Zeitschritt.
        {
                //Do we need this, if BDF for FSI is used correctly? We still need it to save the mass matrices
            if(couplingType=="explicit")// || structureModel=="SCI_sophisticated")
                this->problemTime_->assemble("UpdateChemInTime");
        }        

        MatrixPtr_Type massmatrix;
        sci->setChemMassmatrix( massmatrix );
        //massmatrix->print();
        //if(couplingType=="explicit" ) //|| structureModel=="SCI_sophisticated")
        //    this->problemTime_->systemMass_->addBlock( massmatrix, 1, 1);

        // RHS nach BDF2
        //if(couplingType=="explicit" )//|| structureModel=="SCI_sophisticated")
        //    this->problemTime_->assemble( "ComputeChemRHSInTime" ); // hier ist massmatrix nicht relevant
        //this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 0 );

      
        // ######################
        // System loesen
        // ######################
        // Use BDF1 Parameters for first system
        if (timeSteppingTool_->currentTime() == 0.) {
            if(!chemistryExplicit_){
                for (int i = 0; i < sizeChem; i++)
                {
                    for (int j = 0; j < sizeChem; j++){
                        if (massCoeffSCI[i+sizeStructure][j+sizeStructure] != 0.){
                            massCoeffSCI[i+sizeStructure][j+sizeStructure] = 1./dt ;
                        }
                    }
                }
            
            this->problemTime_->setTimeParameters(massCoeffSCI, problemCoeffSCI);
            }
        }

        
        double time = timeSteppingTool_->currentTime() +  timeSteppingTool_->dt_;
        problemTime_->updateTime ( time );
        
        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));

        if(this->comm_->getRank()==0)
            cout << " ----- Nonlinear System Info|| Number of rows:" << problemTime_->getSystem()->size() << " || number of rhs:" << problemTime_->getRhs()->size() << " || solution: " << problemTime_->getSolution()->size() << endl; 
        //massCoeffSCI.print();
        //problemCoeffSCI.print();
        if("SCI_Linear" != parameterList_->sublist("Parameter").get("Structure Model","SCI_Linear"))
            nlSolver.solve(*this->problemTime_, time, its);
        else{
            problemTime_->combineSystems();
            problemTime_->setBoundaries(time); 
            (*its)[0]=problemTime_->solve();
        }


        if(!chemistryExplicit_)
        {
            if (timeSteppingTool_->currentTime() <= dt+1.e-10) 
            {
                for (int i = 0; i < sizeChem; i++)
                {
                    for (int j = 0; j < sizeChem; j++){
                        massCoeffSCI[i+sizeStructure][j+sizeStructure] = massCoeffChem[i][j];
                    }
                }
                this->problemTime_->setTimeParameters(massCoeffSCI, problemCoeffSCI);
            }
        }
        
        //this->problemTime_->computeValuesOfInterestAndExport();

        timeSteppingTool_->advanceTime(false);//output info);

        // Should be some place else
        //if(couplingType=="explicit" )
        //    this->problemTime_->assemble("UpdateCoupling");

        if (printData) {
            exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
            exporterIterations->exportData( (*its)[0] );
            exporterNewtonIterations->exportData( (*its)[1] );

            vec_dbl_Type d_s(0);
            double norm=0.;

            string name = parameterList_->sublist("General").get("Physic","Structure");             
            if(valueCorner != -1){

                if(name == "Structure"){
                    for(int i=0; i< problemTime_->dimension_ ; i++)
                        d_s.push_back(problemTime_->getSolution()->getBlock(0)->getDataNonConst(0)[problemTime_->dimension_*valueCorner+i]);
                    
                    for(int i=0; i< problemTime_->dimension_ ; i++)
                        norm += pow(d_s[i],2);
                    norm = sqrt(norm);
                }
                else{
                    norm = problemTime_->getSolution()->getBlock(1)->getDataNonConst(0)[valueCorner];
                }


            }

			Teuchos::reduceAll<int, double> ( *this->comm_, Teuchos::REDUCE_MAX, norm , Teuchos::outArg (norm));

            exporterCornerValue->exportData(norm);


        }
        if (printExtraData) {
            vec_dbl_Type v(3,-9999.);
            this->problemTime_->getValuesOfInterest(v);
            
            exporterDisplXTxt->exportData( v[0] );
            exporterDisplYTxt->exportData( v[1] );
        }
        if (print)
        {
            exportTimestep();

        }
        if (printStress){
            timeStep = timeStep + 1.;
            bool heartbeat= false;
            double heartbeatStart1 = parameterList_->sublist("Parameter").get("Heart Beat Start 1",0.) ;
            double heartbeatStart2 = parameterList_->sublist("Parameter").get("Heart Beat Start 2",0.) ;
            double heartbeatEnd1 = parameterList_->sublist("Parameter").get("Heart Beat End 1",0.) ;
            double heartbeatEnd2 = parameterList_->sublist("Parameter").get("Heart Beat End 2",0.) ;
            double modValue = parameterList_->sublist("Parameter").get("Every X Second",1.) ;
            double modValueHeartBeat = parameterList_->sublist("Parameter").get("Every X Second Heartbeat",0.1) ;

            double time = timeSteppingTool_->currentTime();

            if(time >= heartbeatStart1 && time <= heartbeatEnd1)
                heartbeat=true;
            else if(time >= heartbeatStart2 && time <= heartbeatEnd2)
                heartbeat=true;

            if(heartbeat)
                modValue= modValueHeartBeat; // smaller post Processing steps in heart beat phases



            if(fabs(remainder(timeStep,modValue)) < 0. + 1.e-8 ){
                BlockMultiVectorPtr_Type stressVecTmp= sci->getPostProcessingData();
                stressVec = stressVecTmp;
                this->exportPostprocess(stressVec,problemTime_->getDomain(0),sci->getPostprocessingNames()); 
            }

        }

    }

    comm_->barrier();
    if (printExtraData) {
        exporterTimeTxt->closeExporter();
        exporterIterations->closeExporter();
        exporterNewtonIterations->closeExporter();
    }
    if (printExtraData) {
        exporterDisplXTxt->closeExporter();
        exporterDisplYTxt->closeExporter();        
    }
    if (print)
    {
        closeExporter();
        closeExporterPostprocess();
    }
}


template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceWithLoadStepping()
{
    // problemCoeff vor A (= komplettes steady-System)
    // massCoeff vor M (= Massematrix)
    // coeffSourceTerm vor f (= rechte Seite der DGL)
    cout << " Advance with Loadstepping " << endl;
    NonLinElasProblemPtr_Type nonLinElas = Teuchos::rcp_dynamic_cast<NonLinElasProblem_Type>( this->problemTime_->getUnderlyingProblem() );

    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    bool printStress = parameterList_->sublist("General").get("Export Stress",false);
    bool printData = parameterList_->sublist("General").get("Export Data",true);
    bool printExtraData = parameterList_->sublist("General").get("Export Extra Data",false);

    BlockMultiVectorPtr_Type stressVec;
    if (print)
    {
        exportTimestep();
    }

    vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations
    ExporterTxtPtr_Type exporterTimeTxt;
    ExporterTxtPtr_Type exporterDisplXTxt;
    ExporterTxtPtr_Type exporterDisplYTxt;
    ExporterTxtPtr_Type exporterIterations;
    ExporterTxtPtr_Type exporterNewtonIterations;
    ExporterTxtPtr_Type exporterCornerValue;

    GO idExport = parameterList_->sublist("General").get("ExportID",0);
    LO valueCorner =-1;

    
    if (printData) {
        exporterTimeTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplXTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplYTxt = Teuchos::rcp(new ExporterTxt());
        exporterTimeTxt->setup( "time", this->comm_ );

        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterNewtonIterations = Teuchos::rcp(new ExporterTxt());
        exporterNewtonIterations->setup( "newtonIterations" + suffix, this->comm_ );
        
        exporterIterations = Teuchos::rcp(new ExporterTxt());
        exporterIterations->setup( "linearIterations" + suffix, this->comm_ );

        exporterCornerValue = Teuchos::rcp(new ExporterTxt());
        exporterCornerValue->setup( "cornerValue" + suffix, this->comm_ );

        MapConstPtr_Type map = problemTime_->getDomain(0)->getMapUnique();
        vec2D_dbl_ptr_Type points = problemTime_->getDomain(0)->getPointsUnique();

        valueCorner = map->getLocalElement(idExport);
        if(valueCorner != -1)
            cout <<" Value corner local ID " << valueCorner  << " with node values " << points->at(valueCorner).at(0) << " " << points->at(valueCorner).at(1) << " " << points->at(valueCorner).at(2)  << endl;


    }
   
    vec2D_dbl_Type timeParametersVec(0,vec_dbl_Type(2));
    
    int numSegments = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments; i++){

        double startTime = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec.push_back(segment);
    }
    double loadStepSize = parameterList_->sublist("Parameter").get("Load Step Size",1.);

    if(numSegments > 0 ){
        TEUCHOS_TEST_FOR_EXCEPTION( loadStepSize != timeParametersVec[0][1], std::runtime_error, "Load Step Size and First Time Interval Size appear different" );
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( loadStepSize != timeSteppingTool_->dt_, std::runtime_error, "Load Step Size and dt appear different" );
    }
    double dt;
    for(int i=0; i<numSegments ; i++){
        if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0]){
            dt=timeParametersVec[i][1];
            timeSteppingTool_->dt_ = dt;
        }

    }
    
   
    int sizeStructure = 1; // d_s

    dt = timeSteppingTool_->get_dt();
  // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    SmallMatrix<double> massCoeff(sizeStructure);
    SmallMatrix<double> problemCoeff(sizeStructure);
    double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    massCoeff[0][0] = 1.;
    problemCoeff[0][0] =  1.0;
    coeffSourceTermStructure = 0.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!

    this->problemTime_->setTimeParameters(massCoeff, problemCoeff);
    
    if (printExtraData) {
        exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
        vec_dbl_Type v(3,0.);
        this->problemTime_->getValuesOfInterest( v );

        exporterDisplXTxt->exportData( v[0] );
        exporterDisplYTxt->exportData( v[1] );
    }

    // ######################
    // Time loop
    // ######################
    while(timeSteppingTool_->continueTimeStepping())
    {
        for(int i=0; i<numSegments ; i++){
            if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0])
                dt=timeParametersVec[i][1];
        }
        timeSteppingTool_->dt_= dt;
        nonLinElas->timeSteppingTool_->dt_ = dt;
        if(timeSteppingTool_->currentTime() <= 0. + 1e-12){
            timeSteppingTool_->dt_prev_= dt;        
            nonLinElas->timeSteppingTool_->dt_prev_= dt;        
        }
        else{
            timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_;

            this->problemTime_->assemble("UpdateTime"); // Updates to next timestep

            nonLinElas->timeSteppingTool_->dt_prev_ = timeSteppingTool_->dt_;

        }

        

        timeSteppingTool_->printInfo();

       
        problemTime_->updateTime ( timeSteppingTool_->currentTime() );
           
        {
            
            // Hier wird auch direkt ein Update der Loesung bei der Struktur gemacht.
            // Aehnlich zu "UpdateFluidInTime".
            
            /*if(timeSteppingTool_->currentTime() == 0.0)
            {
                // We extract the underlying FSI problem
                MatrixPtr_Type massmatrix;
                sci->setSolidMassmatrix( massmatrix );
                this->problemTime_->systemMass_->addBlock( massmatrix, 0, 0 );
            }*/
            // this should be done automatically rhs will not be used here
            //  this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 2 );
            this->problemTime_->assemble("ComputeSolidRHSLoad");
        }

        double time = timeSteppingTool_->currentTime() +  timeSteppingTool_->dt_;
        problemTime_->updateTime ( time );
        
        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));
        //massCoeffSCI.print();
        //problemCoeffSCI.print();
       // if("SCI_Linear" != parameterList_->sublist("Parameter").get("Structure Model","SCI_Linear"))
            nlSolver.solve(*this->problemTime_, time, its);
        /*else{
            problemTime_->combineSystems();
            problemTime_->setBoundaries(time); 
            (*its)[0]=problemTime_->solve();
            }
        */  
        //this->problemTime_->computeValuesOfInterestAndExport();

        timeSteppingTool_->advanceTime(false);//output info);

       

        if (printData) {
            exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
            exporterIterations->exportData( (*its)[0] );
            exporterNewtonIterations->exportData( (*its)[1] );

            vec_dbl_Type d_s(0);
            double norm=0.;

            string name = parameterList_->sublist("General").get("Physic","Structure");             
            if(valueCorner != -1){

                if(name == "Structure"){
                    for(int i=0; i< problemTime_->dimension_ ; i++)
                        d_s.push_back(problemTime_->getSolution()->getBlock(0)->getDataNonConst(0)[problemTime_->dimension_*valueCorner+i]);
                    
                    for(int i=0; i< problemTime_->dimension_ ; i++)
                        norm += pow(d_s[i],2);
                    norm = sqrt(norm);
                }
                else{
                    norm = problemTime_->getSolution()->getBlock(1)->getDataNonConst(0)[valueCorner];
                }


            }

			Teuchos::reduceAll<int, double> ( *this->comm_, Teuchos::REDUCE_MAX, norm , Teuchos::outArg (norm));

            exporterCornerValue->exportData(norm);


        }
        if (printExtraData) {
            vec_dbl_Type v(3,-9999.);
            this->problemTime_->getValuesOfInterest(v);
            
            exporterDisplXTxt->exportData( v[0] );
            exporterDisplYTxt->exportData( v[1] );
        }
        if (print)
        {
            exportTimestep();

        }
      

    }

    comm_->barrier();
    if (printExtraData) {
        exporterTimeTxt->closeExporter();
        exporterIterations->closeExporter();
        exporterNewtonIterations->closeExporter();
    }
    if (printExtraData) {
        exporterDisplXTxt->closeExporter();
        exporterDisplYTxt->closeExporter();        
    }
    if (print)
    {
        closeExporter();
        closeExporterPostprocess();
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeFSCI()
{
    // problemCoeff vor A (= komplettes steady-System)
    // massCoeff vor M (= Massematrix)
    // coeffSourceTerm vor f (= rechte Seite der DGL)
    //cout << " ###### advanceInTimeFSCI #########" << endl;
  
    
    FSCIProblemPtr_Type fsci = Teuchos::rcp_dynamic_cast<FSCIProblem_Type>( this->problemTime_->getUnderlyingProblem() );
    
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    bool printData = parameterList_->sublist("General").get("Export Data",false);
    bool printStress = parameterList_->sublist("General").get("Export Stress",false);
    bool printFlowRate = parameterList_->sublist("General").get("Export Flow Rate",true);
    bool printExtraData = parameterList_->sublist("General").get("Export Extra Data",false);

    BlockMultiVectorPtr_Type stressVec;
  
    if (print)
    {
        exportTimestep();
    }


    vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations
    ExporterTxtPtr_Type exporterTimeTxt;
    ExporterTxtPtr_Type exporterDisplXTxt;
    ExporterTxtPtr_Type exporterDisplYTxt;
    ExporterTxtPtr_Type exporterIterations;
    ExporterTxtPtr_Type exporterNewtonIterations;
    ExporterTxtPtr_Type exporterFlowRateInlet;
    ExporterTxtPtr_Type exporterFlowRateOutlet;
    ExporterTxtPtr_Type exporterAreaInlet;
    ExporterTxtPtr_Type exporterAreaOutlet;
    ExporterTxtPtr_Type exporterPressureOutlet;


    
    if (printData) {
        exporterTimeTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplXTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplYTxt = Teuchos::rcp(new ExporterTxt());
        exporterTimeTxt->setup( "time", this->comm_ );

        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterNewtonIterations = Teuchos::rcp(new ExporterTxt());
        exporterNewtonIterations->setup( "newtonIterations" + suffix, this->comm_ );
        
        exporterIterations = Teuchos::rcp(new ExporterTxt());
        exporterIterations->setup( "linearIterations" + suffix, this->comm_ );

    }
    if (printFlowRate) {
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");

        exporterFlowRateInlet = Teuchos::rcp(new ExporterTxt());
        exporterFlowRateInlet->setup( "flowRateInlet" + suffix, this->comm_ );

        exporterFlowRateOutlet = Teuchos::rcp(new ExporterTxt());
        exporterFlowRateOutlet->setup( "flowRateOutlet" + suffix, this->comm_ );

        exporterPressureOutlet = Teuchos::rcp(new ExporterTxt());
        exporterPressureOutlet->setup( "pressureOutlet" + suffix, this->comm_ );

        exporterAreaInlet = Teuchos::rcp(new ExporterTxt());
        exporterAreaInlet->setup( "areaInlet" + suffix, this->comm_ );

        exporterAreaOutlet = Teuchos::rcp(new ExporterTxt());
        exporterAreaOutlet->setup( "areaOutlet" + suffix, this->comm_ );


    }
    if (printExtraData) {

        vec_dbl_Type v(3,-9999.);
        this->problemTime_->getValuesOfInterest(v);
        vec_dbl_Type vGathered(this->comm_->getSize());
        Teuchos::gatherAll<int,double>( *this->comm_, 1, &v[0], vGathered.size(), &vGathered[0] );
        int targetRank=0;
        while (vGathered[targetRank] < 0){
            targetRank++;
            TEUCHOS_TEST_FOR_EXCEPTION( targetRank == vGathered.size(), std::runtime_error, "No targetRank for export of displacements was found!" );
        }
        
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterDisplXTxt->setup( "displ_x" + suffix, this->comm_ , targetRank);
        exporterDisplYTxt->setup( "displ_y" + suffix, this->comm_ , targetRank);
        
    }
    vec2D_dbl_Type timeParametersVec(0,vec_dbl_Type(2));
    
    int numSegments = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments; i++){

        double startTime = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec.push_back(segment);
    }
    double dt;
    for(int i=0; i<numSegments ; i++){
        if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0]){
            dt=timeParametersVec[i][1];
            timeSteppingTool_->dt_ = dt;
        }

    }

    // Notwendige Parameter
    bool geometryExplicit = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);
    bool chemistryExplicit_ =    parameterList_->sublist("Parameter").get("Chemistry Explicit",false);

    //std::string couplingType = parameterList_->sublist("Parameter").get("Coupling Type","explicit");

    int sizeFSI = timeStepDef_.size();

    // ACHTUNG
    int sizeFluid = 2; // u_f  + p
    int sizeStructure = 1; // d_s
    int sizeChem =1; // c

    dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

   
    if (printExtraData) {
        exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
        vec_dbl_Type v(3,0.);
        this->problemTime_->getValuesOfInterest( v );

        exporterDisplXTxt->exportData( v[0] );
        exporterDisplYTxt->exportData( v[1] );
    }

    NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));


//    {
//        // Initialize mass matrix for fluid problem
//        MatrixPtr_Type massmatrix;
//        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Fix reassembly of mass matrix for FSI.");
//        MatrixPtr_Type massmatrix;
//        fsi->setSolidMassmatrix( massmatrix );
////        this->problemTime_->assemble( massmatrix, "GetFluidMassmatrix" );
//    }
    // ######################
    // Time loop
    // ######################
#ifdef FEDD_TIMER
    TimeMonitor_Type solveTM(*solveProblemTimer_);
#endif
    while(timeSteppingTool_->continueTimeStepping())
    {
        //cout << " ###### Timeloop #########" << endl;
        for(int i=0; i<numSegments ; i++){
            if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0])
                dt=timeParametersVec[i][1];
        }
        timeSteppingTool_->dt_= dt;
        fsci->timeSteppingTool_->dt_ = dt;
        fsci->problemSCI_->timeSteppingTool_->dt_ = dt;
        if(timeSteppingTool_->currentTime() <= 0. + 1e-12){
            timeSteppingTool_->dt_prev_= dt;        
            fsci->timeSteppingTool_->dt_prev_= dt;  
            fsci->problemSCI_->timeSteppingTool_->dt_prev_ = dt; 
            //fsci->problemFluid_->timeSteppingTool_->dt_prev_= dt; 
        }
        else{
            timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_;

            this->problemTime_->assemble("UpdateTime"); // Updates to next timestep previously:problemTime_->updateTime ( timeSteppingTool_->currentTime() );

            fsci->timeSteppingTool_->dt_prev_ = timeSteppingTool_->dt_;
            fsci->problemSCI_->timeSteppingTool_->dt_prev_ = timeSteppingTool_->dt_;
            //fsci->problemFluid_->timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_; 

        }

        timeSteppingTool_->printInfo();

        string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");

        problemTime_->updateTime ( timeSteppingTool_->currentTime() );


        // ######################
        // Fluid: Mass-, Problem, SourceTerm Koeffizienten
        // ######################
        SmallMatrix<double> massCoeffFluid(sizeFluid);
        SmallMatrix<double> problemCoeffFluid(sizeFluid);
        double coeffSourceTermFluid = 0.0;

        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if (timeStepDef_[i][j]>0 && i==j) {
                    massCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffFluid[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if (timeStepDef_[i][j]>0){
                    problemCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermFluid = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffFluid[i][j] = 1.;
                }
            }
        }


        // ######################
        // Struktur: Mass-, Problem, SourceTerm Koeffizienten
        // ######################
        // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
        SmallMatrix<double> massCoeffStructure(sizeStructure);
        SmallMatrix<double> problemCoeffStructure(sizeStructure);
        double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

        // Koeffizient vor der Massematrix
        for(int i = 0; i < sizeStructure; i++)
        {
            for(int j = 0; j < sizeStructure; j++)
            {
                // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
                // i == j, da vektorwertige Massematrix blockdiagonal ist
                if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0  && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSI
                {
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
                if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0 )
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

            // ######################
        // Chem: Mass-, Problem, SourceTerm Koeffizienten
        // ######################
        SmallMatrix<double> massCoeffChem(sizeChem);
        SmallMatrix<double> problemCoeffChem(sizeChem);
        double coeffSourceTermChem = 0.0;

        /*for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if (timeStepDef_[i][j]>0 && i==j) {
                    massCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffChem[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if (timeStepDef_[i][j]>0){
                    problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffChem[i][j] = 1.;
                }
            }
        }*/
        
        massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt; // 3/(2\Delta t)
        problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1); // 1
        coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1); // 1

        // ######################
        // FSI: Mass-, Problem-Koeffizienten
        // ######################
        SmallMatrix<double> massCoeffFSI(sizeFSI);
        SmallMatrix<double> problemCoeffFSI(sizeFSI);
        for (int i = 0; i < sizeFluid; i++)
        {
            for (int j = 0; j < sizeFluid; j++)
            {
                massCoeffFSI[i][j] = massCoeffFluid[i][j];
                problemCoeffFSI[i][j] = problemCoeffFluid[i][j];
            }
        }

        for (int i = 0; i < sizeStructure; i++)
        {
            for (int j = 0; j < sizeStructure; j++)
            {
                massCoeffFSI[i + sizeFluid][j + sizeFluid] = massCoeffStructure[i][j];
                problemCoeffFSI[i + sizeFluid][j + sizeFluid] = problemCoeffStructure[i][j];
            }
        }

        if(!chemistryExplicit_){
            /*
               massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt; // 3/(2\Delta t)
                    problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1); // 1
                    coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1); // 1
             */
            massCoeffFSI[4][4] = massCoeffChem[0][0];
            problemCoeffFSI[4][4] = problemCoeffChem[0][0];
            problemCoeffFSI[2][4] = 1.; // SCI Coupling 1
            problemCoeffFSI[4][2] = 1.; // SCI Coupling 2
        }
        // Setze noch Einsen an die Stellen, wo Eintraege (Kopplungsbloecke) vorhanden sind.
        problemCoeffFSI[0][3] = 1.0; // C1_T
        problemCoeffFSI[2][3] = 1.0; // C3_T
        problemCoeffFSI[3][0] = 1.0; // C1
        problemCoeffFSI[3][2] = 1.0; // C2
        
        /*if(!geometryExplicit)
        {
            problemCoeffFSI[4][2] = 1.0; // C4
            problemCoeffFSI[4][4] = 1.0; // H (Geometrie)
            string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");
            if(linearization == "Newton" || linearization == "NOX")
            {
                problemCoeffFSI[0][4] = 1.0; // Shape-Derivatives Velocity
                problemCoeffFSI[1][4] = 1.0; // Shape-Derivatives Div-Nebenbedingung
            }
        }*/

        this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        // Ist noetig, falls wir extrapolieren, damit wir
        // immer die korrekten previousSolution_ haben.
        // TODO: Vermutlich reicht lediglich (da erstmal nur BDF2):
        // this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        if(nmbBDF<2 && !parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation"))
        {// we need the last two solution for a second order extrapolation.
            if (timeSteppingTool_->currentTime() != 0.0)
            {
                this->problemTime_->updateSolutionMultiPreviousStep(2);
            }
            else
            {
                this->problemTime_->updateSolutionMultiPreviousStep(1);
            }
        }
        else
        {
            this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        }
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleAddInterfaceRHSTimer_);
#endif
            // Den Block C2*d_s^n in der RHS im Interface-Block setzen.
            this->problemTime_->assemble("AddInterfaceBlockRHS");
        }
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleUpdateMeshDisplacementTimer_);
#endif
            // Alte Gitterbewegung mit der Geometrieloesung ueberschreiben.
            this->problemTime_->assemble("UpdateMeshDisplacement");
        }
        // Das Geometry-Problem separat loesen, falls GE.
        if(geometryExplicit)
        {
            {
#ifdef FEDD_DETAIL_TIMER
                TimeMonitor_Type reassmbleTM(*reassmbleSolveGeometryTimer_);
#endif
                this->problemTime_->assemble("SolveGeometryProblem");
            }
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleMoveMeshTimer_);
#endif
            this->problemTime_->assemble("MoveMesh");
        }

        if(chemistryExplicit_)
        {

            fsci->problemSCI_->assemble("UpdateChemInTime");

            fsci->problemSCI_->assemble("MoveMesh");
            
            this->problemTime_->assemble("SolveChemistryProblem");
        
        }


        // ######################
        // Struktur Zeitsystem
        // ######################
        // In jedem Zeitschritt die RHS der Struktur holen.
        // Die Massematrix wird in FSI jedoch nur fuer t = 0 berechnet, da Referenzkonfiguration
        // in der Struktur.
        {
                
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleSolidMassAndRHSTimer_);
#endif
            // Hier wird auch direkt ein Update der Loesung bei der Struktur gemacht.
            // Aehnlich zu "UpdateFluidInTime".
            
            if(timeSteppingTool_->currentTime() == 0.0)
            {
                // We extract the underlying FSI problem
                MatrixPtr_Type massmatrix;
                fsci->setSolidMassmatrix( massmatrix );
                this->problemTime_->systemMass_->addBlock( massmatrix, 2, 2 );
            }
            // this should be done automatically rhs will not be used here
//            this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 2 );
            this->problemTime_->assemble("ComputeSolidRHSInTime");
        }

        if(geometryExplicit) //  && linearization != "Extrapolation"
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleForTimeTimer_);
#endif
            this->problemTime_->assemble("ForTime");
        }

        // ######################
        // Fluid Zeitsystem
        // ######################
        // Fluid-Loesung aktualisieren fuer die naechste(n) BDF2-Zeitintegration(en)
        // in diesem Zeitschritt.
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleUpdateFluidInTimeTimer_);
#endif
            //Do we need this, if BDF for FSI is used correctly? We still need it to save the mass matrices
            this->problemTime_->assemble("UpdateFluidInTime");

            this->problemTime_->assemble("ComputePressureRHSInTime");

        }
        // Aktuelle Massematrix auf dem Gitter fuer BDF2-Integration und
        // fuer das FSI-System (bei GI wird die Massematrix weiterhin in TimeProblem.reAssemble() assembliert).
        // In der ersten nichtlinearen Iteration wird bei GI also die Massematrix zweimal assembliert.
        // Massematrix fuer FSI holen und fuer timeProblemFluid setzen (fuer BDF2)
        MatrixPtr_Type massmatrix;
        fsci->setFluidMassmatrix( massmatrix );
        this->problemTime_->systemMass_->addBlock( massmatrix, 0, 0 );
        

        // RHS nach BDF2
        this->problemTime_->assemble( "ComputeFluidRHSInTime" ); // hier ist massmatrix nicht relevant
        // this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 0 );

        // Aktuelle Massematrix auf dem Gitter fuer BDF2-Integration und
        // fuer das FSI-System (bei GI wird die Massematrix weiterhin in TimeProblem.reAssemble() assembliert).
        // In der ersten nichtlinearen Iteration wird bei GI also die Massematrix zweimal assembliert.
        // Massematrix fuer FSI holen und fuer timeProblemFluid setzen (fuer BDF2)
        MatrixPtr_Type massmatrixC;
        fsci->setChemMassmatrix( massmatrixC );
       // if(couplingType=="explicit" ) //|| structureModel=="SCI_sophisticated")
       //     this->problemTime_->systemMass_->addBlock( massmatrixC, 4, 4);

        // RHS nach BDF2
        //if(couplingType=="explicit" ) //|| structureModel=="SCI_sophisticated")
        //    this->problemTime_->assemble( "ComputeChemRHSInTime" ); // hier ist mas smatrix nicht relevant
        //this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 0 );


        // ######################
        // System loesen
        // ######################
        // Use BDF1 Parameters for first system
        if (timeSteppingTool_->currentTime() == 0.) {
            for (int i = 0; i < sizeFluid; i++)
            {
                for (int j = 0; j < sizeFluid; j++){
                    if (massCoeffFSI[i][j] != 0.)
                        massCoeffFSI[i][j] = 1./dt ;
                }
            }
            this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        }
        
        double time = timeSteppingTool_->currentTime() +  timeSteppingTool_->dt_;
        problemTime_->updateTime ( time );      

        nlSolver.solve(*this->problemTime_, time, its);
        
        if (timeSteppingTool_->currentTime() <= dt+1.e-10) {
            for (int i = 0; i < sizeFluid; i++)
            {
                for (int j = 0; j < sizeFluid; j++){
                    massCoeffFSI[i][j] = massCoeffFluid[i][j];
                }
            }
            this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        }
        
        this->problemTime_->computeValuesOfInterestAndExport();

        timeSteppingTool_->advanceTime(false/*output info*/);
        //this->problemTime_->assemble("UpdateTime"); // Zeit in FSI inkrementieren

        if (printData) {
            exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
            exporterIterations->exportData( (*its)[0] );
            exporterNewtonIterations->exportData( (*its)[1] );
        }
        if(printFlowRate){
            FE<SC,LO,GO,NO> fe;
		    fe.addFE(problemTime_->getDomain(0));
            double flowRateInlet;
            double flowRateOutlet;

            int flagInlet = this->parameterList_->sublist("General").get("Flag Inlet Fluid", 4);
            int flagOutlet = this->parameterList_->sublist("General").get("Flag Outlet Fluid", 5);

            MultiVectorPtr_Type u_rep = Teuchos::rcp(new MultiVector_Type ( problemTime_->getDomain(0)->getMapVecFieldRepeated() ) );   
    	    u_rep->importFromVector(problemTime_->getSolution()->getBlock(0),false,"Insert");
            fe.assemblyFlowRate(problemTime_->getDomain(0)->getDimension(), flowRateInlet, problemTime_->getDomain(0)->getFEType() , problemTime_->getDomain(0)->getDimension(), flagInlet , u_rep);
            fe.assemblyFlowRate(problemTime_->getDomain(0)->getDimension(), flowRateOutlet, problemTime_->getDomain(0)->getFEType() , problemTime_->getDomain(0)->getDimension(), flagOutlet , u_rep);

            exporterFlowRateInlet->exportData( timeSteppingTool_->currentTime() , flowRateInlet );
            exporterFlowRateOutlet->exportData(  timeSteppingTool_->currentTime() , flowRateOutlet );

            exporterPressureOutlet->exportData(  timeSteppingTool_->currentTime() ,fsci->getPressureOutlet() );

            double areaInlet=0.;
            fe.assemblyArea(problemTime_->getDomain(0)->getDimension(), areaInlet, flagInlet);

            double areaOutlet=0.;
            fe.assemblyArea(problemTime_->getDomain(0)->getDimension(), areaOutlet, flagOutlet);

            exporterAreaInlet->exportData(  timeSteppingTool_->currentTime() ,areaInlet);
            exporterAreaOutlet->exportData(  timeSteppingTool_->currentTime() ,areaOutlet );

        }
        if (printExtraData) {
            vec_dbl_Type v(3,-9999.);
            this->problemTime_->getValuesOfInterest(v);
            
            exporterDisplXTxt->exportData( v[0] );
            exporterDisplYTxt->exportData( v[1] );
        }
        if (printStress){
            bool heartbeat= false;
            double heartbeatStart1 = parameterList_->sublist("Parameter").get("Heart Beat Start 1",0.) ;
            double heartbeatStart2 = parameterList_->sublist("Parameter").get("Heart Beat Start 2",0.) ;
            double heartbeatEnd1 = parameterList_->sublist("Parameter").get("Heart Beat End 1",0.) ;
            double heartbeatEnd2 = parameterList_->sublist("Parameter").get("Heart Beat End 2",0.) ;

            double time = timeSteppingTool_->currentTime();

            if(time >= heartbeatStart1 && time <= heartbeatEnd1)
                heartbeat=true;
            else if(time >= heartbeatStart2 && time <= heartbeatEnd2)
                heartbeat=true;

            double modValue = 10.;

            if(heartbeat)
                modValue= 0.1; // smaller post Processing steps in heart beat phases

            if(fmod(timeSteppingTool_->currentTime(),modValue) < 0. + 5.e-2 ){
                BlockMultiVectorPtr_Type stressVecTmp= fsci->problemSCI_->getPostProcessingData();
                stressVec = stressVecTmp;
                this->exportPostprocess(stressVec,problemTime_->getDomain(2),fsci->problemSCI_->getPostprocessingNames());
            }
        }
        if (print)
        {
            exportTimestep();
        }

    }

    comm_->barrier();
    if (printExtraData) {
        exporterTimeTxt->closeExporter();
        exporterIterations->closeExporter();
        exporterNewtonIterations->closeExporter();
    }
    if (printExtraData) {
        exporterDisplXTxt->closeExporter();
        exporterDisplYTxt->closeExporter();        
    }
    if(printFlowRate){
        exporterFlowRateOutlet->closeExporter();
        exporterFlowRateInlet->closeExporter();
        exporterAreaOutlet->closeExporter();
        exporterAreaInlet->closeExporter();
        exporterPressureOutlet->closeExporter();

    }
    if (print)
    {
        closeExporter();
    }
}



template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeFSI()
{
    // problemCoeff vor A (= komplettes steady-System)
    // massCoeff vor M (= Massematrix)
    // coeffSourceTerm vor f (= rechte Seite der DGL)
    
    
    FSIProblemPtr_Type fsi = Teuchos::rcp_dynamic_cast<FSIProblem_Type>( this->problemTime_->getUnderlyingProblem() );
    
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    bool printData = parameterList_->sublist("General").get("Export Data",false);
    bool printFlowRate = parameterList_->sublist("General").get("Export Flow Rate",true);
    bool printExtraData = parameterList_->sublist("General").get("Export Extra Data",false);
        
    if (print)
    {
        exportTimestep();
    }


    vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations
    ExporterTxtPtr_Type exporterTimeTxt;
    ExporterTxtPtr_Type exporterDisplXTxt;
    ExporterTxtPtr_Type exporterDisplYTxt;
    ExporterTxtPtr_Type exporterIterations;
    ExporterTxtPtr_Type exporterNewtonIterations;
    ExporterTxtPtr_Type exporterFlowRateInlet;
    ExporterTxtPtr_Type exporterFlowRateOutlet;
    ExporterTxtPtr_Type exporterAreaInlet;
    ExporterTxtPtr_Type exporterAreaOutlet;
    ExporterTxtPtr_Type exporterPressureOutlet;

    
    if (printData) {
        exporterTimeTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplXTxt = Teuchos::rcp(new ExporterTxt());
        exporterDisplYTxt = Teuchos::rcp(new ExporterTxt());
        exporterTimeTxt->setup( "time", this->comm_ );

        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterNewtonIterations = Teuchos::rcp(new ExporterTxt());
        exporterNewtonIterations->setup( "newtonIterations" + suffix, this->comm_ );
        
        exporterIterations = Teuchos::rcp(new ExporterTxt());
        exporterIterations->setup( "linearIterations" + suffix, this->comm_ );

    }
    if (printFlowRate) {
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");

        exporterFlowRateInlet = Teuchos::rcp(new ExporterTxt());
        exporterFlowRateInlet->setup( "flowRateInlet" + suffix, this->comm_ );

        exporterFlowRateOutlet = Teuchos::rcp(new ExporterTxt());
        exporterFlowRateOutlet->setup( "flowRateOutlet" + suffix, this->comm_ );

        exporterPressureOutlet = Teuchos::rcp(new ExporterTxt());
        exporterPressureOutlet->setup( "pressureOutlet" + suffix, this->comm_ );

        exporterAreaInlet = Teuchos::rcp(new ExporterTxt());
        exporterAreaInlet->setup( "areaInlet" + suffix, this->comm_ );

        exporterAreaOutlet = Teuchos::rcp(new ExporterTxt());
        exporterAreaOutlet->setup( "areaOutlet" + suffix, this->comm_ );


    }

    if (printExtraData) {

        vec_dbl_Type v(3,-9999.);
        this->problemTime_->getValuesOfInterest(v);
        vec_dbl_Type vGathered(this->comm_->getSize());
        Teuchos::gatherAll<int,double>( *this->comm_, 1, &v[0], vGathered.size(), &vGathered[0] );
        int targetRank=0;
        while (vGathered[targetRank] < 0){
            targetRank++;
            TEUCHOS_TEST_FOR_EXCEPTION( targetRank == vGathered.size(), std::runtime_error, "No targetRank for export of displacements was found!" );
        }
        
        std::string suffix = parameterList_->sublist("General").get("Export Suffix","");
        
        exporterDisplXTxt->setup( "displ_x" + suffix, this->comm_ , targetRank);
        exporterDisplYTxt->setup( "displ_y" + suffix, this->comm_ , targetRank);
        
    }
    
    vec2D_dbl_Type timeParametersVec(0,vec_dbl_Type(2));
    
    int numSegments = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments; i++){

        double startTime = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = parameterList_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec.push_back(segment);
    }
    double dt;
    for(int i=0; i<numSegments ; i++){
        if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0]){
            dt=timeParametersVec[i][1];
            timeSteppingTool_->dt_ = dt;
        }

    }
    // Notwendige Parameter
    bool geometryExplicit = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);

    int sizeFSI = timeStepDef_.size();

    // ACHTUNG
    int sizeFluid = 2; // u_f  + p
    int sizeStructure = 1; // d_s

    dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    // ######################
    // Fluid: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffFluid(sizeFluid);
    SmallMatrix<double> problemCoeffFluid(sizeFluid);
    double coeffSourceTermFluid = 0.0;

    for (int i=0; i<sizeFluid; i++) {
        for (int j=0; j<sizeFluid; j++) {
            if (timeStepDef_[i][j]>0 && i==j) {
                massCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
            }
            else{
                massCoeffFluid[i][j] = 0.0;
            }
        }
    }
    for (int i=0; i<sizeFluid; i++) {
        for (int j=0; j<sizeFluid; j++){
            if (timeStepDef_[i][j]>0){
                problemCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(1);
                coeffSourceTermFluid = timeSteppingTool_->getInformationBDF(1);
            }
            else{
                problemCoeffFluid[i][j] = 1.;
            }
        }
    }


    // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    SmallMatrix<double> massCoeffStructure(sizeStructure);
    SmallMatrix<double> problemCoeffStructure(sizeStructure);
    double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    // Koeffizient vor der Massematrix
    for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
            // i == j, da vektorwertige Massematrix blockdiagonal ist
            if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0  && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSI
            {
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
            if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0 )
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


    // ######################
    // FSI: Mass-, Problem-Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffFSI(sizeFSI);
    SmallMatrix<double> problemCoeffFSI(sizeFSI);
    for (int i = 0; i < sizeFluid; i++)
    {
        for (int j = 0; j < sizeFluid; j++)
        {
            massCoeffFSI[i][j] = massCoeffFluid[i][j];
            problemCoeffFSI[i][j] = problemCoeffFluid[i][j];
        }
    }

    for (int i = 0; i < sizeStructure; i++)
    {
        for (int j = 0; j < sizeStructure; j++)
        {
            massCoeffFSI[i + sizeFluid][j + sizeFluid] = massCoeffStructure[i][j];
            problemCoeffFSI[i + sizeFluid][j + sizeFluid] = problemCoeffStructure[i][j];
        }
    }

    // Setze noch Einsen an die Stellen, wo Eintraege (Kopplungsbloecke) vorhanden sind.
    problemCoeffFSI[0][3] = 1.0; // C1_T
    problemCoeffFSI[2][3] = 1.0; // C3_T
    problemCoeffFSI[3][0] = 1.0; // C1
    problemCoeffFSI[3][2] = 1.0; // C2
    if(!geometryExplicit)
    {
        problemCoeffFSI[4][2] = 1.0; // C4
        problemCoeffFSI[4][4] = 1.0; // H (Geometrie)
        string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");
        if(linearization == "Newton" || linearization == "NOX")
        {
            problemCoeffFSI[0][4] = 1.0; // Shape-Derivatives Velocity
            problemCoeffFSI[1][4] = 1.0; // Shape-Derivatives Div-Nebenbedingung
        }
    }

    this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
    
    if (printExtraData) {
        exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
        vec_dbl_Type v(3,0.);
        this->problemTime_->getValuesOfInterest( v );

        exporterDisplXTxt->exportData( v[0] );
        exporterDisplYTxt->exportData( v[1] );
    }

//    {
//        // Initialize mass matrix for fluid problem
//        MatrixPtr_Type massmatrix;
//        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Fix reassembly of mass matrix for FSI.");
//        MatrixPtr_Type massmatrix;
//        fsi->setSolidMassmatrix( massmatrix );
////        this->problemTime_->assemble( massmatrix, "GetFluidMassmatrix" );
//    }
    // ######################
    // Time loop
    // ######################
#ifdef FEDD_TIMER
    TimeMonitor_Type solveTM(*solveProblemTimer_);
#endif
    while(timeSteppingTool_->continueTimeStepping())
    {
         //cout << " ###### Timeloop #########" << endl;
        for(int i=0; i<numSegments ; i++){
            if(timeSteppingTool_->currentTime()+1.0e-12 > timeParametersVec[i][0])
                dt=timeParametersVec[i][1];
        }
        timeSteppingTool_->dt_= dt;
        fsi->timeSteppingTool_->dt_ = dt;
        if(timeSteppingTool_->currentTime() <= 0. + 1e-12){
            timeSteppingTool_->dt_prev_= dt;        
            fsi->timeSteppingTool_->dt_prev_= dt;  
            //fsci->problemFluid_->timeSteppingTool_->dt_prev_= dt; 
        }
        else{
            timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_;

            this->problemTime_->assemble("UpdateTime"); // Updates to next timestep previously:problemTime_->updateTime ( timeSteppingTool_->currentTime() );

            fsi->timeSteppingTool_->dt_prev_ = timeSteppingTool_->dt_;
            //fsci->problemFluid_->timeSteppingTool_->dt_prev_= timeSteppingTool_->dt_; 

        }

        timeSteppingTool_->printInfo();

        string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");

        problemTime_->updateTime ( timeSteppingTool_->currentTime() );


        // ######################
        // Fluid: Mass-, Problem, SourceTerm Koeffizienten
        // ######################
        SmallMatrix<double> massCoeffFluid(sizeFluid);
        SmallMatrix<double> problemCoeffFluid(sizeFluid);
        double coeffSourceTermFluid = 0.0;

        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if (timeStepDef_[i][j]>0 && i==j) {
                    massCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffFluid[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if (timeStepDef_[i][j]>0){
                    problemCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermFluid = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffFluid[i][j] = 1.;
                }
            }
        }


        // ######################
        // Struktur: Mass-, Problem, SourceTerm Koeffizienten
        // ######################
        // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
        SmallMatrix<double> massCoeffStructure(sizeStructure);
        SmallMatrix<double> problemCoeffStructure(sizeStructure);
        double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

        // Koeffizient vor der Massematrix
        for(int i = 0; i < sizeStructure; i++)
        {
            for(int j = 0; j < sizeStructure; j++)
            {
                // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
                // i == j, da vektorwertige Massematrix blockdiagonal ist
                if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0  && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSI
                {
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
                if(timeStepDef_[i + sizeFluid][j + sizeFluid] > 0 )
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


        // ######################
        // FSI: Mass-, Problem-Koeffizienten
        // ######################
        SmallMatrix<double> massCoeffFSI(sizeFSI);
        SmallMatrix<double> problemCoeffFSI(sizeFSI);
        for (int i = 0; i < sizeFluid; i++)
        {
            for (int j = 0; j < sizeFluid; j++)
            {
                massCoeffFSI[i][j] = massCoeffFluid[i][j];
                problemCoeffFSI[i][j] = problemCoeffFluid[i][j];
            }
        }

        for (int i = 0; i < sizeStructure; i++)
        {
            for (int j = 0; j < sizeStructure; j++)
            {
                massCoeffFSI[i + sizeFluid][j + sizeFluid] = massCoeffStructure[i][j];
                problemCoeffFSI[i + sizeFluid][j + sizeFluid] = problemCoeffStructure[i][j];
            }
        }

        // Setze noch Einsen an die Stellen, wo Eintraege (Kopplungsbloecke) vorhanden sind.
        problemCoeffFSI[0][3] = 1.0; // C1_T
        problemCoeffFSI[2][3] = 1.0; // C3_T
        problemCoeffFSI[3][0] = 1.0; // C1
        problemCoeffFSI[3][2] = 1.0; // C2
        if(!geometryExplicit)
        {
            problemCoeffFSI[4][2] = 1.0; // C4
            problemCoeffFSI[4][4] = 1.0; // H (Geometrie)
            string linearization = this->parameterList_->sublist("General").get("Linearization","Extrapolation");
            if(linearization == "Newton" || linearization == "NOX")
            {
                problemCoeffFSI[0][4] = 1.0; // Shape-Derivatives Velocity
                problemCoeffFSI[1][4] = 1.0; // Shape-Derivatives Div-Nebenbedingung
            }
        }

        this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        // Ist noetig, falls wir extrapolieren, damit wir
        // immer die korrekten previousSolution_ haben.
        // TODO: Vermutlich reicht lediglich (da erstmal nur BDF2):
        // this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        if(nmbBDF<2 && !parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation"))
        {// we need the last two solution for a second order extrapolation.
            if (timeSteppingTool_->currentTime() != 0.0)
            {
                this->problemTime_->updateSolutionMultiPreviousStep(2);
            }
            else
            {
                this->problemTime_->updateSolutionMultiPreviousStep(1);
            }
        }
        else
        {
            this->problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        }
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleAddInterfaceRHSTimer_);
#endif
            // Den Block C2*d_s^n in der RHS im Interface-Block setzen.
            this->problemTime_->assemble("AddInterfaceBlockRHS");
        }
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleUpdateMeshDisplacementTimer_);
#endif
            // Alte Gitterbewegung mit der Geometrieloesung ueberschreiben.
            this->problemTime_->assemble("UpdateMeshDisplacement");
        }
        // Das Geometry-Problem separat loesen, falls GE.
        if(geometryExplicit)
        {
            {
#ifdef FEDD_DETAIL_TIMER
                TimeMonitor_Type reassmbleTM(*reassmbleSolveGeometryTimer_);
#endif
                this->problemTime_->assemble("SolveGeometryProblem");
            }
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleMoveMeshTimer_);
#endif
            this->problemTime_->assemble("MoveMesh");
        }


        // ######################
        // Struktur Zeitsystem
        // ######################
        // In jedem Zeitschritt die RHS der Struktur holen.
        // Die Massematrix wird in FSI jedoch nur fuer t = 0 berechnet, da Referenzkonfiguration
        // in der Struktur.
        {
            
            
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleSolidMassAndRHSTimer_);
#endif
            // Hier wird auch direkt ein Update der Loesung bei der Struktur gemacht.
            // Aehnlich zu "UpdateFluidInTime".
            
            if(timeSteppingTool_->currentTime() == 0.0)
            {
                // We extract the underlying FSI problem
                MatrixPtr_Type massmatrix;
                fsi->setSolidMassmatrix( massmatrix );
                this->problemTime_->systemMass_->addBlock( massmatrix, 2, 2 );
            }
            // this should be done automatically rhs will not be used here
//            this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 2 );
            this->problemTime_->assemble("ComputeSolidRHSInTime");
        }


        if(geometryExplicit) //  && linearization != "Extrapolation"
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleForTimeTimer_);
#endif
            this->problemTime_->assemble("ForTime");
        }

        // ######################
        // Fluid Zeitsystem
        // ######################
        // Fluid-Loesung aktualisieren fuer die naechste(n) BDF2-Zeitintegration(en)
        // in diesem Zeitschritt.
        {
#ifdef FEDD_DETAIL_TIMER
            TimeMonitor_Type reassmbleTM(*reassmbleUpdateFluidInTimeTimer_);
#endif
            //Do we need this, if BDF for FSI is used correctly? We still need it to save the mass matrices
            this->problemTime_->assemble("UpdateFluidInTime");

            this->problemTime_->assemble("ComputePressureRHSInTime");

        }
        // Aktuelle Massematrix auf dem Gitter fuer BDF2-Integration und
        // fuer das FSI-System (bei GI wird die Massematrix weiterhin in TimeProblem.reAssemble() assembliert).
        // In der ersten nichtlinearen Iteration wird bei GI also die Massematrix zweimal assembliert.
        // Massematrix fuer FSI holen und fuer timeProblemFluid setzen (fuer BDF2)
        MatrixPtr_Type massmatrix;
        fsi->setFluidMassmatrix( massmatrix );
        this->problemTime_->systemMass_->addBlock( massmatrix, 0, 0 );
        

        // RHS nach BDF2
        this->problemTime_->assemble( "ComputeFluidRHSInTime" ); // hier ist massmatrix nicht relevant
//        this->problemTime_->getRhs()->addBlock( Teuchos::rcp_const_cast<MultiVector_Type>(rhs->getBlock(0)), 0 );


        // ######################
        // System loesen
        // ######################
        // Use BDF1 Parameters for first system
        if (timeSteppingTool_->currentTime() == 0.) {
            for (int i = 0; i < sizeFluid; i++)
            {
                for (int j = 0; j < sizeFluid; j++){
                    if (massCoeffFSI[i][j] != 0.)
                        massCoeffFSI[i][j] = 1./dt ;
                }
            }
            this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        }
            
        double time = timeSteppingTool_->currentTime() +  timeSteppingTool_->dt_;
        problemTime_->updateTime ( time );            
        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));

        nlSolver.solve(*this->problemTime_, time, its);
        
        if (timeSteppingTool_->currentTime() <= dt+1.e-10) {
            for (int i = 0; i < sizeFluid; i++)
            {
                for (int j = 0; j < sizeFluid; j++){
                    massCoeffFSI[i][j] = massCoeffFluid[i][j];
                }
            }
            this->problemTime_->setTimeParameters(massCoeffFSI, problemCoeffFSI);
        }
        
        this->problemTime_->computeValuesOfInterestAndExport();

        //timeSteppingTool_->advanceTime(true/*output info*/);
        //this->problemTime_->assemble("UpdateTime"); // Zeit in FSI inkrementieren
        timeSteppingTool_->advanceTime(false/*output info*/);

        if (printData) {
            exporterTimeTxt->exportData( timeSteppingTool_->currentTime() );
            exporterIterations->exportData( (*its)[0] );
            exporterNewtonIterations->exportData( (*its)[1] );
        }
        if(printFlowRate){
            FE<SC,LO,GO,NO> fe;
		    fe.addFE(problemTime_->getDomain(0));
            double flowRateInlet;
            double flowRateOutlet;

            int flagInlet = this->parameterList_->sublist("General").get("Flag Inlet Fluid", 4);
            int flagOutlet = this->parameterList_->sublist("General").get("Flag Outlet Fluid", 5);

            MultiVectorPtr_Type u_rep = Teuchos::rcp(new MultiVector_Type ( problemTime_->getDomain(0)->getMapVecFieldRepeated() ) );   
    	    u_rep->importFromVector(problemTime_->getSolution()->getBlock(0),false,"Insert");
            fe.assemblyFlowRate(problemTime_->getDomain(0)->getDimension(), flowRateInlet, problemTime_->getDomain(0)->getFEType() , problemTime_->getDomain(0)->getDimension(), flagInlet , u_rep);
            fe.assemblyFlowRate(problemTime_->getDomain(0)->getDimension(), flowRateOutlet, problemTime_->getDomain(0)->getFEType() , problemTime_->getDomain(0)->getDimension(), flagOutlet , u_rep);

            exporterFlowRateInlet->exportData(  timeSteppingTool_->currentTime() , flowRateInlet );
            exporterFlowRateOutlet->exportData(  timeSteppingTool_->currentTime() ,flowRateOutlet );

            exporterPressureOutlet->exportData(  timeSteppingTool_->currentTime() , fsi->getPressureOutlet() );

            double areaInlet=0.;
            fe.assemblyArea(problemTime_->getDomain(0)->getDimension(), areaInlet, flagInlet);

            double areaOutlet=0.;
            fe.assemblyArea(problemTime_->getDomain(0)->getDimension(), areaOutlet, flagOutlet);

            exporterAreaInlet->exportData( timeSteppingTool_->currentTime() , areaInlet);
            exporterAreaOutlet->exportData(  timeSteppingTool_->currentTime() ,areaOutlet );

        }
        if (printExtraData) {
            vec_dbl_Type v(3,-9999.);
            this->problemTime_->getValuesOfInterest(v);
            
            exporterDisplXTxt->exportData( v[0] );
            exporterDisplYTxt->exportData( v[1] );
        }
        if (print)
        {
            exportTimestep();
        }

    }

    comm_->barrier();
    if (printExtraData) {
        exporterTimeTxt->closeExporter();
        exporterIterations->closeExporter();
        exporterNewtonIterations->closeExporter();
    }
    if (printExtraData) {
        exporterDisplXTxt->closeExporter();
        exporterDisplYTxt->closeExporter();        
    }
    if(printFlowRate){
        exporterFlowRateOutlet->closeExporter();
        exporterFlowRateInlet->closeExporter();
        exporterAreaOutlet->closeExporter();
        exporterAreaInlet->closeExporter();
        exporterPressureOutlet->closeExporter();
    }
    if (print)
    {
        closeExporter();
    }
}


/* UEBERARBEITEN!!!!!!!!!!!!!!!!! */
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeLinearMultistep(){

    //TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "advanceInTimeLinearMultistep rework.");
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print) {
        exportTimestep();
    }
    int size = timeStepDef_.size();
    double dt = timeSteppingTool_->get_dt();
    int nmbBDF = timeSteppingTool_->getBDFNumber();
    vec_dbl_Type coeffPrevSteps(nmbBDF);
    for (int i=0; i<coeffPrevSteps.size(); i++) {
       coeffPrevSteps.at(i) = timeSteppingTool_->getInformationBDF(i+2) / dt;
    }

    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    double coeffSourceTerm = 0.;

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (timeStepDef_[i][j]==1 && i==j) {
                massCoeff[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
            }
            else{
                massCoeff[i][j] = 0.;
            }
        }
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++){
            if (timeStepDef_[i][j]==1){
                problemCoeff[i][j] =  timeSteppingTool_->getInformationBDF(1);
                coeffSourceTerm = timeSteppingTool_->getInformationBDF(1); // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!
            }
            else{
                problemCoeff[i][j] = 1.;
            }
        }
    }
    problemTime_->setTimeParameters(massCoeff, problemCoeff);

	//#########
    //time loop
    //#########
    while (timeSteppingTool_->continueTimeStepping()) {

        problemTime_->updateSolutionMultiPreviousStep(nmbBDF);

        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateMultistepRhs(coeffPrevSteps,nmbBDF);/*apply mass matrix to u_t*/
        if (problemTime_->hasSourceTerm()) {
            problemTime_->assembleSourceTerm(time);
            addSourceTermToRHS(coeffSourceTerm);
        }

        problemTime_->combineSystems();
        // Uebergabeparameter fuer BC noch hinzunehmen!
        problemTime_->setBoundaries(time);
        problemTime_->solve();

        timeSteppingTool_->advanceTime(true/*output info*/);

        if (print) {
            exportTimestep();
        }

    }
    comm_->barrier();
    if (print) {
        closeExporter();
    }

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeNonLinearMultistep(){

    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print) {
        exportTimestep();
    }

    int size = timeStepDef_.size();
    double dt = timeSteppingTool_->get_dt();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    vec_dbl_Type coeffPrevSteps(nmbBDF);
    for (int i=0; i<coeffPrevSteps.size(); i++) {
        coeffPrevSteps.at(i) = timeSteppingTool_->getInformationBDF(i+2) / dt;
    }

    SmallMatrix<double> massCoeff(size);
    SmallMatrix<double> problemCoeff(size);
    double coeffSourceTerm = 0.;

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (timeStepDef_[i][j]>0 && i==j) {
                massCoeff[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
            }
            else{
                massCoeff[i][j] = 0.;
            }
        }
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++){
            if (timeStepDef_[i][j]>0){
                problemCoeff[i][j] = timeSteppingTool_->getInformationBDF(1);
                coeffSourceTerm = timeSteppingTool_->getInformationBDF(1);
            }
            else{
                problemCoeff[i][j] = 1.;
            }
        }
    }
    problemTime_->setTimeParameters(massCoeff, problemCoeff);
    //#########
    //time loop
    //#########
    while (timeSteppingTool_->continueTimeStepping()) {

        // For the first time step we use BDF1
        if (timeSteppingTool_->currentTime()==0.) {
            SmallMatrix<double> tmpmassCoeff(size);
            SmallMatrix<double> tmpproblemCoeff(size);
            for (int i=0; i<size; i++) {
                for (int j=0; j<size; j++) {
                    if (timeStepDef_[i][j]>0 && i==j) {
                        tmpmassCoeff[i][j] = 1. / dt;
                    }
                    else{
                        tmpmassCoeff[i][j] = 0.;
                    }
                }
            }
            for (int i=0; i<size; i++) {
                for (int j=0; j<size; j++){
                    if (timeStepDef_[i][j]>0){
                        tmpproblemCoeff[i][j] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat
                    }
                    else{
                        tmpproblemCoeff[i][j] = 1.;
                    }
                }
            }
            problemTime_->setTimeParameters(tmpmassCoeff, tmpproblemCoeff);
        }
        if(nmbBDF<2 && !parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation")) {
            if (timeSteppingTool_->currentTime()!=0.){
                problemTime_->updateSolutionMultiPreviousStep(2);
            }
            else{
                problemTime_->updateSolutionMultiPreviousStep(1);
            }
        }
        else{
            problemTime_->updateSolutionMultiPreviousStep(nmbBDF);
        }
        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateTime ( time );
        
        if (timeSteppingTool_->currentTime()==0.) {
            vec_dbl_Type tmpcoeffPrevSteps(1, 1. / dt);
            problemTime_->updateMultistepRhs(tmpcoeffPrevSteps,1);/*apply (mass matrix / dt) to u_t*/
        }
        else{
            problemTime_->updateMultistepRhs(coeffPrevSteps,nmbBDF);/*apply (mass matrix / dt) to u_t*/
        }
        if (problemTime_->hasSourceTerm()) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Check sourceterm.");
//            problemTime_->AssembleSourceTerm(time);
//            if (timeSteppingTool_->CurrentTime()==0.) {
//               AddSourceTermToRHS(1.);
//            }
//            else{
//                AddSourceTermToRHS(coeffSourceTerm); //ACHTUNG
//            }
        }

        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","FixedPoint"));
        nlSolver.solve(*problemTime_,time);

        // After the first time step we can use the desired BDF Parameters
        if (timeSteppingTool_->currentTime()==0.) {
            problemTime_->setTimeParameters(massCoeff, problemCoeff);
        }

        timeSteppingTool_->advanceTime(true/*output info*/);

        if (print) {
            exportTimestep();
        }
    }

    if (print) {
        closeExporter();
    }

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeLinearExternal(){
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print) {
        exportTimestep();
    }
    
    double dt = timeSteppingTool_->get_dt();
    int counter=0;
    //#########
    //time loop
    //#########
    while (timeSteppingTool_->continueTimeStepping()) {
        
        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateTime ( time );

        problemTime_->assemble( "Assemble" );
        
        if (problemTime_->hasSourceTerm()){
            problemTime_->assembleSourceTerm(time); // assemble source term and add to rhs
            problemTime_->addToRhs( problemTime_->getSourceTerm() );
        }
        
        problemTime_->setBoundaries( time );
        
        BlockMultiVectorPtr_Type tmpSol = Teuchos::rcp(new BlockMultiVector_Type(problemTime_->getSolution()) );
        
        //problemTime_->getSystemCombined()->print();
        //problemTime_->getRhs()->print();
        
        problemTime_->solve();
        
        //AceGen systems and solutions are always for a Newton method, i.e., we solve for an update even in the linear case.
                    // ?? Alpha a + beta this?
        problemTime_->getSolution()->update( 1., tmpSol, -1. );
        
        problemTime_->assemble( "SetSolutionNewton" ); // Newton solution is actually the solution at t+1 and time solution is still at time t
    
        problemTime_->assemble( "OnlyUpdate" ); // updates AceGEN history
        
        problemTime_->assemble( "SetSolutionInTime" ); // Here, we update the time solution with the Newton solution, which basically means that we can go to the next timestep
        
        timeSteppingTool_->advanceTime(true/*output info*/);
        counter++;
        if (print)
            exportTimestep();
    }
    
    if (print) {
        closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::advanceInTimeNonLinearExternal(){
        
    bool print = parameterList_->sublist("General").get("ParaViewExport",false);
    if (print) {
        exportTimestep();
    }
    
    double dt = timeSteppingTool_->get_dt();
    int counter=0;
    //#########
    //time loop
    //#########
    while (timeSteppingTool_->continueTimeStepping()) {
        
        double time = timeSteppingTool_->currentTime() + dt;
        problemTime_->updateTime ( time );
            
        if (problemTime_->hasSourceTerm()){
            problemTime_->assembleSourceTerm(time); // assemble source term; added to residual vector in calculateNonLinResidualVec()
        }
        NonLinearSolver<SC, LO, GO, NO> nlSolver(parameterList_->sublist("General").get("Linearization","Newton"));
        nlSolver.solve(*problemTime_,time);
        problemTime_->assemble( "OnlyUpdate" ); // updates AceGEN history
        
        problemTime_->assemble( "SetSolutionInTime" ); // Here, we update the time solution with the Newton solution, which basically means that we can go to the next timestep
        
        timeSteppingTool_->advanceTime(true/*output info*/);
        counter++;
        if (print) {
            exportTimestep();
        }
    }
    
    if (print) {
        closeExporter();
    }
}
    
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::addRhsDAE(SmallMatrix<double> coeff, BlockMatrixPtr_Type bMat, BlockMultiVectorPtr_Type vec){

    BlockMultiVectorPtr_Type mv = Teuchos::rcp( new BlockMultiVector_Type( vec ) );    
    bMat->apply( *vec, *mv , coeff );
    problemTime_->getRhs()->update( 1., *mv, 1. );
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::addRhsDAE(SmallMatrix<double> coeff, BlockMultiVectorPtr_Type vec){
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Is this correct? addRhsDAE!");
    for (int i=0; i<vec->size(); i++)
        problemTime_->getRhs()->getBlockNonConst(i)->update( coeff[i][i], vec->getBlock(i), 1. );
}
    
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::addSourceTermToRHS(double coeff)
{
   BlockMultiVectorPtr_Type tmpMV = problemTime_->getSourceTerm();
    
   problemTime_->getRhs()->update(coeff, *tmpMV, 1.);
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::exportTimestep(){

    if (!boolExporterSetup_) {
        setupExporter();
    }
    if (verbose_) {
        cout << "-- Exporting..."<< flush;
    }
    for (int i=0; i<exporter_vector_.size(); i++) {
        
        exporter_vector_[i]->save( timeSteppingTool_->currentTime() );

    }

    if (verbose_) {
        cout << "done! --"<< endl;
    }

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::exportPostprocess(BlockMultiVectorPtr_Type postProcessVec, DomainConstPtr_Type domain, vec_string_Type exportNames){

    if (!boolExporterSetupPostprocess_) {
        setupExporter(postProcessVec,domain,exportNames);
    }
    if (verbose_) {
        cout << "-- Exporting Postprocessing Data..."<< flush;
    }
    for (int i=0; i<postProcessVec->size(); i++) {

        MultiVectorConstPtr_Type exportMV= postProcessVec->getBlock(i);
        exporter_vector_postprocess_[0]->updateVariables(exportMV, exportNames[i]);

    }    
    exporter_vector_postprocess_[0]->save( timeSteppingTool_->currentTime() );

    if (verbose_) {
        cout << "done! --"<< endl;
    }

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setupExporter(BlockMultiVectorPtr_Type postProcessVec, DomainConstPtr_Type domain, vec_string_Type exportNames){

    if(!boolExporterSetupPostprocess_){

        ExporterPtr_Type exporterPtr(new Exporter_Type());
        MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>(domain->getMesh());

        exporterPtr->setup("PostprocessingData", meshNonConst, domain->getFEType(), this->parameterList_);
        for(int i=0; i< postProcessVec->size(); i++){
            
            MultiVectorConstPtr_Type exportVector = postProcessVec->getBlock(i);

            std::string suffix = exportNames[i];
            
            UN dofsPerNode = 1;
            
            if (dofsPerNode == 1)
                exporterPtr->addVariable( exportVector, suffix , "Scalar", dofsPerNode, domain->getMapUnique() );
            else
                exporterPtr->addVariable( exportVector, suffix, "Vector", dofsPerNode, domain->getMapUnique() );

            
        }
        exporter_vector_postprocess_.push_back(exporterPtr);
                
        boolExporterSetupPostprocess_ = true;
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::exportTimestep(BlockMultiVectorPtr_Type& solShort){

    if (!boolExporterSetup_) {
        setupExporter(solShort);
    }
    if (verbose_) {
        cout << "-- Exporting..."<< flush;
    }
    for (int i=0; i<exporter_vector_.size(); i++) {
        
        exporter_vector_[i]->save( timeSteppingTool_->currentTime() );
    }

    if (verbose_) {
        cout << "done! --"<< endl;
    }

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::closeExporter(){

    if (boolExporterSetup_) {
        for (int i=0; i<exporter_vector_.size(); i++)
            exporter_vector_[i]->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::closeExporterPostprocess(){

    if (boolExporterSetupPostprocess_) {
        for (int i=0; i<exporter_vector_postprocess_.size(); i++)
            exporter_vector_postprocess_[i]->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setupExporter(){
    for (int i=0; i<timeStepDef_.size(); i++) {
        // \lambda in FSI, koennen wir nicht exportieren, weil keine Elementliste dafuer vorhanden
        bool exportThisBlock  = true;
        //bool chemistryExplicit = this->parameterList_->sublist("Parameter").get("Chemistry Explicit",false);
        if(this->parameterList_->sublist("Parameter").get("FSI",false) == true  )
            exportThisBlock = (i != 3);
       else if(this->parameterList_->sublist("Parameter").get("FSCI",false) == true){
            exportThisBlock = (i != 3);
            //exportThisBlock = (i != 4);
        //if(chemistryExplicit)   
        //   exportThisBlock = false;
       }

        if(exportThisBlock)
        {
            ExporterPtr_Type exporterPtr(new Exporter_Type());
            MultiVectorConstPtr_Type exportVector = problemTime_->getSolution()->getBlock(i);

            DomainConstPtr_Type dom = problemTime_->getDomain(i);

            if(this->parameterList_->sublist("Parameter").get("FSCI",false) == true)
                if(i==4 )
                    dom = problemTime_->getDomain(5);

            int exportEveryXTimesteps = parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string plSuffix = "Suffix variable" + to_string(i+1);
            std::string suffix = parameterList_->sublist("Exporter").get(plSuffix, "" );
            std::string varName = problemTime_->getVariableName(i) + suffix;
            if(this->parameterList_->sublist("Parameter").get("FSCI",false) == true)
                if(i==4 )
                    varName = problemTime_->getVariableName(5) + suffix;

            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>(dom->getMesh());
            exporterPtr->setup(varName, meshNonConst, dom->getFEType(), parameterList_);
            
//            exporterPtr->setup(dom->getDimension(), dom->getNumElementsGlobal(), dom->getElements(), dom->getPointsUnique(), dom->getMapUnique(), dom->getMapRepeated(), dom->getFEType(), varName, exportEveryXTimesteps, comm_ , parameterList_);


            UN dofsPerNode = problemTime_->getDofsPerNode(i);
            if(this->parameterList_->sublist("Parameter").get("FSCI",false) == true)
                if(i==4 )
                    dofsPerNode = problemTime_->getDofsPerNode(5);

            if (dofsPerNode == 1)
                exporterPtr->addVariable( exportVector, varName, "Scalar", dofsPerNode, dom->getMapUnique() );
            else
                exporterPtr->addVariable( exportVector, varName, "Vector", dofsPerNode, dom->getMapUnique() );

            exporter_vector_.push_back(exporterPtr);
            export_solution_vector_.push_back(exportVector);
        }
    }
    boolExporterSetup_ = true;

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setupExporter(BlockMultiVectorPtr_Type& solShort){
    for (int i=0; i<timeStepDef_.size(); i++) {
        // \lambda in FSI, koennen wir nicht exportieren, weil keine Elementliste dafuer vorhanden
        bool exportThisBlock  = true;
        exportThisBlock = !(this->parameterList_->sublist("Parameter").get("FSI",false) && i == 3);
        if(exportThisBlock)
        {
            ExporterPtr_Type exporterPtr(new Exporter_Type());
            MultiVectorConstPtr_Type exportVector = problemTime_->getSolution()->getBlock(i);
            MultiVectorConstPtr_Type exportVectorShort = solShort->getBlock(i);

            DomainConstPtr_Type dom = problemTime_->getDomain(i);

            int exportEveryXTimesteps = parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string plSuffix = "Suffix variable" + to_string(i+1);
            std::string suffix = parameterList_->sublist("Exporter").get(plSuffix, "" );
            std::string varName = problemTime_->getVariableName(i) + suffix;
            std::string varNameShort = problemTime_->getVariableName(i) + "_short_" + suffix;
            
            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>(dom->getMesh());
            exporterPtr->setup(varName, meshNonConst, dom->getFEType(), parameterList_);
            
//            exporterPtr->setup(dom->getDimension(), dom->getNumElementsGlobal(), dom->getElements(), dom->getPointsUnique(), dom->getMapUnique(), dom->getMapRepeated(), dom->getFEType(), varName, exportEveryXTimesteps, comm_ , parameterList_);

            UN dofsPerNode = problemTime_->getDofsPerNode(i);
            
            if (dofsPerNode == 1){
                exporterPtr->addVariable( exportVector, varName, "Scalar", dofsPerNode, dom->getMapUnique() );
                exporterPtr->addVariable( exportVectorShort, varNameShort, "Scalar", dofsPerNode, dom->getMapUnique() );
            }
            else {
                exporterPtr->addVariable( exportVector, varName, "Vector", dofsPerNode, dom->getMapUnique() );
                exporterPtr->addVariable( exportVectorShort, varNameShort, "Vector", dofsPerNode,     dom->getMapUnique() );
            }


            exporter_vector_.push_back(exporterPtr);
            export_solution_vector_.push_back(exportVector);
        }
    }
    boolExporterSetup_ = true;

}


   //ch 21.03.19: Funktion aendern!
template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setupTimeStepping(){

    int size = this->problem_->getSystem()->size();
    
    problemTime_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problem_,comm_));
    
    // Fuer FSI
    if(this->parameterList_->sublist("Parameter").get("FSI",false) || this->parameterList_->sublist("Parameter").get("SCI",false) || this->parameterList_->sublist("Parameter").get("FSCI",false))
    {
        // Beachte: Massematrix ist schon vektorwertig!
        // Reset auf Massesystem von problemTime_ (=FSI), da auf problemTime_ kein assemble() bzw. assembleMassSystem() aufgerufen wird.
        this->problemTime_->systemMass_.reset(new BlockMatrix_Type(this->problemTime_->getSystem()->size()));
        this->problemTime_->systemCombined_.reset( new BlockMatrix_Type( this->problemTime_->getSystem()->size() ));
    }
    else
    {
        size = timeStepDef_.size();
        SmallMatrix<double> massCoeff(size);
        SmallMatrix<double> problemCoeff(size); //hier egal, nur Massenmatrix wichtig. Deswegen auf Null belassen

        for (int i=0; i<size; i++) {
            for (int j=0; j<size; j++) {
                if (timeStepDef_[i][j]>1 && i==j) {
                    massCoeff[i][j] = 1.0;
                }
                else{
                    massCoeff[i][j] = 0.;
                }
            }
        }

        // Setze massCoeff und problemCoeff in TimeProblem
        // ACHTUNG: problemCoeff ist Null. Es wird also nur die Massematrix in systemCombined_ geschrieben.        
        problemTime_->setTimeParameters(massCoeff,problemCoeff);

        // Assembliere Massematrix -> ggf. ReAssemble() falls nichtlinear ->
        // -> das linearProblem (bzw. nonlinearProblem) von oben mit der Massematrix verbinden
        problemTime_->setTimeDef( timeStepDef_ );
        problemTime_->assemble( "MassSystem" );
    }


}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setTimeStep(double dt){

    checkTimeSteppingDef();


}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::setFinalTime(double T){

    checkTimeSteppingDef();

}

template<class SC,class LO,class GO,class NO>
void DAESolverInTime<SC,LO,GO,NO>::checkTimeSteppingDef(){
#ifdef ASSERTS_WARNINGS
    MYASSERT(isTimeSteppingDefined_,"Time stepping not defined!");
#endif

}
}
#endif
