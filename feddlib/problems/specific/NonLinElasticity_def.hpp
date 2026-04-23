#ifndef NonLinElasticity_def_hpp
#define NonLinElasticity_def_hpp

#ifndef NONLINELAS_START
#define NONLINELAS_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Assemble NonLinElasticity:") + std::string(S))));
#endif

#ifndef NONLINELAS_STOP
#define NONLINELAS_STOP(A) A.reset();
#endif

#include "NonLinElasticity_decl.hpp"
#include <sys/resource.h>
#include <iostream>

/*!
 Definition of NonLinElasticity
 
 @brief NonLinElasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::outArg;

namespace FEDD {





template<class SC,class LO,class GO,class NO>
NonLinElasticity<SC,LO,GO,NO>::NonLinElasticity(const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domain->getComm() ),
u_rep_()
{
    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();
    this->addVariable( domain , FEType , "u" , domain->getDimension());
    this->dim_ = this->getDomain(0)->getDimension();
    
    C_ = this->parameterList_->sublist("Parameter").get("C",1.);
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
//    poissonRatio_ = this->parameterList_->sublist("Parameter").get("Poisson Ratio",0.4);
//    mue_ = this->parameterList_->sublist("Parameter").get("Mu",2.0e+6);
//    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstante \lambda
//    E_ = mue_*2.*(1. + poissonRatio_);
//    lambda_ = (poissonRatio_*E_)/((1 + poissonRatio_)*(1 - 2*poissonRatio_));
    loadStepping_ =    !(parameterList->sublist("Timestepping Parameter").get("Class","Singlestep")).compare("Loadstepping");
    externalForce_ =   parameterList->sublist("Parameter").get("External Force",false);
    nonlinearExternalForce_ = parameterList->sublist("Parameter").get("Nonlinear External Force",false);

    timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));

    postProcessingnames_.resize(5);
    postProcessingnames_ = {"vonMisesStress", "SCirc","SAxial","SRadial","W"};

}

template<class SC,class LO,class GO,class NO>
NonLinElasticity<SC,LO,GO,NO>::~NonLinElasticity(){

}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}
    
template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::assemble(std::string type) const{
    
    if (type == ""){
        if (this->verbose_)
            std::cout << "-- Assembly nonlinear elasticity ... " << std::flush;

        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyEmptyMatrix(A);

        this->system_.reset(new BlockMatrix_Type(1));
        this->system_->addBlock( A, 0, 0 );
                
        double density = this->parameterList_->sublist("Parameter").get("Density",1000.);
        std::string sourceType = 	this->parameterList_->sublist("Parameter").get("Source Type","volume");

        this->assembleSourceTerm( 0. );
        if(sourceType == "volume")
            this->sourceTerm_->scale(density);
        
        this->addToRhs( this->sourceTerm_ );

        this->setBoundariesRHS();
                    
        this->solution_->putScalar(0.);
        
        u_rep_ = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);
        
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
        
        this->reAssemble("Newton-Residual");
    }
    else if(type == "UpdateTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateTime)" << '\n';

        updateTime();
        return;
    }
    else
        this->reAssemble(type);
}

// Damit die richtige timeSteppingTool_->currentTime() genommen wird.
template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssemble(std::string type) const {
    std::string material_model = this->parameterList_->sublist("Parameter").get("Material model","Neo-Hooke");
    if (this->verbose_)
        std::cout << "-- Reassembly nonlinear elasticity with material model " << material_model <<" ("<<type <<") ... " << std::flush;
    
    if (type=="Newton-Residual") {
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);
        
        MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
        
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
        this->getComm()->barrier();
        this->getComm()->barrier();

        { 
            NONLINELAS_START(Assembly," Assembling Jacobian and Residual");    
        #ifdef FEDD_HAVE_ACEGENINTERFACE
            bool useInterface = this->parameterList_->sublist("Parameter").get("Use AceGen Interface", true);
            // std::cout << " Checking useInterface parameter for assembly: " << useInterface << std::endl;
            if(this->getFEType(0) =="P2" && useInterface && this->dim_ == 3 && material_model != "Saint-Venant-Kirchhoff"){
                this->system_->addBlock( W, 0, 0 );  
            
                f->putScalar(0.);   
                this->residualVec_->addBlock( f, 0 ); 

                // std::cout << " Checking sci paramters for assembly: " << this->parameterList_->sublist("Parameter").get("SCI",false) << std::endl;
                

                if(this->parameterList_->sublist("Parameter").get("SCI",false) == true || this->parameterList_->sublist("Parameter").get("FSCI",false) == true )
                    this->feFactory_->assemblyAceDeformDiffuBlock(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,concentration_,u_rep_,this->system_,0,0,this->residualVec_,0, this->parameterList_, "Jacobian", true/*call fillComplete*/);
                else 
                    this->feFactory_->assemblyNonLinearElasticity(this->dim_, this->getDomain(0)->getFEType(),2, this->dim_, u_rep_, this->system_, this->residualVec_, this->parameterList_,true);

            }
            else 
                this->feFactory_->assemblyElasticityJacobianAndStressAceFEM(this->dim_, this->getDomain(0)->getFEType(), W, f, u_rep_, this->parameterList_, C_);
        #else
            this->feFactory_->assemblyElasticityJacobianAndStressAceFEM(this->dim_, this->getDomain(0)->getFEType(), W, f, u_rep_, this->parameterList_, C_);
        #endif
            this->getComm()->barrier();
            this->getComm()->barrier();

            NONLINELAS_STOP(Assembly);
        }
        
        MultiVectorPtr_Type fUnique = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
        fUnique->putScalar(0.);
        fUnique->exportFromVector( f, true, "Add" );

        this->residualVec_->addBlock( fUnique, 0 );

        this->system_->addBlock( W, 0, 0 );


        if(loadStepping_) 
            assembleSourceTermLoadstepping();
 
    }
    else if(type=="Newton"){ //we already assemble the new tangent when we calculate the stresses above
        
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;


    // double current_mem = get_mem_usage_mb();
    // std::cout << "Processor " << this->getComm()->getRank() << " ΔMem: " << current_mem - prev_mem << " MB" << std::endl;
    // print_mem_usage();
    if(this->parameterList_->sublist("General").get("Track Memory",false)){
        this->globalMemoryVector_.push_back(this->memoryLogger_->get_global_memory_usage());
        this->memoryLogger_->log(this->newtonStep_, this->timeSteppingTool_->currentTime());
    }
}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{
    
}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){


    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only Newton/NOX implemented for nonlinear material models!");

}


template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::assembleSourceTermLoadstepping(double time) const
{
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    
    
    if(  loadStepping_ == true){
        // The if condition resets the rhs. If we skip it when we attemp loadstepping, the rhs will be updated continously and wrongly increase with each timestep
        this->getRhs()->scale(0.0);
    }

   
    if (this->hasSourceTerm())
    {
        if(externalForce_){

            MultiVectorPtr_Type FERhs = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));
            vec_dbl_Type funcParameter(4,0.);
            funcParameter[0] = timeSteppingTool_->t_;            
            // how can we use different parameters for different blocks here?
            funcParameter[1] =this->parameterList_->sublist("Parameter").get("Volume force",0.00211);

            funcParameter[3] =this->parameterList_->sublist("Parameter").get("Final time force",1.0);
            funcParameter[4] =this->parameterList_->sublist("Parameter").get("dt",0.1);


            if(nonlinearExternalForce_){
                MatrixPtr_Type A( new Matrix_Type (this->system_->getBlock(0,0)));
                //A->print();
                MatrixPtr_Type AKext(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );          
                MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow()*2 ) );          
                MultiVectorPtr_Type Kext_vec;
                this->feFactory_->assemblyNonlinearSurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, u_rep_,Kext, funcParameter, this->rhsFuncVec_[0],this->parameterList_);
               
                A->addMatrix(1.,AKext,0.);
                Kext->addMatrix(1.,AKext,1.);

                AKext->fillComplete(this->getDomain(0)->getMapVecFieldUnique(),this->getDomain(0)->getMapVecFieldUnique());

                this->system_->addBlock(AKext,0,0);

            }
            else            
                this->feFactory_->assemblySurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, u_rep_, funcParameter, this->rhsFuncVec_[0],this->parameterList_);


            this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs, false, "Add" );
        }
        else
            this->assembleSourceTerm( timeSteppingTool_->t_ );
        //this->problemTimeStructure_->getSourceTerm()->scale(density);
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
        
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpPtr= this->sourceTerm_;

        this->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);
        this->rhs_->addBlock( this->getRhs()->getBlockNonConst(0), 0 );

    }
}


template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    this->reAssemble("Newton-Residual");
    std::string material_model = this->parameterList_->sublist("Parameter").get("Material model","Neo-Hooke");
    bool useInterface = this->parameterList_->sublist("Parameter").get("Use AceGen Interface", true);

#ifdef FEDD_HAVE_ACEGENINTERFACE
    if(this->getFEType(0) =="P2" && useInterface && this->dim_ == 3 && material_model!= "Saint-Venant-Kirchhoff"){
        if(this->parameterList_->sublist("Parameter").get("SCI",false) == true || this->parameterList_->sublist("Parameter").get("FSCI",false) == true ){
            this->feFactory_->assemblyAceDeformDiffuBlock(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,concentration_,u_rep_,this->system_,0,0,this->residualVec_,0, this->parameterList_, "Rhs", true/*call fillComplete*/);
        }            
    }
 #endif
    if (!type.compare("standard")){
        this->residualVec_->update(-1.,*this->rhs_,1.);
        //if ( !this->sourceTerm_.is_null() )
        //    this->residualVec_->update(-1.,*this->sourceTerm_,1.);
    }
    else if(!type.compare("reverse")){
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs
        //if ( !this->sourceTerm_.is_null() )
        //    this->residualVec_->update(1.,*this->sourceTerm_,1.);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown type for residual computation.");
    }
    
    // this might be set again by the TimeProblem after adding of M*u
    this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );

}

template<class SC,class LO,class GO,class NO>
typename NonLinElasticity<SC,LO,GO,NO>::BlockMultiVectorPtr_Type NonLinElasticity<SC,LO,GO,NO>::getPostProcessingData() const
{
    BlockMultiVectorPtr_Type postProcess =Teuchos::rcp(new BlockMultiVector_Type(5)) ;
        
    /*
    0 -- "Volume","
    1 -- Sxx",
    2 -- "Sxy"
    3 -- "Sxz" 
    4 -- "Syx"
    5 -- "Syy"
    6 -- "Syz" 
    7 -- "Szx" 
    8 -- "Szy" 
    9 -- "Szz" 
    10 -- "MisesStress" 
    11 -- "SCirc"
    12 -- "SAxial",
    13 -- "SRadial"
    14 -- "Exx"
    15 -- "Exy"
    16 -- "Exz" 
    17 -- "Eyx"
    18 -- "Eyy"
    19 -- "Eyz" 
    20 -- "Ezx"
    21 -- "Ezy"
    22 -- "Ezz"
    23 -- "W"
    24 -- "Growth1"
    25 -- "Growth2" 
    26 -- "Growth3"
    27 -- "Stretch1"
    28 -- "Stretch2"
    29 -- "DetF",
    30 - 38 "Ag1n1","Ag1n2","Ag1n3","Ag2n1","Ag2n2","Ag2n3","Ag3n1","Ag3n2","Ag3n3"
    39 -- "a11"
    40 -- "a12"
    41 -- "a13"
    42 -- "a21"
    43 -- "a22"
    44 -- "a23"
    45 -- "nC1" <----- !!
    46 -- "nC2" <----- !!
    47 -- "nD1" <----- !! 
    48 -- "nD2" <----- !!
    49 -- "ScDir1"
    50 -- "ScDir2" 
    51 -- "ScDir3" 
    52 -- "SaDir1"
    53 -- "SaDir2"
    54 -- "SaDir3"
    55 -- "SrDir1"
    56 -- "SrDir2"
    57 -- "SrDir3"*/
    std::cout << "Assembling postprocessing data..." << std::endl;
    if(this->parameterList_->sublist("Parameter").get("SCI",false) == true || this->parameterList_->sublist("Parameter").get("FSCI",false) == true ){
    
        MultiVectorPtr_Type vonMisesStress = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapUnique() ));
        this->feFactory_->postProcessing(10, vonMisesStress);

        MultiVectorPtr_Type SCirc = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapUnique() ));
        this->feFactory_->postProcessing(11, SCirc);

        MultiVectorPtr_Type SAxial = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapUnique() ));
        this->feFactory_->postProcessing(12, SAxial);

        MultiVectorPtr_Type SRadial = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapUnique() ));
        this->feFactory_->postProcessing(13, SRadial);

        MultiVectorPtr_Type W = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapUnique() ));
        this->feFactory_->postProcessing(23, W);


        postProcess->addBlock(vonMisesStress,0);
        postProcess->addBlock(SCirc,1);
        postProcess->addBlock(SAxial,2);
        postProcess->addBlock(SRadial,3);
        postProcess->addBlock(W,4);
    }
    
    return postProcess;
}

template<class SC,class LO,class GO,class NO>
vec_string_Type NonLinElasticity<SC,LO,GO,NO>::getPostprocessingNames()
{
    return postProcessingnames_;
}

}
#endif
