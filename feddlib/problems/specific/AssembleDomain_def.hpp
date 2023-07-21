#ifndef AssembleDomain_def_hpp
#define AssembleDomain_def_hpp
#include "AssembleDomain_decl.hpp"
/*!
 Definition of NonLinElasticity
 
 @brief NonLinElasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template<class SC,class LO,class GO,class NO>
AssembleDomain<SC,LO,GO,NO>::AssembleDomain(const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList, string problemType):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domain->getComm() ),
sol_rep_()
{
    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();
    this->addVariable( domain , FEType , "u" , domain->getDofs());

    this->dim_ = this->getDomain(0)->getDimension();
       
    problemType_=problemType;   
    domainVec_.push_back(domain);
    problemSize_=1;

    sol_rep_ = Teuchos::rcp( new BlockMultiVector_Type(1));
    for(int probRow = 0; probRow <  problemSize_; probRow++){
        MapConstPtr_Type rowMap;
        if(domainVec_.at(probRow)->getDofs() == 1)
            rowMap = domainVec_.at(probRow)->getMapRepeated();
        else
            rowMap = domainVec_.at(probRow)->getMapVecFieldRepeated();

        MultiVectorPtr_Type sol_rep = Teuchos::rcp( new MultiVector_Type( rowMap ) );
        sol_rep_->addBlock(sol_rep,probRow);
    }

    loadStepping_ =    !(parameterList->sublist("Timestepping Parameter").get("Class","Singlestep")).compare("Loadstepping");
    externalForce_ =   parameterList->sublist("Parameter").get("External Force",false);
    nonlinearExternalForce_ = parameterList->sublist("Parameter").get("Nonlinear External Force",false);

    timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));

}

template<class SC,class LO,class GO,class NO>
AssembleDomain<SC,LO,GO,NO>::~AssembleDomain(){

}

template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}
    
template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::assemble(std::string type) const{
    
    if (type == ""){
        if (this->verbose_)
            std::cout << "-- Assembly Domain ... " << std::flush;     
        
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
        
        this->reAssemble("Newton");

        this->assembleSourceTerm(0.);

        //if(sourceType == "volume")
        //    this->sourceTerm_->scale(density);

        this->addToRhs( this->sourceTerm_ );
        
        this->setBoundariesRHS();
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
void AssembleDomain<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
}

template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::reAssemble(std::string type) const {

    if (this->verbose_)
        std::cout << "-- Reassembly " << problemType_ << " " << std::flush;
    
    for(int probRow = 0; probRow <  problemSize_; probRow++){
        for(int probCol = 0; probCol <  problemSize_; probCol++){
            MapConstPtr_Type rowMap;

            if(domainVec_.at(probRow)->getDofs() == 1)
                rowMap = domainVec_.at(probRow)->getMapUnique();
            else
                rowMap = domainVec_.at(probRow)->getMapVecFieldUnique();

            MatrixPtr_Type A = Teuchos::rcp(new Matrix_Type( rowMap, domainVec_[probRow]->getDimension() * domainVec_[probRow]->getApproxEntriesPerRow() ) );

            this->system_->addBlock(A,probRow,probCol);

        }
	}

    for(int probRow = 0; probRow <  problemSize_; probRow++){
        MapConstPtr_Type rowMap;
        if(domainVec_.at(probRow)->getDofs() == 1)
            rowMap = domainVec_.at(probRow)->getMapUnique();
        else
            rowMap = domainVec_.at(probRow)->getMapVecFieldUnique();

        MultiVectorPtr_Type resVecUnique = Teuchos::rcp( new MultiVector_Type( rowMap, 1 ) );
        resVecUnique->putScalar(0.);
        this->residualVec_->addBlock(resVecUnique,probRow);

        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        sol_rep_->getBlock(probRow)->importFromVector(u, true);

    }
   
    this->feFactory_->globalAssembly(problemType_, this->dim_, 0, sol_rep_, this->system_, this->residualVec_,this->getParameterList(),"Jacobian");
  
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}
    


template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                            const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                            ) const
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;
    
    TEUCHOS_TEST_FOR_EXCEPTION( this->solution_->getBlock(0)->getMap()->getUnderlyingLib() != "Tpetra", std::runtime_error, "Use of NOX only supports Tpetra. Epetra support must be implemented.");
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");
    
    RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);
    
    this->solution_->fromThyraMultiVector(vecThyraNonConst);
    
    const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();
    
    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
    typedef Xpetra::Matrix<SC,LO,GO,NO> XpetraMatrix_Type;
    typedef RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;
    typedef RCP<const XpetraMatrix_Type> XpetraMatrixConstPtr_Type;
    
    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);
    
    if ( fill_f || fill_W || fill_W_prec ) {
        
        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {
            
            this->calculateNonLinResidualVec("standard");
            
            RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }
        
        XpetraMatrixPtr_Type W;
        if (fill_W) {
            
            this->reAssemble("Newton");
            
            this->setBoundariesSystem();
            
            RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            RCP<TpetraMatrix_Type> W_tpetraMat = rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
            
            XpetraMatrixConstPtr_Type W_systemXpetra = this->getSystem()->getBlock( 0, 0 )->getXpetraMatrix();
            
            XpetraMatrixPtr_Type W_systemXpetraNonConst = rcp_const_cast<XpetraMatrix_Type>(W_systemXpetra);
            Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
            Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatXpetra = xTpetraMat.getTpetra_CrsMatrixNonConst();
            
            W_tpetraMat->resumeFill();
            
            for (auto i=0; i<tpetraMatXpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                tpetraMatXpetra->getLocalRowView( i, indices, values);
                W_tpetraMat->replaceLocalValues( i, indices, values);
            }
            W_tpetraMat->fillComplete();
            
        }
        
        if (fill_W_prec) {
            this->setupPreconditioner( "Monolithic" );
            
            // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
            // If this is the case, we need to pre apply the coarse level to the residual(f_out).
            
            std::string levelCombination = this->parameterList_->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive");
            if (!levelCombination.compare("Multiplicative")) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX. In general we need to pre-apply the coarse problem. This must be implemented here.");
            }
            
        }
    }
}

template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::assembleSourceTermLoadstepping(double time) const
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
                this->feFactory_->assemblyNonlinearSurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, sol_rep_->getBlock(0),Kext, funcParameter, this->rhsFuncVec_[0],this->parameterList_);
               
                A->addMatrix(1.,AKext,0.);
                Kext->addMatrix(1.,AKext,1.);

                AKext->fillComplete(this->getDomain(0)->getMapVecFieldUnique(),this->getDomain(0)->getMapVecFieldUnique());

                this->system_->addBlock(AKext,0,0);

            }
            else            
                this->feFactory_->assemblySurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, sol_rep_->getBlock(0), funcParameter, this->rhsFuncVec_[0],this->parameterList_);


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
Teuchos::RCP<Thyra::LinearOpBase<SC> > AssembleDomain<SC,LO,GO,NO>::create_W_op() const
{
    
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > AssembleDomain<SC,LO,GO,NO>::create_W_prec() const
{
    this->initializeSolverBuilder();
    this->initializePreconditioner();
    
    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst= Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);
    
    return thyraPrecNonConst;
}

template<class SC,class LO,class GO,class NO>
void AssembleDomain<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    
    this->feFactory_->globalAssembly(problemType_, this->dim_, 0, sol_rep_, this->system_, this->residualVec_,this->getParameterList(),"Rhs");

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
}
#endif
