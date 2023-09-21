#ifndef PrecOpFaCSCI_DEF_hpp
#define PrecOpFaCSCI_DEF_hpp
#include "PrecOpFaCSCI_decl.hpp"
#include <Thyra_TpetraMultiVector_decl.hpp>
#include <Teuchos_VerboseObject.hpp>
/*!
 Definition of PrecOpFaCSCI

 @brief  PrecOpFaCSCI
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
using namespace Thyra;
        
// Constructors

template<class SC, class LO, class GO, class NO>
PrecOpFaCSCI<SC,LO,GO,NO>::PrecOpFaCSCI()
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(false),
useSolidPreconditioner_(true)
{
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSCI, cant be used yet.");
}

template<class SC, class LO, class GO, class NO>
PrecOpFaCSCI<SC,LO,GO,NO>::PrecOpFaCSCI(CommConstPtr_Type comm, bool fluidPrecMonolithic, bool useFluidPreconditioner, bool useSolidPreconditioner, bool onlyDiagonal)
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(fluidPrecMonolithic),
useFluidPreconditioner_(useFluidPreconditioner),
useSolidPreconditioner_(useSolidPreconditioner),
onlyDiagonal_(onlyDiagonal)
{
    comm_=comm;
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSCI, cant be used yet.");
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setGE(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type sciInv,
                                     ThyraLinOpPtr_Type sciS,
                                     ThyraLinOpPtr_Type sciC,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT){

    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setSCIC(sciC);
    setSCIS(sciS);
    setSCIInv(sciInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);

    initialize();
}

    
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setC1(ThyraLinOpPtr_Type C1){
    C1_ = C1;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setC1T(ThyraLinOpPtr_Type C1T){
    C1T_ = C1T;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setC2(ThyraLinOpPtr_Type C2){
    C2_ = C2;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setC4(ThyraLinOpPtr_Type C4){
    C4_ = C4;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setShapeDeriv(ThyraLinOpPtr_Type shape_v, ThyraLinOpPtr_Type shape_p){
    shape_v_ = shape_v;
    shape_p_ = shape_p;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setSCIInv(ThyraLinOpPtr_Type sciInv){
    sciInv_ = sciInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setSCIC(ThyraLinOpPtr_Type sciC){
    sciC_ = sciC;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setSCIS(ThyraLinOpPtr_Type sciS){
    sciS_ = sciS;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setFluidInv(ThyraLinOpPtr_Type fInv){
    fInv_ = fInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setGeoInv(ThyraLinOpPtr_Type gInv){
    gInv_ = gInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setFluidF(ThyraLinOpPtr_Type fF){
    fF_ = fF;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::setFluidBT(ThyraLinOpPtr_Type fBT){
    fBT_ = fBT;
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::initialize(){

    std::cout << "  ########## Init PrecOpFaCSCI ########### " << std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION(fInv_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: Fluid preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(sciInv_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: Structure preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(C1_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: C1 not set.");
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( 5 );
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( 5 );
    vectorSpacesRange[0] = fF_->range();
    vectorSpacesRange[1] = fBT_->domain();
    vectorSpacesRange[2] = sciS_->range();
    vectorSpacesRange[3] = C1_->range();
    vectorSpacesRange[4] = sciC_->range();
    
    vectorSpacesDomain[0] = fF_->domain();
    vectorSpacesDomain[1] = fBT_->domain();
    vectorSpacesDomain[2] = sciS_->domain();
    vectorSpacesDomain[3] = C1T_->domain();
    vectorSpacesDomain[4] = sciC_->domain();

    if ( !gInv_.is_null() ) {
        vectorSpacesRange.push_back( gInv_->range() );
        vectorSpacesDomain.push_back( gInv_->domain() );
    }
    //     defaultProductRange_;
    //    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductDomain_;
    
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSR = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pR);
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSD = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pD);
//    
//    this->defaultProductRange_ = Thyra::multiVectorProductVectorSpace<SC>( pVSR , 1);
//    this->defaultProductDomain_ = Thyra::multiVectorProductVectorSpace<SC>( pVSD, 1);
    
    this->defaultProductRange_ = pR;
    this->defaultProductDomain_ = pD;
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::applyIt(
                                         const EOpTransp M_trans,
                                         const MultiVectorBase<SC> &X_in,
                                         const Ptr<MultiVectorBase<SC> > &Y_inout,
                                         const SC alpha,
                                         const SC beta
                                         ) const
{
    applyImpl(M_trans, X_in, Y_inout, alpha, beta);

}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::applyImpl(
                                   const EOpTransp M_trans,
                                   const MultiVectorBase<SC> &X_in,
                                   const Ptr<MultiVectorBase<SC> > &Y_inout,
                                   const SC alpha,
                                   const SC beta
                                   ) const
{
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    // X_in : Input Vector for applying to system. Defined operator as system being applied. Now we partition X into the different components. 
    // Different systems contain the preconditioners then.

   // std::cout << "  ########## Apply PrecOpFaCSCI ########### " << std::endl;
    using Teuchos::rcpFromRef;
    typedef Teuchos::ScalarTraits<SC> ST;
    typedef RCP<MultiVectorBase<SC> > MultiVectorPtr;
    typedef RCP<const MultiVectorBase<SC> > ConstMultiVectorPtr;
    typedef RCP<const LinearOpBase<SC> > ConstLinearOpPtr;

    int rank = comm_->getRank();
    
    // Getting the input and output vector of the whole FSI system. This is then split into the different components of fluid, sci, and geometry
    Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
        = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( rcpFromRef(X_in) );

    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
        = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );

    Y_inout->assign(0.);


    // SOLID : Solid at block component 2 
    Teuchos::RCP< const MultiVectorBase< SC > > X_s = X->getMultiVectorBlock(2);
    Teuchos::RCP< MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);
    assign(Y_s.ptr(), *X_s); 

    // Chemistry: Chemistry at block component 5
    Teuchos::RCP< const MultiVectorBase< SC > > X_chem = X->getMultiVectorBlock(4);
    Teuchos::RCP< MultiVectorBase< SC > > Y_chem = Y->getNonconstMultiVectorBlock(4);
    assign(Y_chem.ptr(), *X_chem);

   // Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XsTpetra =
   // Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_s );

    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_sci( 2 );
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_sci( 2 );

    X_sci[0] = Y_s;
    X_sci[1] = Y_chem;

    Y_sci[0] = Y_s;
    Y_sci[1] = Y_chem;
    
    //std::cout << "  ########## Apply PrecOpFaCSCI: SCI  " << std::endl;

    /*if (useSolidPreconditioner_){
    
        Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > sciMonoVS = sciInv_->domain();
        // std::cout << " Initi X_scimono_ " << std::endl;
        if ( X_scimono_.is_null() ){
            X_scimono_ = createMembers( sciMonoVS, X_sci[0]->domain()->dim() );
            Y_scimono_ = createMembers( sciMonoVS, Y_sci[0]->domain()->dim() );
        }
        //We can/should speedup this process
       // std::cout << " Copy to mono " << std::endl;

        //copyToMonoSCI(X_sci);
        //std::cout << " Apply " << std::endl;

        sciInv_->apply(NOTRANS, *X_scimono_, Y_scimono_.ptr(), 1., 0.);
        //copyFromMonoSCI(Y_sci); 
    }*/
    
    // apply solid preconditioner
    /*if (useSolidPreconditioner_)
        sInv_->apply(NOTRANS, *X_s, Y_s.ptr(), 1., 0.);
    else
        assign(Y_s.ptr(), *X_s);*/

    //std::cout << "  ########## Apply PrecOpFaCSCI: Geomertry  " << std::endl;

    Teuchos::RCP< const MultiVectorBase< SC > > X_g;
    Teuchos::RCP< MultiVectorBase< SC > > Y_g;
    // apply geometry preconditioner
    /*if ( !gInv_.is_null() ) {
        X_g = X->getMultiVectorBlock(4);
        Y_g = Y->getNonconstMultiVectorBlock(4);

        assign(Y_g.ptr(), *X_g);

        C4_->apply(NOTRANS, *Y_s, Y_g.ptr(), -1., 1.);
        
        gInv_->apply(NOTRANS, *Y_g, Y_g.ptr(), 1., 0.);
    }*/
    
    Teuchos::RCP< const MultiVectorBase< SC > > X_l = X->getMultiVectorBlock(3);
    Teuchos::RCP< MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);

    assign(Y_l.ptr(), *X_l);
    
    C2_->apply( NOTRANS, *Y_s, Y_l.ptr(), -1., 1. );

    
    Teuchos::RCP< const MultiVectorBase< SC > > X_fv = X->getMultiVectorBlock(0);
    Teuchos::RCP< MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
    assign(Y_fv.ptr(), *X_fv); 

    
    Teuchos::RCP< const MultiVectorBase< SC > > X_fp = X->getMultiVectorBlock(1);
    Teuchos::RCP< MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
    assign(Y_fp.ptr(), *X_fp);

    if (!shape_v_.is_null() && !shape_p_.is_null()) {
        shape_v_->apply(NOTRANS, *Y_g, Y_fv.ptr(), -1., 1.);
        shape_p_->apply(NOTRANS, *Y_g, Y_fp.ptr(), -1., 1.);
    }
    
    
    if (Z_fv_.is_null())
        Z_fv_ = Y_fv->clone_mv();
    else
        assign(Z_fv_.ptr(), *Y_fv);
    

    //std::cout << "  ########## Apply PrecOpFaCSCI: Fluid Condensation  " << std::endl;

    // fluid condensation
    if (tmp_l_.is_null())
        tmp_l_ = Y_l->clone_mv();

    C1_->apply(NOTRANS, *Y_fv, tmp_l_.ptr(), 1., 0);
    C1T_->apply(NOTRANS, *tmp_l_, Y_fv.ptr(), -1., 1.);
    
    C1T_->apply(NOTRANS, *Y_l, Y_fv.ptr(), 1., 1.);

    
    if (productRangeFluid_.is_null()) {
        Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRangeFluid( 2 );
        
        vectorSpacesRangeFluid[0] = X_fv->range();
        vectorSpacesRangeFluid[1] = X_fp->range();
        
        productRangeFluid_ = Thyra::productVectorSpace<SC>( vectorSpacesRangeFluid );
    }
    
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid( 2 );
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid( 2 );

    X_fluid[0] = Y_fv;
    X_fluid[1] = Y_fp;

    Y_fluid[0] = Y_fv;
    Y_fluid[1] = Y_fp;

   // std::cout << "  ########## Apply PrecOpFaCSCI: Fluid  " << std::endl;

    if (useFluidPreconditioner_){
    
        if (fluidPrecMonolithic_) {
            Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > fMonoVS = fInv_->domain();

            if ( X_fmono_.is_null() ){
                X_fmono_ = createMembers( fMonoVS, X_fluid[0]->domain()->dim() );
                Y_fmono_ = createMembers( fMonoVS, Y_fluid[0]->domain()->dim() );
            }
            //We can/should speedup this process
            //copyToMono(X_fluid);
            fInv_->apply(NOTRANS, *X_fmono_, Y_fmono_.ptr(), 1., 0.);
            //copyFromMono(Y_fluid);
        }
        else{
            Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodX_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, X_fluid );
            Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodY_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, Y_fluid );
            
            Teuchos::RCP< MultiVectorBase<SC> > X_f = Teuchos::rcp_dynamic_cast<MultiVectorBase< SC > >(prodX_f);
            Teuchos::RCP< MultiVectorBase<SC> > Y_f = Teuchos::rcp_dynamic_cast<MultiVectorBase< SC > >(prodY_f);
            
            fInv_->apply(NOTRANS, *X_f, Y_f.ptr(), 1., 0.);
        }
    }
    else{
        assign(Y_fv.ptr(), *X_fv);
        assign(Y_fp.ptr(), *X_fp);
    }

    //std::cout << "  ########## Apply PrecOpFaCSCI: Others  " << std::endl;


    fBT_->apply(NOTRANS, *Y_fp, Z_fv_.ptr(), -1., 1.);
 
    fF_->apply(NOTRANS, *Y_fv, Z_fv_.ptr(), -1., 1.);
    
    C1_->apply(NOTRANS, *Z_fv_, Y_l.ptr(), 1., 0.);
    
        
        
    Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_fvT =
    Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fv );
    Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_fpT =
    Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fp );
    Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_sT =
    Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_s );
    Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_lT =
    Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_l );
    Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_chemT =
    Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_chem );
        
        
    
}
// private
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::copyToMono( Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid ) const{
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[1]->range());
    

    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fmono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );

    std::cout << "Local Offset mono " << localOffset_mono << " local Sub dim mono " << localSubDim_mono << std::endl;


    for (int j=0; j<X_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_fluid)(i,j) = (*thyData_v)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_fluid)( i + localSubDim_v, j ) = (*thyData_p)(i,j);

    }

}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::copyToMonoSCI( Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid ) const{
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[1]->range());
    
    std::cout << "Copy to mono 1 " << std::endl;

    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    std::cout << "Local Offset v " << localOffset_v << " local Sub dim v " << localSubDim_v << std::endl;

    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    std::cout << "Local Offset p " << localOffset_p << " local Sub dim p " << localSubDim_p << std::endl;


        std::cout << "Copy to mono 2 " << std::endl;
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
        std::cout << "Copy to mono 3 " << std::endl;

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_scimono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    std::cout << "Local Offset mono " << localOffset_mono << " local Sub dim mono " << localSubDim_mono << std::endl;


    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_scimono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );

    std::cout << "Copy to mono 4 " << std::endl;

    for (int j=0; j<X_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_fluid)(i,j) = (*thyData_v)(i,j);
        
        std::cout << "Copy to mono 5 " << std::endl;


        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_fluid)( i + localSubDim_v, j ) = (*thyData_p)(i,j);

    }
        std::cout << "Copy to mono 6 " << std::endl;

}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::copyFromMono(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid) const{
    

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fmono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );
    
    for (int j=0; j<Y_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_v)(i,j) = (*thyData_fluid)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_p)(i,j) = (*thyData_fluid)( i + localSubDim_v, j );
        
    }    
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSCI<SC,LO,GO,NO>::copyFromMonoSCI(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid) const{
    

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_scimono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_scimono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );
    
    for (int j=0; j<Y_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_v)(i,j) = (*thyData_fluid)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_p)(i,j) = (*thyData_fluid)( i + localSubDim_v, j );
        
    }    
}

}

#endif
