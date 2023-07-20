#ifndef Mesh_decl_hpp
#define Mesh_decl_hpp

#define FULL_Mesh_TIMER

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Elements.hpp"
#include "feddlib/core/FE/Elements.hpp"
#include "feddlib/core/Mesh/AABBTree.hpp"

/*!
Defintion of Mesh

@brief  Mesh
@author Christian Hochmuth
@version 1.0
@copyright CH
*/

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Mesh {
    
public:
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    typedef Teuchos::RCP<Mesh> Mesh_ptr_Type;
    
    typedef Teuchos::RCP<Teuchos::Comm<int> > CommPtr_Type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;
    typedef const CommConstPtr_Type CommConstPtrConst_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtrConst_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    typedef AABBTree<SC,LO,GO,NO> AABBTree_Type;
    typedef Teuchos::RCP<AABBTree_Type > AABBTreePtr_Type;
    
    /* ###################################################################### */
    //
    Mesh();

    Mesh(CommConstPtrConst_Type& comm);
    
    ~Mesh();
    
    
    /*!
     Delete all member variables
     */
    void deleteData();
    
    void setParameterList( ParameterListPtr_Type& pL );

    ParameterListConstPtr_Type getParameterList( ) const;

    /// @brief Unique node map 
    /// @return mapUnique_
    MapConstPtr_Type getMapUnique() const;
    
    /// @brief Repeated node map 
    /// @return mapRepeated_
    MapConstPtr_Type getMapRepeated() const;

    /// @brief Unique P2 node map
    /// @return 
    MapConstPtr_Type getMapUniqueP2() const;
    
    /// @brief Repeated P2 node Map
    /// @return 
    MapConstPtr_Type getMapRepeatedP2() const;
    
    /// @brief Element Map. Always unique.
    /// @return elementMap_
    MapConstPtr_Type getElementMap();
	
    /// @brief Map of edges. Not unique. Always repeated
    /// @return edgeMap_
    MapConstPtr_Type getEdgeMap(); // Edge Map
    
    /// @brief 2D vector of repeated points. Full coordinates in each row
    /// @return pointsRep_
    vec2D_dbl_ptr_Type getPointsRepeated() const;

    /// @brief 2D vector of unique points. Full coordinates in each row
    /// @return pointsUnique_
    vec2D_dbl_ptr_Type getPointsUnique() const;
    
    /// @brief Vector of flags of repeated nodes. Corresponts to mapRepeated_ numbering
    /// @return bcFlagsRep_
    vec_int_ptr_Type getBCFlagRepeated() const;
    
    /// @brief Vector of flags of unique nodes. Corresponts to mapUnique_ numbering
    /// @return bcFlagsUnique_
    vec_int_ptr_Type getBCFlagUnique() const;
        
    virtual void dummy() = 0;
    
    /// @brief List of FiniteElements.cpp defined in Elements_Type.  
    /// @return elementsC_
    ElementsPtr_Type getElementsC();

    /// @brief  List of FiniteElements.cpp defined in Elements_Type. Here, FiniteElements represent dim-1 surfaces.
    /// @return surfaceElements_
    ElementsPtr_Type getSurfaceElements();
    
    /// @brief Dimension 
    /// @return dim
    int getDimension();
    
    /// @brief Global number of elements
    /// @return 
    GO getNumElementsGlobal();
    
    /// @brief Local number of elements (on the rank executing function)
    /// @return 
    LO getNumElements();
    
    /// @brief Number of points depending on unique or repeated view - default return is unique
    /// @param type Unqiue or Repeated
    /// @return 
    LO getNumPoints(std::string type="Unique");
    
    int getOrderElement();

    /// @brief 
    /// @return 
    CommConstPtrConst_Type getComm(){return comm_;};
    
    /// @brief 
    /// @param flags 
    /// @return 
    int setStructuredMeshFlags(int flags){return 0;};
    
    /// @brief 
    /// @param type 
    void setElementFlags(std::string type="");
    
    /// @brief Defining reference configuration. This needs to be executed before a 'moveMesh()' operation can be performed
    void setReferenceConfiguration();
    
    /// @brief 
    /// @param displacementUnique 
    /// @param displacementRepeated 
    void moveMesh( MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated );
    
    // Creates an AABBTree from own vertice- and elementlist.
    void create_AABBTree();
    
    /// @brief 
    /// @param query_points 
    /// @return 
    vec_int_ptr_Type  findElemsForPoints(vec2D_dbl_ptr_Type query_points);
    
    /// @brief 
    /// @param point 
    /// @param element 
    /// @return 
    vec_dbl_Type getBaryCoords(vec_dbl_Type point, int element);
    
    /// @brief 
    /// @param point 
    /// @param element 
    /// @return 
    bool isPointInElem(vec_dbl_Type point, int element);

    /// @brief 
    /// @return 
    tuple_intint_Type getRankRange() const {return rankRange_;};
    
    /// @brief 
    void deleteSurfaceElements(){ surfaceElements_.reset(); };
    

	/// @brief Returns elements as a vector type
	/// @return 
	vec2D_int_ptr_Type getElements();
    
    /* ###################################################################### */
    
    int                     dim_;
    long long               numElementsGlob_;

    std::string 			FEType_;
    MapPtr_Type             mapUnique_;
    MapPtr_Type 			mapRepeated_;
    vec2D_dbl_ptr_Type		pointsRep_;
    vec2D_dbl_ptr_Type 		pointsUni_;
    vec_int_ptr_Type 		bcFlagRep_;
    vec_int_ptr_Type		bcFlagUni_;

    ElementsPtr_Type        surfaceElements_;

    ElementsPtr_Type        elementsC_;
    MapPtr_Type				elementMap_;
    MapPtr_Type				edgeMap_;


    CommConstPtrConst_Type  comm_;
    
	vec2D_int_ptr_Type  elementsVec_;
    
    vec2D_dbl_ptr_Type		pointsRepRef_; // Repeated Referenzkonfiguration
    vec2D_dbl_ptr_Type		pointsUniRef_; // Unique Referenzkonfiguration
    
    //vec_int_ptr_Type 		elementFlag_;

    MapPtr_Type mapUniqueP2Map_;
    MapPtr_Type mapRepeatedP2Map_;
    
    ParameterListPtr_Type pList_;

    int elementOrder_;
    int surfaceElementOrder_;
    int edgesElementOrder_;

    AABBTreePtr_Type AABBTree_;
    
    tuple_intint_Type rankRange_;
    
    /* ###################################################################### */
private:

};
}

#endif
