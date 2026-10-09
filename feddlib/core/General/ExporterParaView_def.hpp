#ifndef ExporterParaView_DEF_hpp
#define ExporterParaView_DEF_hpp

#include "OutputHistory.hpp"
#include "HDF5VectorInfo.hpp"
#include <iomanip>

/*!
 Definition of ExporterParaView

 @brief  ExporterParaView
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template<class SC,class LO,class GO,class NO>
ExporterParaView<SC,LO,GO,NO>::ExporterParaView():
hdf5exporter_(),
comm_(),
closingLinesPosition_(),
closingLinesPositionTimes_(),
closingLines_(),
xmf_out_(),
xmf_times_out_(),
filename_(),
outputFilename_(),
postfix_(),
FEType_(),
variables_(0),
uniqueMaps_(0),
varNames_(0),
varTypes_(0),
varDofPerNode_(0),
pointsHDF_(),
elementsHDF_(),
dim_(0),
nmbElementsGlob_(0),
nmbPointsGlob_(0),
nmbExportValuesGlob_(0),
nmbPointsPerElement_(0),
timeIndex_(0),
writeDt_(0),
saveTimestep_(1),
verbose_(false),
parameterList_(),
pointsUnique_()
{

}
    
template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::setup(std::string filename,
                                          MeshPtr_Type mesh,
                                          std::string FEType,
                                          int saveTimestep,
                                          ParameterListPtr_Type parameterList){

    
    setup( filename, mesh, FEType, parameterList);
    saveTimestep_ = saveTimestep;
    TEUCHOS_TEST_FOR_EXCEPTION(saveTimestep_ < 1, std::logic_error, "Export interval must be positive");
}
  
template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::setup(std::string filename,
                                          MeshPtr_Type mesh,
                                          std::string FEType,
                                          ParameterListPtr_Type parameterList){
    
    
    parameterList_ = parameterList;
    comm_ = mesh->getComm();
    verbose_ = (comm_->getRank() == 0);
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( mesh->getComm() );
    filename_ = filename;
    outputFilename_ = filename_ + ".h5";
    FEType_ = FEType;
    dim_ = mesh->getDimension();
    nmbElementsGlob_ = mesh->getNumElementsGlobal();
    
    pointsUnique_ = mesh->getPointsUnique();
    
    if (FEType_ == "P0") {
        if (dim_ == 2 ) {
            nmbPointsPerElement_ = 3 ;
            
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 4 ;
            
        }
    }
    else if (FEType_ == "P1") {
        if (dim_ == 2 ) {
            nmbPointsPerElement_ = 3 ;
            
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 4 ;
            
        }
    }
    else if (FEType_ == "P1-disc") {
        if (dim_ == 2 ) {
            nmbPointsPerElement_ = 3 ;
            
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 4 ;
            
        }
    }
    else if(FEType_ == "P2"){
        
        if (dim_ == 2 ) {
            nmbPointsPerElement_ = 6 ;
            
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 10 ;
        }
    }
    else if(FEType_ == "P2-CR"){
        if (dim_ == 2 ) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong dimension for P2-CR.");
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 10 ; // only export P2 points. Paraview 5.6.0 and prior versions ignore face centered values.
        }
    }
    else if(FEType_ == "Q1"){
        if (dim_ == 2 ) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong dimension for Q1.");
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 8 ;
        }
    }
    
    else if(FEType_ == "Q2"){
        if (dim_ == 2 ) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong dimension for Q2.");
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 20 ; // only export Q20 points.
        }
    }
    else if(FEType_ == "Q2-20"){
        if (dim_ == 2 ) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong dimension for Q2-20.");
        }
        else if(dim_ == 3){
            nmbPointsPerElement_ = 20 ; // only export Q20 points.
        }
    }
    else  {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong FEType, choose either P0 or P1 or P1-disc or P2 or P2-CR or Q1 or Q2 or Q2-20. Subdomain export with P0 - Export stopped");
    }
    
    nmbPointsGlob_ = mesh->getMapUnique()->getGlobalNumElements();
    
    closingLines_ = "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n";
    
    timeIndex_ = 0;
    
    // Something different happens to the element List and the elements
    // Probably they have the following form:
    // ElementMap :   0     1       2       3       4       5  
    // ->            0 1 2  3 4 5   6 7 8   ...
    // Element GIDs:  0 1 2 3 ...
    // Where the map tells us which IDs belong to which element
    MapConstPtr_Type elementMap = mesh->getElementMap();
    Teuchos::ArrayView<const GO> nodeElementList = elementMap->getNodeElementList(); // Global Ids on this processor
    vec_GO_Type nodeElementListInteger( nmbPointsPerElement_ * nodeElementList.size() );
    int counter=0;
    for (int i=0; i<nodeElementList.size(); i++) { // Number of elements
        for (int j=0; j<nmbPointsPerElement_; j++){ // number of points per element
            nodeElementListInteger[counter] = (int) nmbPointsPerElement_*nodeElementList[i] + j;
            counter++; // from 0 ... to nmbPointsPerElement_*localElements
        }
    }
    Teuchos::ArrayView<GO> globalMapIDs = Teuchos::arrayViewFromVector( nodeElementListInteger);
    MapPtr_Type	mapElements = Teuchos::rcp( new Map_Type((int) (nmbPointsPerElement_*nmbElementsGlob_), globalMapIDs,elementMap->getIndexBase()*nmbPointsPerElement_, comm_));

    // They contain global IDs of nodes corresponding to 'elements'
    elementsHDF_.reset(new MultiVector_Type(mapElements,1));
    
    ElementsPtr_Type elements = mesh->getElementsC();
    counter = 0;
    for (int i=0; i<elements->numberElements(); i++) {
        for (int j=0; j<nmbPointsPerElement_; j++) {
            int globalIndex = (int) mesh->getMapRepeated()->getGlobalElement( elements->getElement(i).getNode(j) );
            (elementsHDF_->getDataNonConst(0))[counter] = globalIndex;
            counter++;
        }
    }

    
    pointsHDF_.reset(new MultiVector_Type(mesh->getMapUnique(),dim_));

    updatePoints();
    
    writeDt_ = false;
    
    saveTimestep_ = parameterList_.is_null() ? 1 :
        parameterList_->sublist("Exporter").get("Export every X timesteps", 1);
    TEUCHOS_TEST_FOR_EXCEPTION(saveTimestep_ < 1, std::logic_error, "Export interval must be positive");
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::initializeOutput()
{
    if (outputInitialized_) return;
    output::archive(parameterList_, *comm_, {outputFilename_, filename_ + ".xmf", filename_ + "_times.xmf"});
    int existing = 0, nextIndex = 0, hasTimes = 0;
    double lastTime = -std::numeric_limits<double>::infinity();
    std::string retained, retainedTimes;
    if (output::resume(parameterList_)) {
        output::onRank(*comm_, 0, [&] {
            const bool h5Exists = std::filesystem::exists(outputFilename_);
            const bool xmfExists = std::filesystem::exists(filename_ + ".xmf");
            if (h5Exists != xmfExists) throw std::runtime_error("Both HDF5 and XMF files are required to resume " + filename_);
            if (!h5Exists) return; // Output may have been disabled in the original run.
            const double time = output::restartTime(parameterList_);
            const output::Collection collection(output::read(filename_ + ".xmf"));
            const double dt = parameterList_->sublist("Timestepping Parameter").get("dt", 0.);
            if (!(dt > 0.)) throw std::runtime_error("Resume output requires a positive dt");
            std::smatch cadence;
            const std::regex interval("Name=\"FEDD export interval\" Value=\"([0-9]+)\"");
            double oldInterval = saveTimestep_;
            if (std::regex_search(collection.header, cadence, interval))
                oldInterval = std::stoi(cadence[1]);
            else if (collection.frames.size() > 1) {
                const auto& first = collection.frames.front();
                const auto& last = collection.frames.back();
                oldInterval = (last.time - first.time) / ((last.index - first.index) * dt);
            }
            if (std::abs(oldInterval - saveTimestep_) > 1.e-7)
                throw std::runtime_error("Export interval differs from the existing output");
            if (collection.frames.front().xml.find("<Grid Name=\"Mesh" + FEType_ + " ") == std::string::npos)
                throw std::runtime_error("FE type differs from the existing output");
            const double origin = collection.frames.front().time - collection.frames.front().index * saveTimestep_ * dt;
            const double index = (time - origin) / dt;
            if (index < 0. || std::abs(index - std::round(index)) > 1.e-7 || index > std::numeric_limits<int>::max())
                throw std::runtime_error("Restart time does not match the existing output timestep grid");
            nextIndex = static_cast<int>(std::round(index));
            std::set<std::string> keep;
            for (const auto& frame : collection.frames) {
                if (frame.time > time + 1.e-12) continue;
                lastTime = frame.time;
                keep.insert(frame.groups.begin(), frame.groups.end());
            }
            // Check the old layout before deleting any superseded frame.
            const auto& layout = collection.frames.front().xml;
            std::set<std::string> savedVariables;
            const std::regex attribute("<Attribute[^>]*Name=\"([^\"]+)\"");
            for (std::sregex_iterator i(layout.begin(), layout.end(), attribute), end; i != end; ++i)
                savedVariables.insert((*i)[1]);
            if (savedVariables != std::set<std::string>(varNames_.begin(), varNames_.end()))
                throw std::runtime_error("Existing ParaView fields do not match the resumed exporter");
            const bool movingMesh = parameterList_->sublist("Exporter").get("Write new mesh", false);
            if ((collection.frames.front().groups.count("PointsX") == 0) != movingMesh)
                throw std::runtime_error("Write new mesh differs from the existing output");
            retained = collection.through(time);
            hasTimes = std::filesystem::exists(filename_ + "_times.xmf");
            if (hasTimes) retainedTimes = output::Collection(output::read(filename_ + "_times.xmf")).through(time);
            checkpoint::H5Handle file(H5Fopen(outputFilename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
            if (file < 0) throw std::runtime_error("Cannot open output " + outputFilename_);
            const auto checkShape = [&](const std::string& group, hsize_t entries) {
                checkpoint::H5Handle dataset(H5Dopen(file, (group + "/Values").c_str(), H5P_DEFAULT), H5Dclose);
                if (dataset < 0) throw std::runtime_error("Missing output dataset " + group);
                checkpoint::H5Handle space(H5Dget_space(dataset), H5Sclose);
                if (space < 0 || H5Sget_simple_extent_npoints(space) != static_cast<hssize_t>(entries))
                    throw std::runtime_error("Existing output layout differs for " + group);
            };
            checkShape("Connections", static_cast<hsize_t>(nmbElementsGlob_) * nmbPointsPerElement_);
            for (const auto& frame : collection.frames) {
                if (frame.time > time + 1.e-12) continue;
                std::ostringstream index;
                index << '.' << std::setfill('0') << std::setw(5) << frame.index;
                for (const auto& group : frame.groups) {
                    hsize_t entries = nmbPointsGlob_;
                    if (group == "Connections") entries = static_cast<hsize_t>(nmbElementsGlob_) * nmbPointsPerElement_;
                    for (std::size_t i = 0; i < varNames_.size(); ++i)
                        if (group == varNames_[i] + index.str())
                            entries = static_cast<hsize_t>(uniqueMaps_[i]->getGlobalNumElements()) *
                                (varTypes_[i] == "Vector" ? 3 : 1);
                    checkShape(group, entries);
                }
            }
            // Remove orphaned frame groups too: a failed write may have flushed
            // HDF5 datasets without publishing their frame in the XMF collection.
            std::vector<std::string> groups;
            const auto collect = [](hid_t, const char* name, void* data) -> herr_t {
                static_cast<std::vector<std::string>*>(data)->push_back(name); return 0;
            };
            if (H5Giterate(file, "/", nullptr, collect, &groups) < 0)
                throw std::runtime_error("Cannot enumerate output frames");
            const std::regex frameGroup(".*\\.[0-9]+|Points[XYZ][0-9]+|Connections[0-9]+");
            for (const auto& group : groups)
                if (!keep.count(group) && std::regex_match(group, frameGroup) &&
                    H5Ldelete(file, group.c_str(), H5P_DEFAULT) < 0)
                    throw std::runtime_error("Cannot remove superseded output frame " + group);
            output::replace(filename_ + ".xmf", retained + closingLines_);
            if (hasTimes) output::replace(filename_ + "_times.xmf", retainedTimes + closingLines_);
            existing = 1;
        });
    }
    Teuchos::broadcast(*comm_, 0, 1, &existing);
    if (existing) {
        Teuchos::broadcast(*comm_, 0, 1, &nextIndex);
        Teuchos::broadcast(*comm_, 0, 1, &lastTime);
        Teuchos::broadcast(*comm_, 0, 1, &hasTimes);
        timeIndex_ = nextIndex;
        lastOutputTime_ = lastTime;
        resumeOutput_ = true;
        writeDt_ = hasTimes != 0;
        hdf5exporter_.reset(new HDF5_Type(comm_));
        hdf5exporter_->open(outputFilename_);
        output::onRank(*comm_, 0, [&] {
            xmf_out_.open(filename_ + ".xmf", std::ios::in | std::ios::out);
            if (!xmf_out_) throw std::runtime_error("Cannot reopen XMF output");
            xmf_out_ << std::setprecision(17);
            closingLinesPosition_ = static_cast<std::streamoff>(retained.size());
            if (hasTimes) {
                xmf_times_out_.open(filename_ + "_times.xmf", std::ios::in | std::ios::out);
                if (!xmf_times_out_) throw std::runtime_error("Cannot reopen XMF time output");
                xmf_times_out_ << std::setprecision(17);
                closingLinesPositionTimes_ = static_cast<std::streamoff>(retainedTimes.size());
            }
        });
    } else {
        initHDF5();
        output::onRank(*comm_, 0, [&] { initXmf(); });
    }
    outputInitialized_ = true;
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::addVariable(MultiVecConstPtr_Type &u,
                                                  std::string varName,
                                                  std::string varType,
                                                  int dofPerNode,
                                                  MapConstPtrConst_Type& mapUnique){

    variables_.push_back(u);
    varNames_.push_back(varName);
    varTypes_.push_back(varType);
    varDofPerNode_.push_back(dofPerNode);

    nmbExportValuesGlob_ = mapUnique->getGlobalNumElements();

    uniqueMaps_.push_back(mapUnique);

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::save(double time){

    initializeOutput();
    if (resumeOutput_ && time <= lastOutputTime_ + 1.e-12) { ++timeIndex_; return; }

    if (timeIndex_ % saveTimestep_ == 0) {
        makePostfix();

        if (!parameterList_.is_null()){
            if ( parameterList_->sublist("Exporter").get("Write new mesh",false) )
                writeMeshPointsHDF5();
        }

        writeVariablesHDF5();

        writeXmf(time);
        if (writeDt_) writeXmfTime(time, parameterList_.is_null() ? lastDt_ :
            parameterList_->sublist("Timestepping Parameter").get("dt", lastDt_));
        lastOutputTime_ = time;
    }
    else{
        if (this->verbose_)
            std::cout << "\n \t ### Export criterion not satisfied for current time step - no export of time step! ###" << std::endl;
    }

    timeIndex_++;

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::save(double time, double dt){

    lastDt_ = dt;
    initializeOutput();
    if (resumeOutput_ && time <= lastOutputTime_ + 1.e-12) { ++timeIndex_; return; }

    if (timeIndex_ % saveTimestep_ == 0) {

        makePostfix();
        if (!parameterList_.is_null()){
            if ( parameterList_->sublist("Exporter").get("Write new mesh",false) )
                writeMeshPointsHDF5();
        }
        writeVariablesHDF5();

        writeXmf(time);

        if (!writeDt_) {
            initXmfTimes();
        }
        writeXmfTime(time, dt);
        lastOutputTime_ = time;
    }
    else{
        if (this->verbose_)
            std::cout << "\n \t ### Export criterion not satisfied for current time step - no export of time step! ###" << std::endl;
    }

    timeIndex_++;

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::closeExporter(){

    if (!hdf5exporter_.is_null()) hdf5exporter_->close();
    xmf_out_.close();
    if (writeDt_) {
        xmf_times_out_.close();
    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::initHDF5(){

    hdf5exporter_.reset( new HDF5_Type(comm_) );
    hdf5exporter_->create(outputFilename_);
    std::string nameConn = "Connections";

    writeMeshElements(nameConn);
    // Mesh is only written once. If we have a new mesh for a new export, we call writeMeshPointsHDF5() in function save(...), when all the data is updated.
    if (parameterList_.is_null())
        writeMeshPointsHDF5();
    else if ( parameterList_->sublist("Exporter").get("Write new mesh",false) == false )
        writeMeshPointsHDF5();

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeMeshPointsHDF5(){
    if (!parameterList_.is_null()){
        if ( parameterList_->sublist("Exporter").get("Write new mesh",false) ) {
            updatePoints();
            std::string nameP_X = "PointsX" + std::to_string(timeIndex_);
            std::string nameP_Y = "PointsY" + std::to_string(timeIndex_);
            std::string nameP_Z = "PointsZ" + std::to_string(timeIndex_);

            writeMeshPoints(nameP_X, nameP_Y, nameP_Z );
        }
    
        else{
            std::string nameP_X = "PointsX";
            std::string nameP_Y = "PointsY";
            std::string nameP_Z = "PointsZ";
            writeMeshPoints( nameP_X, nameP_Y, nameP_Z );
        }
    }
    else{
        std::string nameP_X = "PointsX";
        std::string nameP_Y = "PointsY";
        std::string nameP_Z = "PointsZ";
        writeMeshPoints( nameP_X, nameP_Y, nameP_Z );
    }
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeMeshElements( std::string nameConn ){
    //Triangle/Tetrahedron Connections
    hdf5exporter_->createGroup(nameConn);

    hdf5exporter_->write(nameConn,elementsHDF_);
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeMeshPoints(std::string nameP_X,
                                                    std::string nameP_Y,
                                                    std::string nameP_Z){



    hdf5exporter_->createGroup(nameP_X);
    hdf5exporter_->createGroup(nameP_Y);

    hdf5exporter_->write(nameP_X,pointsHDF_->getVector(0));
    hdf5exporter_->write(nameP_Y,pointsHDF_->getVector(1));

    if (dim_ == 3) {
        hdf5exporter_->createGroup(nameP_Z);
        hdf5exporter_->write(nameP_Z,pointsHDF_->getVector(2));
    }
    hdf5exporter_->flush();

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::updatePoints(){
    int dim = -1;
    if (pointsUnique_->size()>0)
        dim = pointsUnique_->at(0).size();

    for (int i=0; i<pointsUnique_->size(); i++) {
        for (int j = 0; j < dim; j++) {
            pointsHDF_->getDataNonConst(j)[i] = (*pointsUnique_)[i][j];
        }
    }
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeVariablesHDF5(){

    for (int i=0; i<variables_.size(); i++) {

        if (varTypes_.at(i)=="Vector") {
            MultiVectorPtr_Type u_export(new MultiVector_Type(uniqueMaps_.at(i),3)); // ParaView always uses 3D Data. for 2D Data the last entries (for the 3rd Dim) are all zero.

            hdf5exporter_->createGroup(varNames_[i]+postfix_);

            prepareVectorField(variables_.at(i),u_export, varDofPerNode_.at(i));

            hdf5exporter_->write(varNames_[i]+postfix_,u_export,true);
        }
        else if(varTypes_.at(i)=="Scalar"){
            MultiVectorPtr_Type u_export(new MultiVector_Type(uniqueMaps_.at(i),1)); // ParaView always uses 3D Data. for 2D Data the last entries (for the 3rd Dim) are all zero.

            hdf5exporter_->createGroup(varNames_[i]+postfix_);

            prepareScalar(variables_.at(i),u_export); //conversion to int

            hdf5exporter_->write(varNames_[i]+postfix_,u_export);
        }
        hdf5exporter_->flush();
    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::initXmf(){

    if (comm_->getRank()==0) {

        xmf_out_.open((filename_ + ".xmf").c_str(),std::ios_base::out);
        if (!xmf_out_) throw std::runtime_error("Cannot create XMF output");
        xmf_out_ << std::setprecision(17);
        xmf_out_ 	<< "<?xml version=\"1.0\" ?>\n"
                    << "<!DOCTYPE Xdmf SYSTEM \""
                    << filename_
                    << ".xdmf\" [\n"
                    << "<!ENTITY DataFile \""
                    << filename_
                    << ".h5\">\n"
                    << "]>\n"
                    << "<!-- "
                    << filename_
                    << ".h5  -->\n"
                    << "<Xdmf>\n"
                    << "  <Domain Name=\""
                    << filename_
                    << "\">\n"
                    << "    <Grid Name=\""
                    << filename_
                    << "Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
                    << "    <Information Name=\"FEDD export interval\" Value=\"" << saveTimestep_ << "\" />\n"
                    << "\n";

        closingLinesPosition_ = xmf_out_.tellp();
     // write closing lines
        xmf_out_ << closingLines_;
        xmf_out_.flush();
    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::initXmfTimes(){
    writeDt_ = true;
    if (comm_->getRank()==0) {

        xmf_times_out_.open((filename_ + "_times.xmf").c_str(),std::ios_base::out);
        xmf_times_out_ 	<< "<?xml version=\"1.0\" ?>\n"
        << "<!DOCTYPE Xdmf SYSTEM \""
        << filename_
        << ".xdmf\" []>\n"
        << "<Xdmf>\n"
        << "  <Domain Name=\""
        << filename_
        << "\">\n"
        << "    <Grid Name=\""
        << filename_
        << "Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
        << "\n";

        closingLinesPositionTimes_ = xmf_times_out_.tellp();
        // write closing lines
        xmf_times_out_ << closingLines_;
        xmf_times_out_.flush();
    }

}
template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeXmf(double time){

    std::string nameP_X;
    std::string nameP_Y;
    std::string nameP_Z;
    std::string nameConn = "Connections";
	if(redo_ == true)
    	nameConn = "Connections" + std::to_string(timeIndex_);
    if (!parameterList_.is_null()){
        if ( parameterList_->sublist("Exporter").get("Write new mesh",false) ) {
            nameP_X = "PointsX" + std::to_string(timeIndex_);
            nameP_Y = "PointsY" + std::to_string(timeIndex_);
            nameP_Z = "PointsZ" + std::to_string(timeIndex_);
        }
        else{
            nameP_X = "PointsX";
            nameP_Y = "PointsY";
            nameP_Z = "PointsZ";
        }
    }
    else {
        nameP_X = "PointsX";
        nameP_Y = "PointsY";
        nameP_Z = "PointsZ";
    }
    writeXmfElements( nameConn, time);
    writeXmfPoints( nameP_X, nameP_Y, nameP_Z );
    writeXmfVariables();

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeXmfElements( std::string nameConn, double time ){

    if (verbose_) {
        xmf_out_.seekp (closingLinesPosition_);

        xmf_out_ <<
        "<!-- Time " << time << " Iteration " << postfix_.substr (1, 5) << " -->\n" <<
        "    <Grid Name=\"Mesh" << FEType_ << " " << time << "\">\n" <<
        "      <Time TimeType=\"Single\" Value=\"" << time << "\" />\n";

        //        writeTopology (M_xdmf);
        std::string FEstring;
        if (FEType_=="P0") {
            if (dim_ == 2) {
                FEstring = "Triangle";
            }
            else if(dim_ == 3){
                FEstring = "Tetrahedron";
            }
        }
        else if (FEType_=="P1-disc") {
            if (dim_ == 2) {
                FEstring = "Triangle";
            }
            else if(dim_ == 3){
                FEstring = "Tetrahedron";
            }
        }
        else if (FEType_=="P1") {
            if (dim_ == 2) {
                FEstring = "Triangle";
            }
            else if(dim_ == 3){
                FEstring = "Tetrahedron";
            }
        }

        else if (FEType_=="P2"){
            if (dim_ == 2) {
                FEstring = "Tri_6";
            }
            else if(dim_ == 3){
                FEstring = "Tet_10";
            }
        }
        else if (FEType_=="P2-CR"){
            //there is no 2D version
            if(dim_ == 3){
                FEstring = "Tet_10";
            }
        }
        else if (FEType_=="Q1"){
            //there is no 2D version
            if(dim_ == 3){
                FEstring = "Hexahedron";
            }
        }
        else if (FEType_=="Q2"){
            //there is no 2D version
            if(dim_ == 3){
                FEstring = "Hex_20";
            }
        }
        else if (FEType_=="Q2-20"){
            //there is no 2D version
            if(dim_ == 3){
                FEstring = "Hex_20";
            }
        }

        xmf_out_	  << "      <Topology\n"
        << "         Type=\""
        << FEstring
        << "\"\n"
        << "         NumberOfElements=\""
        << nmbElementsGlob_
        << "\"\n"
        << "         BaseOffset=\""
        << 0
        << "\">\n"
        << "         <DataStructure Format=\"HDF\"\n"
        << "                        Dimensions=\""
        << nmbElementsGlob_//this->M_mesh->numGlobalElements()
        << " "
        << nmbPointsPerElement_//this->M_mesh->numLocalVertices()
        << "\"\n" << "                        DataType=\"Int\"\n"
        << "                        Precision=\"8\">\n" <<  "             "
        << outputFilename_ << ":/"<< nameConn <<"/Values\n"
        << "         </DataStructure>\n" << "      </Topology>\n";

        closingLinesPosition_ = xmf_out_.tellp();

    }
}
template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeXmfPoints(std::string nameP_X,
                                                   std::string nameP_Y,
                                                   std::string nameP_Z ){

    if (verbose_) {
        xmf_out_.seekp (closingLinesPosition_);
        //        writeGeometry (M_xdmf);
        if (dim_ == 2) {
            xmf_out_ <<
            "      <Geometry Type=\"X_Y\">\n" <<
            "         <DataStructure Format=\"HDF\"\n" <<
            "                        Dimensions=\"" << nmbPointsGlob_ << "\"\n" <<
            "                        DataType=\"Float\"\n" <<
            "                        Precision=\"8\">\n" <<
            "             " << outputFilename_ << ":/" << nameP_X << "/Values\n" <<
            "         </DataStructure>\n" <<
            "         <DataStructure Format=\"HDF\"\n" <<
            "                        Dimensions=\"" << nmbPointsGlob_ << "\"\n" <<
            "                        DataType=\"Float\"\n" <<
            "                        Precision=\"8\">\n" <<
            "             " << outputFilename_ << ":/" << nameP_Y << "/Values\n" <<
            "         </DataStructure>\n" <<
            "      </Geometry>\n" <<
            "\n";
        }
        else if(dim_ ==3){
            xmf_out_ <<
            "      <Geometry Type=\"X_Y_Z\">\n" <<
            "         <DataStructure Format=\"HDF\"\n" <<
            "                        Dimensions=\"" << nmbPointsGlob_ << "\"\n" <<
            "                        DataType=\"Float\"\n" <<
            "                        Precision=\"8\">\n" <<
            "             " << outputFilename_ << ":/" << nameP_X << "/Values\n" <<
            "         </DataStructure>\n" <<
            "         <DataStructure Format=\"HDF\"\n" <<
            "                        Dimensions=\"" << nmbPointsGlob_ << "\"\n" <<
            "                        DataType=\"Float\"\n" <<
            "                        Precision=\"8\">\n" <<
            "             " << outputFilename_ << ":/" << nameP_Y << "/Values\n" <<
            "         </DataStructure>\n" <<
            "         <DataStructure Format=\"HDF\"\n" <<
            "                        Dimensions=\"" << nmbPointsGlob_ << "\"\n" <<
            "                        DataType=\"Float\"\n" <<
            "                        Precision=\"8\">\n" <<
            "             " << outputFilename_ << ":/" << nameP_Z << "/Values\n" <<
            "         </DataStructure>\n" <<
            "      </Geometry>\n" <<
            "\n";
        }
        closingLinesPosition_ = xmf_out_.tellp();

    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeXmfVariables(){

    std::string centerString;
    if (FEType_=="P0") {
        centerString = "Cell";
    }
    else if (FEType_=="P1-disc") {

    }
    else if (FEType_=="P1") {
        centerString = "Node";
    }

    else if (FEType_=="P2"){
        centerString = "Node";
    }
    else if (FEType_=="Q2"){
        centerString = "Node";
    }

    if(verbose_){
        xmf_out_.seekp (closingLinesPosition_);

        for (int i=0; i<varNames_.size(); i++) {
            int dof = varDofPerNode_.at(i);
            if (dof == 2) {
                dof=3; // ParaView always uses 3D Data. for 2D Data the last entries (for the 3rd Dim) are all zero.
            }
            xmf_out_ <<
            "\n      <Attribute\n" <<
            "         Type=\"" << varTypes_.at(i) << "\"\n" <<
            "         Center=\"" << centerString << "\"\n" <<
            "         Name=\"" << varNames_.at(i) << "\">\n";


            xmf_out_ <<
            "         <DataStructure ItemType=\"HyperSlab\"\n" <<
            "                        Dimensions=\"" << nmbExportValuesGlob_ << " " << dof << "\"\n" <<
            "                        Type=\"HyperSlab\">\n" <<
            "           <DataStructure  Dimensions=\"3 2\"\n" <<
            "                           Format=\"XML\">\n" <<
            "               0    0\n" <<
            "               1    1\n" <<
            "               " << nmbExportValuesGlob_ << " " << dof << "\n" <<
            "           </DataStructure>\n" <<

            "           <DataStructure  Format=\"HDF\"\n" <<
            "                           Dimensions=\"" << nmbExportValuesGlob_ << " " << dof << "\"\n" <<
            "                           DataType=\"Float\"\n" <<
            "                           Precision=\"8\">\n" <<
            "               " << outputFilename_ << ":/" << varNames_.at(i)
            << postfix_  << "/Values\n" << // see also in writeVector/scalar
            "           </DataStructure>\n" <<
            "         </DataStructure>\n";


            xmf_out_ <<
            "      </Attribute>\n";

            }

        xmf_out_ << "\n"
        "    </Grid>\n\n";
        closingLinesPosition_ = xmf_out_.tellp();
        // write closing lines
        xmf_out_ << closingLines_;
        xmf_out_.flush();
    }
}
template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::writeXmfTime(double time, double dt){

    if (comm_->getRank()==0) {

        xmf_times_out_.seekp (closingLinesPositionTimes_);

        xmf_times_out_ <<
        "<!-- Time " << time << " Iteration " << postfix_.substr (1, 5) << " -->\n" <<
        "    <Grid Name=\"Mesh Times " << time << "\">\n" <<
        "      <Time TimeType=\"Single\" Value=\"" << time << "\" />\n";
        xmf_times_out_	  <<
        "      <Topology\n"
        << "         Type=\"Polyvertex\"\n"
        << "         NumberOfElements=\""
        << 2
        << "\"\n"
        << "         BaseOffset=\""
        << 0
        << "\">\n"
        << "         <DataStructure Format=\"XML\"\n"
        << "                        Dimensions=\"2\"\n"
        << "                        DataType=\"Int\"\n"
        << "                        Precision=\"8\">\n"
        << "                    0 1 \n"
        << "         </DataStructure>\n" << "      </Topology>\n";


        xmf_times_out_ <<
        "      <Geometry Type=\"X_Y\">\n" <<
        "         <DataStructure Format=\"XML\"\n" <<
        "                        Dimensions=\"" << 2 << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             " << "0.0 1.0 \n" <<
        "         </DataStructure>\n" <<
        "         <DataStructure Format=\"XML\"\n" <<
        "                        Dimensions=\"" << 2 << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             " << "0.0 1.0 \n" <<
        "         </DataStructure>\n" <<        "      </Geometry>\n" <<
        "\n";

        xmf_times_out_ <<
        "\n      <Attribute\n" <<
        "         Type=\"" << "Scalar" << "\"\n" <<
        "         Center=\"" << "Node" << "\"\n" <<
        "         Name=\"" << "t_dt" << "\">\n";


        xmf_times_out_ <<
        "        <DataStructure  Format=\"XML\"\n" <<
        "                        Dimensions=\"" << 2 << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "           " << time << " " << dt << "\n" <<
        "        </DataStructure>\n";

        xmf_times_out_ <<
        "      </Attribute>\n";

        xmf_times_out_ << "\n"
        "    </Grid>\n\n";
        closingLinesPositionTimes_ = xmf_times_out_.tellp();
        // write closing lines
        xmf_times_out_ << closingLines_;
        xmf_times_out_.flush();
    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::prepareVectorField(MultiVecConstPtr_Type &u,
                                                          MultiVectorPtr_Type &u_export,
                                                          int dof) const{

    TEUCHOS_TEST_FOR_EXCEPTION(u->getNumVectors()!=1, std::logic_error, "Can only export single vector");

    Teuchos::ArrayRCP<const SC> tmpData = u->getData(0);
    for (int i=0; i<(u->getLocalLength()/dof); i++) {
        for (int j=0; j < dof; j++) {
            u_export->getDataNonConst(j)[i] =  tmpData[ dof * i + j ];
        }
    }
}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::prepareScalar(MultiVecConstPtr_Type &u,
                                                  MultiVectorPtr_Type &u_export) const{

    TEUCHOS_TEST_FOR_EXCEPTION(u->getNumVectors()!=1, std::logic_error, "Can only export single vector");

    Teuchos::ArrayRCP<const SC> tmpData = u->getData(0);
    for (int i=0; i<tmpData.size(); i++) {
        u_export->getDataNonConst(0)[i] =  tmpData[  i];
    }

}

template<class SC,class LO,class GO,class NO>
void ExporterParaView<SC,LO,GO,NO>::makePostfix(){

    int postfixLength = 5;
    std::ostringstream index;
    index.fill ( '0' );

    if (timeIndex_ % saveTimestep_ == 0){
        index << std::setw (postfixLength) << ( timeIndex_ / saveTimestep_ );
        postfix_ = "." + index.str();
    }
}
}
#endif
