#ifndef HDF5EXPORT_DEF_hpp
#define HDF5EXPORT_DEF_hpp

namespace FEDD {
 
template<class SC,class LO,class GO,class NO>
HDF5Export<SC,LO,GO,NO>::HDF5Export(MapConstPtr_Type writeMap, std::string outputFilename,
                                Teuchos::RCP<Teuchos::ParameterList> parameters):
hdf5exporter_(),
comm_(),
pendingWrites_(false)
{

    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( writeMap->getComm() );

    hdf5exporter_.reset( new HDF5_Type(mpiComm, parameters) ); // Building HDF5 Exporter

    // The output file may lie in a (checkpoint) directory: create it first
    std::string::size_type slash = outputFilename.rfind('/');
    if (slash != std::string::npos) {
        if (mpiComm->getRank() == 0)
            mkdir(outputFilename.substr(0, slash).c_str(), 0777);
        mpiComm->barrier();
    }

    hdf5exporter_->create(outputFilename+".h5"); // Creating output file with the 'outoutFilename'

    outputFilename_ = outputFilename; 
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::writeVariablesHDF5(std::string varName,
    const MultiVectorConstPtr_Type writeVector, bool flushAfterWrite){

    hdf5exporter_->write(varName, writeVector, false, flushAfterWrite);
    pendingWrites_ = !flushAfterWrite;
    
    if(writeVector->getMap()->getComm()->getRank() == 0 )
        std::cout << " HDF5_Export:: Exporting to file " << outputFilename_ << " with variable name " << varName << std::endl;

}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::flush(){
    hdf5exporter_->flush();
    pendingWrites_ = false;
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::closeExporter(){
    if (pendingWrites_)
        flush();
    hdf5exporter_->close();
}


}
#endif
