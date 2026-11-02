#ifndef HDF5TOOLBOX_DEF_hpp
#define HDF5TOOLBOX_DEF_hpp

#include "HDF5VectorInfo.hpp"

/*
  

*/

namespace FEDD {

namespace hdf5_detail {

/** @brief Collectively write a replicated scalar/string dataset using rank zero's value.
 * All ranks supply the same datatype and extent. Non-root ranks select no data,
 * but still participate in H5Dwrite. Return an error on every rank if preparation
 * or writing fails, allowing the caller to close its handles before throwing.
 */
inline int writeReplicatedDataset(hid_t dataset, hid_t datatype, hid_t space,
                                  const void* value, const Teuchos::Comm<int>& comm)
{
    checkpoint::H5Handle fileSpace(H5Scopy(space), H5Sclose);
    checkpoint::H5Handle memorySpace(H5Scopy(space), H5Sclose);
    checkpoint::H5Handle transfer(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
    int localError = dataset < 0 || datatype < 0 || fileSpace < 0 || memorySpace < 0 || transfer < 0;
    if (!localError && comm.getRank() != 0) {
        if (H5Sselect_none(fileSpace) < 0 || H5Sselect_none(memorySpace) < 0)
            localError = 1;
    }
#ifdef HAVE_MPI
    if (transfer >= 0 && H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE) < 0)
        localError = 1;
#endif
    int globalError = 0;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &localError, &globalError);
    if (globalError)
        return globalError;
    localError = H5Dwrite(dataset, datatype, memorySpace, fileSpace, transfer, value) < 0;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &localError, &globalError);
    return globalError;
}

} // namespace hdf5_detail

struct FindDataset_t
{
  std::string name;
  bool found;
};

static herr_t FindDataset(hid_t loc_id, const char *name, void *opdata)
{
  std::string& token = ((FindDataset_t*)opdata)->name;
  if (token == name)
    ((FindDataset_t*)opdata)->found = true;

  return(0);
}

template <class SC, class LO, class GO, class NO>
HDF5Toolbox<SC, LO, GO, NO>::HDF5Toolbox(CommConstPtr_Type comm,
                                      Teuchos::RCP<Teuchos::ParameterList> parameters):
  isOpen_(false),
  comm_(comm),
  parameters_(parameters)
{}

// ----------------------------------------------------------------
/// @brief Write/export a vector X to the HDF5 file under the group name GroupName
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const MultiVectorConstPtr_Type X,
                                      bool writeTranspose, bool flushAfterWrite)
{

  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  hid_t       group_id, dset_id;
  hid_t       filespace_id, memspace_id;
  herr_t      status;

  // need a linear distribution to use hyperslabs
  MultiVectorPtr_Type linearX;

  {
    // Reuse the target buffer/importer while the source distribution and column count match.
    checkpoint::Timer timer(parameters_, "hdf5_redistribution");
    const MapConstPtr_Type sourceMap = X->getMap();
    // Tpetra maps are distributed objects. Reject partial replacement before
    // entering importer communication with inconsistent map-directory lifetimes.
    const int localRebuild = linearWriteBuffer_.is_null() || writeSourceMap_.is_null() ||
        linearWriteBuffer_->getNumVectors() != X->getNumVectors() ||
        writeSourceMap_->getTpetraMap().getRawPtr() != X->getMapTpetra().getRawPtr();
    int rebuild = 0;
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &localRebuild, &rebuild);
    TEUCHOS_TEST_FOR_EXCEPTION(rebuild != 0 && rebuild != comm_->getSize(), std::logic_error,
        "HDF5 write maps and column counts must be replaced collectively on all ranks");
    checkpoint::Timer setupTimer(parameters_, "hdf5_buffer_setup", rebuild != 0);
    if (rebuild) {
      MapPtr_Type linearMap = Teuchos::rcp(new Map_Type(sourceMap->getGlobalNumElements(),
          sourceMap->getNodeNumElements(), sourceMap->getIndexBase(), sourceMap->getComm()));
      linearWriteBuffer_ = Teuchos::rcp(new MultiVector_Type(linearMap, X->getNumVectors()));
      // Retain the immutable underlying map even if the caller replaces its wrapper's contents.
      writeSourceMap_ = Teuchos::rcp(new Map_Type(X->getMapTpetra()));
    }
    linearWriteBuffer_->importFromVector(X, true);
    linearX = linearWriteBuffer_;
  }


  int NumVectors = X->getNumVectors();
  int GlobalLength = X->getMap()->getGlobalNumElements(); //X.GlobalLength();

  // Whether or not we do writeTranspose or not is
  // handled by one of the components of q_dimsf, offset and count.
  // They are determined by indexT
  int indexT(0);
  if (writeTranspose) indexT = 1; 

  hsize_t q_dimsf[] = {static_cast<hsize_t>(GlobalLength), static_cast<hsize_t>(GlobalLength)};
  q_dimsf[indexT] = NumVectors;

  int localWriteError = 0;
  {
    checkpoint::Timer timer(parameters_, "hdf5_dataset_create");
    filespace_id = H5Screate_simple(2, q_dimsf, NULL);

    if (!isContained(GroupName))
      createGroup(GroupName);

    group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

    // Create the dataset with default properties and close filespace_id.
    dset_id = H5Dcreate(group_id, "Values", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    localWriteError = filespace_id < 0 || group_id < 0 || dset_id < 0;
    H5Sclose(filespace_id);
  }

  // Reject a failed collective dataset creation before writing vector properties.
  // In particular, an existing Values dataset must remain a vector-write error.
  int globalWriteError = 0;
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &localWriteError, &globalWriteError);
  if (globalWriteError) {
    if (dset_id >= 0) H5Dclose(dset_id);
    if (group_id >= 0) H5Gclose(group_id);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(globalWriteError, std::runtime_error,
                             "HDF5 vector write failed for group " << GroupName);

  // Create property list for collective dataset write.
  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif


  // Select hyperslab in the file. A rank without entries (e.g. a mesh with fewer subdomains than ranks)
  // has no first global element; it still takes part in the collective write, with an empty selection
  const bool hasLocalEntries = linearX->getLocalLength() > 0;
  const hsize_t firstIndex = hasLocalEntries ? static_cast<hsize_t>(linearX->getMap()->getGlobalElement(0)-X->getMap()->getIndexBase()) : 0;
  hsize_t offset[] = {firstIndex, firstIndex};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {static_cast<hsize_t>(linearX->getLocalLength()),
                     static_cast<hsize_t>(linearX->getLocalLength())};
  hsize_t block[] = {1, 1};
  
  for (int n = 0; n < NumVectors; ++n)
  {
    // Select hyperslab in the file
    offset[indexT] = n;
    count [indexT] = 1;

    // Print local info before selecting hyperslab
    // std::cout << "[Rank " << comm_->getRank() << "] "
    //           << "Writing vector " << n
    //           << " | file offset = (" << offset[0] << ", " << offset[1] << ")"
    //           << " | file count = (" << count[0] << ", " << count[1] << ")"
    //           << " | local length = " << linearX->getLocalLength()
    //           << " | global length = " << 1
    //           << std::endl;

    // Get the dataset's file dataspace and select the correct hyperslab
    filespace_id = H5Dget_space(dset_id);

    herr_t selStatus = hasLocalEntries
        ? H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, stride, count, block)
        : H5Sselect_none(filespace_id);

    if (selStatus < 0) {
      localWriteError = 1;
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR selecting hyperslab: "
                << "offset=(" << offset[0] << "," << offset[1] << "), "
                << "count=("  << count[0]  << "," << count[1]  << ")\n";
    }

    // Each process defines dataset in memory (its local buffer)
    hsize_t dimsm[] = { std::max<hsize_t>(linearX->getLocalLength(), 1) };
    memspace_id = H5Screate_simple(1, dimsm, NULL);
    if (!hasLocalEntries)
      H5Sselect_none(memspace_id);

    if (memspace_id < 0) {
      localWriteError = 1;
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR creating memory dataspace of size "
                << dimsm[0] << std::endl;
    }

    // Print detailed info before writing
    // std::cout << "[Rank " << comm_->getRank() << "] "
    //           << "Calling H5Dwrite for vector " << n
    //           << " | memspace dims = " << dimsm[0]
    //           << " | plist_id = " << plist_id_
    //           << std::endl;

    // Measure the vector-data transfer separately from dataset preparation.
    herr_t writeStatus;
    {
      checkpoint::Timer timer(parameters_, "hdf5_dataset_write");
      writeStatus = H5Dwrite(
        dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id,
        plist_id_, linearX->getData(n).get());
    }

    if (writeStatus < 0) {
      localWriteError = 1;
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR during H5Dwrite for vector " << n << std::endl;
    }


    // Optional: flush stdout to preserve rank ordering in logs
    std::cout.flush();
    std::cerr.flush();
  }

  hid_t filespace = H5Dget_space(dset_id);
  hsize_t totalElems = H5Sget_simple_extent_npoints(filespace);

  // std::cout << "[Rank " << comm_->getRank()
  //         << "] Dataset total elements = " << totalElems << std::endl;

  hssize_t fileCount = H5Sget_select_npoints(filespace_id);
  hssize_t memCount  = H5Sget_select_npoints(memspace_id);

  // std::cout << "[Rank " << comm_->getRank() << "] "
  //         << "Hyperslab file selection = " << fileCount
  //         << ", memory selection = " << memCount << std::endl;
          
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id_);
  // Close local memory dataspace
  H5Sclose(memspace_id);
  
  {
    // These HDF5 records describe a vector; the checkpoint manifest is timed elsewhere.
    checkpoint::Timer timer(parameters_, "hdf5_vector_metadata_write");
    write(GroupName, "GlobalLength", GlobalLength);
    write(GroupName, "NumVectors", NumVectors);
    write(GroupName, "__type__", "Tpetra_MultiVector");
  }
  // Do not let checkpoint publication treat a logged vector-write error as
  // success. All healthy ranks finish the write sequence before failing together.
  if (flushAfterWrite) {
    // Standalone writes flush immediately; checkpoint batches flush at completion.
    checkpoint::Timer timer(parameters_, "hdf5_flush_after_write");
    if (H5Fflush(file_id_, H5F_SCOPE_GLOBAL) < 0) localWriteError = 1;
  }
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &localWriteError, &globalWriteError);
  TEUCHOS_TEST_FOR_EXCEPTION(globalWriteError != 0, std::runtime_error,
                             "HDF5 vector write failed for group " << GroupName);
}
// ---------------------------------------------------------
// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const MapConstPtr_Type Map,
                        MultiVectorPtr_Type X)
{
  TEUCHOS_TEST_FOR_EXCEPTION(X->getNumVectors() != 1, std::logic_error,
                            "HDF5 vector import currently supports one vector only.");
  checkpoint::validateVectorDataset(file_id_, GroupName, Map->getGlobalNumElements(), X->getNumVectors());
  // gets the length of the std::vector
  int GlobalLength;
  readIntVectorProperties(GroupName, GlobalLength);

  MapPtr_Type linearMap = Teuchos::rcp(new Map_Type(X->getMap()->getGlobalNumElements(),X->getMap()->getNodeNumElements(), X->getMap()->getIndexBase(), X->getMap()->getComm()));

  // we need to first create a linear map, read the std::vector,
  // then import it to the actual nonlinear map
  MultiVectorPtr_Type linearX = Teuchos::rcp(new MultiVector_Type(linearMap, X->getNumVectors()));

  read(GroupName, "Values", linearMap->getNodeNumElements(), linearMap->getGlobalNumElements(),
        H5T_NATIVE_DOUBLE, linearX->getDataNonConst(0).get());

  X->importFromVector(linearX);

}

// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::readIntVectorProperties(const std::string& GroupName,
                                           int& GlobalLength)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " <<  GroupName << " not found");

  std::string Label;
  read(GroupName, "__type__", Label);
                    
  read(GroupName, "GlobalLength", GlobalLength);
}

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const std::string& DataSetName,
                        GO MySize, int GlobalSize,
                        const hid_t type, void* data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  // Convert MySize to hsize_t for HDF5 calls
hsize_t MySize_t = MySize;

// Compute the prefix sum across ranks to get the global offset
GO itmp = 0;
tpetraScanSum(comm_, &MySize, &itmp, 1);

// Exclusive prefix (start index for this rank)
hsize_t Offset_t = static_cast<hsize_t>(itmp - MySize);

// === DEBUG OUTPUT ===
// std::cout << "[Rank " << comm_->getRank() << "] "
//           << "Preparing to read dataset '" << DataSetName << "' in group '" << GroupName << "'\n"
//           << "  local size (MySize)     = " << MySize << "\n"
//           << "  computed offset (Offset)= " << Offset_t << "\n"
//           << "  total elements (sum)    = " << itmp << std::endl;

// Open group and dataset
hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
if (group_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Gopen failed for group '" 
            << GroupName << "'\n";
}

hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);
if (dataset_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Dopen failed for dataset '" 
            << DataSetName << "'\n";
}

// Get the dataset’s file dataspace
hid_t filespace_id = H5Dget_space(dataset_id);
if (filespace_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Dget_space failed\n";
}

// Select hyperslab in the file for this process
hsize_t offset[2] = { 0, Offset_t };   // row index, element offset
hsize_t count [2] = { 1, MySize_t };  // one vector, 29 elements

// A rank without entries selects nothing
herr_t selStatus = MySize_t > 0
    ? H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL)
    : H5Sselect_none(filespace_id);

if (selStatus < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR selecting hyperslab: offset=" << Offset_t
            << " count=" << MySize_t << std::endl;
} 
// else {
//   std::cout << "[Rank " << comm_->getRank()
//             << "] Selected hyperslab: offset=" << Offset_t
//             << " count=" << MySize_t << std::endl;
// }

// Create memory dataspace for local buffer
hsize_t memSize_t = std::max<hsize_t>(MySize_t, 1);
hid_t mem_dataspace = H5Screate_simple(1, &memSize_t, NULL);
if (MySize_t == 0)
  H5Sselect_none(mem_dataspace);
if (mem_dataspace < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR creating memory dataspace (size=" << MySize_t << ")\n";
}

// Print before reading
// std::cout << "[Rank " << comm_->getRank()
//           << "] Reading from file hyperslab offset=" << Offset_t
//           << " count=" << MySize_t
//           << " into local buffer at " << static_cast<void*>(data)
//           << std::endl;

// Perform the read
herr_t status = H5Dread(dataset_id, type, mem_dataspace, filespace_id,
                        H5P_DEFAULT, data);

if (status < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR during H5Dread (offset=" << Offset_t
            << ", count=" << MySize_t << ")\n";
} 
// else {
//   std::cout << "[Rank " << comm_->getRank()
//             << "] Successfully read " << MySize_t << " elements from dataset '"
//             << DataSetName << "'\n";
// }

hid_t filespace = H5Dget_space(dataset_id);
hsize_t totalElems = H5Sget_simple_extent_npoints(filespace);

// std::cout << "[Rank " << comm_->getRank()
//         << "] Dataset total elements = " << totalElems << std::endl;


hssize_t fileCount = H5Sget_select_npoints(filespace_id);
hssize_t memCount  = H5Sget_select_npoints(mem_dataspace);

// std::cout << "[Rank " << comm_->getRank() << "] "
//           << "Hyperslab file selection = " << fileCount
//           << ", memory selection = " << memCount << std::endl;

// Close resources
H5Sclose(mem_dataspace);
H5Sclose(filespace_id);
H5Dclose(dataset_id);
H5Gclose(group_id);

// Optional: ensure ordered output by rank
comm_->barrier();
std::cout.flush();
std::cerr.flush();
//  H5Dclose(filespace_id);
}

// -------------------------------------------------------------------------
/// @brief Create a new file.
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::create(const std::string FileName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isOpen(),std::runtime_error,"HDF5Toolbox:: an HDF5 is already open, first close the current one - using method Close(), then open/create a new one");

  checkpoint::Timer timer(parameters_, "hdf5_file_create");
  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);
// #ifdef HAVE_MPI
  {
    // Tell HDF5 what MPI communicator to use for parallel file access
    // for the above property list.
    //
    // HAVE_MPI is defined, so we know that Trilinos was built with
    // MPI.  However, we don't know whether Comm_ wraps an MPI
    // communicator.  Comm_ could very well be a serial communicator.
    MPI_Comm mpiComm = MPI_COMM_NULL; 

    Teuchos::RCP<const Teuchos::MpiComm<int>> mpiWrapper =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm_, false);

    if (!mpiWrapper.is_null()) {
      mpiComm = *(mpiWrapper->getRawMpiComm());
    }
    else {
      // Try Serial communicator next
      Teuchos::RCP<const Teuchos::SerialComm<int>> serialWrapper =
        Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<int>>(comm_, false);

      if (!serialWrapper.is_null()) {
        // Comm_ is a Teuchos::SerialComm<int>
        // Use MPI_COMM_SELF for serial HDF5 access
        mpiComm = MPI_COMM_SELF;
      }
      else {
        // Unknown communicator subclass
        const char* const errMsg =
          "Tpetra::HDF5::Create: This HDF5 object was created with a "
          "Teuchos::Comm<int> that is neither Teuchos::MpiComm<int> "
          "nor Teuchos::SerialComm<int>. We don't know how to get an "
          "MPI_Comm from it.";
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,errMsg);
      }
    }

    // By this point, mpiComm should be something other than
    // MPI_COMM_NULL.  Otherwise, Comm_ wraps MPI_COMM_NULL.
    if (mpiComm == MPI_COMM_NULL) {
      const char* const errMsg = "TpetraExt::HDF5::Create: The Tpetra_Comm "
        "object with which this HDF5 instance was created wraps MPI_COMM_NULL, "
        "which is an invalid MPI communicator.  HDF5 requires a valid MPI "
        "communicator.";
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,errMsg);
    }

    // Tell HDF5 what MPI communicator to use for parallel file access
    // for the above property list.  For details, see e.g.,
    //
    // http://www.hdfgroup.org/HDF5/doc/UG/08_TheFile.html
    //
    // [last accessed 06 Oct 2011]
    H5Pset_fapl_mpio(plist_id_, mpiComm, MPI_INFO_NULL);

    // Aggregate HDF5's internal metadata-cache writes during collective flushes.
    // Application property datasets are handled separately by writeReplicatedDataset.
    const int localMetadataError = H5Pset_coll_metadata_write(plist_id_, true) < 0;
    int globalMetadataError = 0;
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &localMetadataError, &globalMetadataError);
    if (globalMetadataError)
      H5Pclose(plist_id_);
    TEUCHOS_TEST_FOR_EXCEPTION(globalMetadataError, std::runtime_error,
                               "Cannot enable collective HDF5 metadata writes for " << FileName);
  }
// #endif

#if 0
  unsigned int boh = H5Z_FILTER_MAX;
  H5Pset_filter(plist_id_, H5Z_FILTER_DEFLATE, H5Z_FILTER_MAX, 0, &boh);
#endif

  // create the file collectively and release property list identifier.
  file_id_ = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      plist_id_);
  H5Pclose(plist_id_);

  isOpen_ = true;
}


// -------------------------------------------------------------------------
/// @brief Checking if the dataset Name is contained in the group GroupName

template <class SC, class LO, class GO, class NO>
bool HDF5Toolbox<SC, LO, GO, NO>::isContained(std::string Name, std::string GroupName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  FindDataset_t data;
  data.name = Name;
  data.found = false;

  // recursively look for groups
  size_t pos = Name.find("/");
  if (pos != std::string::npos)
  {
    std::string NewGroupName = Name.substr(0, pos);
    if (GroupName != "")
      NewGroupName = GroupName + "/" + NewGroupName;
    std::string NewName = Name.substr(pos + 1);
    return isContained(NewName, NewGroupName);
  }

  GroupName = "/" + GroupName;

  //int idx_f =
  H5Giterate(file_id_, GroupName.c_str(), NULL, FindDataset, (void*)&data);

  return(data.found);
}
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
//! Create group \c GroupName.
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::createGroup(const std::string& GroupName)
{
  hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}
// -------------------------------------------------------------------------



// ==========================================================================
// ==========================================================================

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const std::string& DataSetName,
                            double what)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_DOUBLE,
                            filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  const int writeError = hdf5_detail::writeReplicatedDataset(
      dset_id, H5T_NATIVE_DOUBLE, filespace_id, &what, *comm_);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  TEUCHOS_TEST_FOR_EXCEPTION(writeError, std::runtime_error,
                             "HDF5 property write failed for " << GroupName << "/" << DataSetName);
}

// ==========================================================================

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const std::string& DataSetName,
                            int what)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_INT,
                      filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  const int writeError = hdf5_detail::writeReplicatedDataset(
      dset_id, H5T_NATIVE_INT, filespace_id, &what, *comm_);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  TEUCHOS_TEST_FOR_EXCEPTION(writeError, std::runtime_error,
                             "HDF5 property write failed for " << GroupName << "/" << DataSetName);
}


// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName,
                            const std::string& DataSetName,
                            const std::string& data)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hsize_t len = 1;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataspace_id = H5Screate_simple(1, &len, NULL);

  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, data.size() + 1);

  hid_t dataset_id = H5Dcreate(group_id, DataSetName.c_str(), atype,
                               dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  const int writeError = hdf5_detail::writeReplicatedDataset(
      dataset_id, atype, dataspace_id, data.c_str(), *comm_);
  
  // Close/release resources.
  H5Dclose(dataset_id);
  H5Gclose(group_id);
  H5Sclose(dataspace_id);
  H5Tclose(atype);
  TEUCHOS_TEST_FOR_EXCEPTION(writeError, std::runtime_error,
                             "HDF5 property write failed for " << GroupName << "/" << DataSetName);
}
// ==========================================================================
// ==========================================================================

// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const std::string& DataSetName, int& data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " + GroupName + " not found");

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName,
                           const std::string& DataSetName,
                           std::string& data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " + GroupName + " not found");

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  hid_t datatype_id = H5Dget_type(dataset_id);
//  size_t typesize_id = H5Tget_size(datatype_id);
  H5T_class_t typeclass_id = H5Tget_class(datatype_id);

  if(typeclass_id != H5T_STRING)
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error, "requested group " + GroupName + " is not a std::string");

  char data2[160];
  H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, data2) ;                

  data = data2;

  H5Dclose(dataset_id);
  H5Gclose(group_id);
}
// ------------------------------------------
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::open(const std::string FileName, int AccessType)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiWrapper =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm_, false);

  if (!mpiWrapper.is_null()) {
    mpiComm = *(mpiWrapper->getRawMpiComm());
  }

  H5Pset_fapl_mpio(plist_id_, mpiComm, MPI_INFO_NULL);
#endif

  // create the file collectively and release property list identifier.
  file_id_ = H5Fopen(FileName.c_str(), AccessType, plist_id_);
  H5Pclose(plist_id_);

  isOpen_ = true;
}


// ----------------------------------------------
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::tpetraScanSum(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                   const GO* sendbuf,
                   GO* recvbuf,
                   int count)
{
#ifdef HAVE_MPI
  // Try to downcast to MPI communicator
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm, false);

  if (!mpiComm.is_null()) {
    // Use MPI_Scan directly
    MPI_Comm rawComm = *(mpiComm->getRawMpiComm());
    MPI_Scan(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, rawComm);
    return;
  }
#endif

  // Serial (or non-MPI) fallback
  for (int i = 0; i < count; ++i) {
    recvbuf[i] = sendbuf[i];
  }
}
}
#endif
