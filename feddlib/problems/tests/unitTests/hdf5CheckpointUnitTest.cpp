#include "feddlib/core/Checkpointing/CheckpointMetadata.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include <Tpetra_Core.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Array.hpp>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

using namespace FEDD;

namespace {
int timerCalls(const char* operation)
{
    Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter(
        std::string("FEDDlib: Checkpoint Feature - ") + operation);
    return timer.is_null() ? 0 : timer->numCalls();
}

Teuchos::RCP<const Map<>> cyclicMap(int length, int shift, bool reverse, int indexBase,
                                   const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
{
    Teuchos::Array<default_go> entries;
    for (int i = 0; i < length; ++i)
        if ((i + shift) % comm->getSize() == comm->getRank())
            entries.push_back(indexBase + i);
    if (reverse)
        std::reverse(entries.begin(), entries.end());
    return Teuchos::rcp(new Map<>(length, entries(), indexBase, comm));
}

Teuchos::RCP<MultiVector<>> fieldValues(const Teuchos::RCP<const Map<>>& map, int columns, double offset)
{
    Teuchos::RCP<MultiVector<>> values = Teuchos::rcp(new MultiVector<>(map, columns));
    for (int column = 0; column < columns; ++column) {
        Teuchos::ArrayRCP<default_sc> data = values->getDataNonConst(column);
        for (int i = 0; i < data.size(); ++i)
            data[i] = offset + map->getGlobalElement(i) + 100. * column;
    }
    return values;
}

struct ExpectedDataset {
    const char* key;
    int length;
    int columns;
    int indexBase;
    double offset;
};
}

// Exercise real parallel HDF5 I/O without a PDE solve or a reference-file dependency.
int main(int argc, char** argv)
{
    Tpetra::ScopeGuard session(&argc, &argv);
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    const std::string directory = "hdf5CheckpointUnitTest-data";
    checkpoint::onRoot(*comm, [&] {
        std::filesystem::remove_all(directory);
        std::filesystem::create_directory(directory);
    });
    const std::string filename = directory + "/properties.h5";
    {
        HDF5Toolbox<> writer(comm);
        writer.create(filename);
        // Non-root payloads deliberately differ, with matching string extents.
        writer.write("properties", "integer", comm->getRank() == 0 ? 17 : -1);
        writer.write("properties", "real", comm->getRank() == 0 ? 2.5 : -1.);
        writer.write("properties", "string", std::string(comm->getRank() == 0 ? "root" : "xxxx"));
        writer.close();
    }
    checkpoint::onRoot(*comm, [&] {
        checkpoint::H5Handle file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        checkpoint::H5Handle integer(H5Dopen(file, "properties/integer", H5P_DEFAULT), H5Dclose);
        checkpoint::H5Handle real(H5Dopen(file, "properties/real", H5P_DEFAULT), H5Dclose);
        checkpoint::H5Handle string(H5Dopen(file, "properties/string", H5P_DEFAULT), H5Dclose);
        checkpoint::H5Handle stringType(H5Dget_type(string), H5Tclose);
        int integerValue = 0;
        double realValue = 0.;
        char stringValue[5] = {};
        if (H5Dread(integer, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &integerValue) < 0 ||
            H5Dread(real, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &realValue) < 0 ||
            H5Dread(string, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, stringValue) < 0 ||
            integerValue != 17 || realValue != 2.5 || std::string(stringValue) != "root")
            throw std::runtime_error("Replicated datasets did not preserve rank zero's values");
        std::cout << "PASS rank-zero integer, real and string property values\n";
    });

    // Rank zero writes the properties even when it owns no vector entries.
    Teuchos::RCP<const Map<>> writeMap = Teuchos::rcp(new Map<>(
        1, comm->getRank() == comm->getSize() - 1 ? 1 : 0, 0, comm));
    Teuchos::RCP<MultiVector<>> values = Teuchos::rcp(new MultiVector<>(writeMap));
    values->putScalar(7.25);
    {
        HDF5Export<> writer(writeMap, directory + "/empty-ranks");
        writer.writeVariablesHDF5("state", values);
        writer.closeExporter();
    }
    Teuchos::RCP<const Map<>> readMap = Teuchos::rcp(new Map<>(
        1, comm->getRank() == 0 ? 1 : 0, 0, comm));
    HDF5Import<> reader(readMap, directory + "/empty-ranks");
    Teuchos::RCP<MultiVector<>> restored = reader.readVariablesHDF5("state");
    int localPassed = restored->getLocalLength() == 0 || restored->getData(0)[0] == 7.25;
    int allPassed = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localPassed, &allPassed);
    if (!allPassed)
        throw std::runtime_error("Vector restoration failed after redistribution with empty ranks");
    reader.closeImporter();
    if (comm->getRank() == 0)
        std::cout << "PASS empty writer/reader ranks and changed distribution\n";

    const auto check = [&](bool passed, const char* label) {
        int local = passed ? 1 : 0, global = 0;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &local, &global);
        if (!global)
            throw std::runtime_error(label);
        if (comm->getRank() == 0)
            std::cout << "PASS " << label << '\n';
    };
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::rcp(new Teuchos::ParameterList);
    parameters->sublist("General").set("Checkpoint timings", true);
    const int automaticBefore = timerCalls("hdf5_flush_after_write");
    const int explicitBefore = timerCalls("hdf5_flush");
    {
        HDF5Export<> writer(writeMap, directory + "/batch", parameters);
        writer.writeVariablesHDF5("first", values, false);
        values->putScalar(8.25);
        writer.writeVariablesHDF5("second", values, false);
        check(timerCalls("hdf5_flush_after_write") == automaticBefore &&
              timerCalls("hdf5_flush") == explicitBefore, "batch defers both vector flushes");
        writer.flush();
        check(timerCalls("hdf5_flush") == explicitBefore + 1, "one checked flush completes the batch");
        values->putScalar(9.25);
        writer.writeVariablesHDF5("automatic", values);
        check(timerCalls("hdf5_flush_after_write") == automaticBefore + 1,
              "ordinary writes still flush immediately");
        writer.closeExporter();
        check(timerCalls("hdf5_flush") == explicitBefore + 1, "closure does not repeat a completed flush");
        bool rejected = false;
        try { writer.flush(); }
        catch (const std::logic_error&) { rejected = true; }
        check(rejected, "flush rejects a closed file on every rank");
    }
    {
        HDF5Import<> batchReader(readMap, directory + "/batch");
        for (const std::pair<std::string, double>& state : {
                 std::make_pair(std::string("first"), 7.25),
                 std::make_pair(std::string("second"), 8.25),
                 std::make_pair(std::string("automatic"), 9.25)}) {
            restored = batchReader.readVariablesHDF5(state.first);
            check(restored->getLocalLength() == 0 || restored->getData(0)[0] == state.second,
                  "batched and automatic states retain their own values");
        }
        batchReader.closeImporter();
    }
    const int beforeClose = timerCalls("hdf5_flush");
    {
        HDF5Export<> writer(writeMap, directory + "/close-batch", parameters);
        writer.writeVariablesHDF5("state", values, false);
        writer.closeExporter();
    }
    check(timerCalls("hdf5_flush") == beforeClose + 1, "explicit closure flushes pending recovery data");
    {
        HDF5Import<> closeReader(readMap, directory + "/close-batch");
        restored = closeReader.readVariablesHDF5("state");
        check(restored->getLocalLength() == 0 || restored->getData(0)[0] == 9.25,
              "deferred data survives explicit closure");
        closeReader.closeImporter();
    }
    const int setupBefore = timerCalls("hdf5_buffer_setup");
    const std::vector<ExpectedDataset> expected = {
        {"initial", 4, 1, 0, 10.}, {"reused", 4, 1, 0, 20.}, {"equivalent", 4, 1, 0, 25.},
        {"redistributed", 4, 1, 0, 30.}, {"columns", 4, 2, 0, 40.},
        {"length", 5, 1, 0, 50.}, {"index-base", 5, 1, 5, 60.},
        {"ordering", 5, 1, 5, 70.}};
    {
        HDF5Toolbox<> writer(comm, parameters);
        writer.create(directory + "/cache.h5");
        Teuchos::RCP<const Map<>> map = cyclicMap(4, 0, false, 0, comm);
        writer.write("initial", fieldValues(map, 1, 10.), false, false);
        Teuchos::RCP<const MultiVector<>> original = fieldValues(map, 1, 20.);
        Teuchos::RCP<MultiVector<>> copied = Teuchos::rcp(new MultiVector<>(original));
        writer.write("reused", copied, false, false);
        check(timerCalls("hdf5_buffer_setup") == setupBefore + 1,
              "copied solution reuses its buffer and importer");
        Teuchos::RCP<const Map<>> equivalent = cyclicMap(4, 0, false, 0, comm);
        if (comm->getRank() == 0)
            map = equivalent;
        bool partialRejected = false;
        try { writer.write("partial", fieldValues(map, 1, 25.), false, false); }
        catch (const std::logic_error&) { partialRejected = true; }
        check(partialRejected || comm->getSize() == 1,
              "partial-rank map replacement is rejected collectively before I/O");
        map = equivalent;
        writer.write("equivalent", fieldValues(map, 1, 25.), false, false);
        check(timerCalls("hdf5_buffer_setup") == setupBefore + 2,
              "collective map replacement rebuilds every importer");
        map = cyclicMap(4, 1, true, 0, comm);
        writer.write("redistributed", fieldValues(map, 1, 30.), false, false);
        check(timerCalls("hdf5_buffer_setup") == setupBefore + 3,
              "changed ownership rebuilds the importer");
        writer.write("columns", fieldValues(map, 2, 40.), false, false);
        map = cyclicMap(5, 0, false, 0, comm);
        writer.write("length", fieldValues(map, 1, 50.), false, false);
        map = cyclicMap(5, 0, false, 5, comm);
        writer.write("index-base", fieldValues(map, 1, 60.), false, false);
        map = cyclicMap(5, 0, true, 5, comm);
        writer.write("ordering", fieldValues(map, 1, 70.), false, false);
        check(timerCalls("hdf5_buffer_setup") == setupBefore + 7,
              "column count, global length, index base and local ordering invalidate the cache");
        writer.flush();
        writer.close();
    }
    checkpoint::onRoot(*comm, [&] {
        checkpoint::H5Handle file(H5Fopen((directory + "/cache.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        if (comm->getSize() > 1 && H5Lexists(file, "partial", H5P_DEFAULT) > 0)
            throw std::runtime_error("Rejected map replacement created an incomplete dataset");
        for (const ExpectedDataset& state : expected) {
            checkpoint::validateVectorDataset(file, state.key, state.length, state.columns);
            checkpoint::H5Handle dataset(H5Dopen(file,
                (std::string(state.key) + "/Values").c_str(), H5P_DEFAULT), H5Dclose);
            std::vector<double> data(state.length * state.columns);
            if (H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
                throw std::runtime_error("Cannot inspect cached-buffer output");
            for (int column = 0; column < state.columns; ++column)
                for (int i = 0; i < state.length; ++i)
                    if (data[column * state.length + i] != state.offset + state.indexBase + i + 100. * column)
                        throw std::runtime_error("Cached buffer retained incorrect values for " + std::string(state.key));
        }
        std::cout << "PASS cached and rebuilt buffers preserve values, columns and global ordering\n";
    });
    return 0;
}
