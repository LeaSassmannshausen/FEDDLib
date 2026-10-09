#ifndef ExporterTxt_hpp
#define ExporterTxt_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/core_config.h"
#include <fstream>
#include <limits>

/*!
 Declaration of ExporterTxt
 
 @brief  ExporterTxt
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
class ExporterTxt {
public:
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;
    
    bool verbose_;

    std::ofstream txt_out_;
    
    ExporterTxt();

    /** @brief Configure fresh, resumed, or archived text output.
     * @param parameters Optional simulation settings containing Exporter options.
     * @param timeColumn Whether the first output column already contains time.
     * Value-only logs keep their format and use a .txt.times timestamp sidecar.
     */
    void setup(std::string filename, CommConstPtr_Type comm, int targetRank=0,
               ParameterListPtr_Type parameters=Teuchos::null, bool timeColumn=false);

    /// Write a value-only row with its physical time recorded separately.
    void exportDataAtTime(double time, double data);
    
    void exportData(double data);

    void writeTxt(double data);

    void exportData(double data1, double data2);

    void exportData(std::string data1, double data2);

    void exportData(std::string data1, std::string data2);

    void writeTxt(double data1, double data2);
    
    void writeTxt(std::string data1, double data2);

    void writeTxt(std::string data1, std::string data2);

    void closeExporter();

    
private:
    bool timeColumn_ = false;
    bool trackTime_ = false;
    bool resumeOutput_ = false;
    double lastTime_ = -std::numeric_limits<double>::infinity();
    std::ofstream time_out_;
    bool acceptTime(double time);
    void recordTime(double time);
    };
}

#endif
