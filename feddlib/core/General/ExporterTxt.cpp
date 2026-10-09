#include "ExporterTxt.hpp"
#include "OutputHistory.hpp"
#include <Teuchos_Array.hpp>
#include <iomanip>
/*!
 Definition of ExporterTxt
 
 @brief  ExporterTxt
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {
ExporterTxt::ExporterTxt():
    verbose_(false),
    txt_out_()
{
    
}


void ExporterTxt::setup(std::string filename, CommConstPtr_Type comm, int targetRank,
                        ParameterListPtr_Type parameters, bool timeColumn){
    verbose_ = comm->getRank() == targetRank;
    timeColumn_ = timeColumn;
    trackTime_ = !parameters.is_null();
    resumeOutput_ = output::resume(parameters);
    const std::string dataFile = filename + ".txt", timeFile = dataFile + ".times";
    output::archive(parameters, *comm, {dataFile, timeFile});
    // Cache the original shared clock before opening/pruning time.txt. This also
    // permits resuming legacy value-only logs with a matching time.txt series.
    Teuchos::Array<double> legacyTimes;
    if (resumeOutput_) {
        auto& options = parameters->sublist("Exporter");
        if (!options.isParameter("_Legacy text times")) {
            output::onRank(*comm, 0, [&] {
                std::string clock = "time" + parameters->sublist("General").get("Export Suffix", std::string("")) + ".txt";
                if (!std::filesystem::exists(clock)) clock = "time.txt";
                if (std::filesystem::exists(clock)) {
                    std::istringstream input(output::read(clock));
                    double time;
                    while (input >> time) legacyTimes.push_back(time);
                    if (!input.eof()) throw std::runtime_error("Invalid legacy time.txt series");
                }
            });
            int size = static_cast<int>(legacyTimes.size());
            Teuchos::broadcast(*comm, 0, 1, &size);
            legacyTimes.resize(size);
            if (size) Teuchos::broadcast(*comm, 0, size, legacyTimes.getRawPtr());
            options.set("_Legacy text times", legacyTimes);
        } else legacyTimes = options.get<Teuchos::Array<double>>("_Legacy text times");
    }
    output::onRank(*comm, targetRank, [&] {
        bool append = false;
        if (resumeOutput_ && std::filesystem::exists(dataFile)) {
            std::vector<std::string> rows;
            std::istringstream input(output::read(dataFile));
            std::string row;
            while (std::getline(input, row)) if (!row.empty()) rows.push_back(row);
            std::vector<double> times;
            if (std::filesystem::exists(timeFile)) {
                std::istringstream inputTimes(output::read(timeFile));
                double time;
                while (inputTimes >> time) times.push_back(time);
                if (!inputTimes.eof()) throw std::runtime_error("Invalid text timestamps: " + timeFile);
                if (times.size() > rows.size() || rows.size() > times.size() + 1)
                    throw std::runtime_error("Text/timestamp row counts disagree: " + dataFile);
                // A crash may leave one value written before its timestamp.
                rows.resize(times.size());
            } else if (timeColumn_) {
                for (const auto& line : rows) {
                    std::istringstream fields(line);
                    double time;
                    if (!(fields >> time)) throw std::runtime_error("Invalid timestamp in " + dataFile);
                    times.push_back(time);
                }
            } else {
                if (legacyTimes.size() != rows.size() && legacyTimes.size() != rows.size() + 1)
                    throw std::runtime_error("Cannot resume value-only log without matching timestamps: " + dataFile);
                // Iteration logs omit the initial row present in some time.txt files.
                times.assign(legacyTimes.end() - rows.size(), legacyTimes.end());
            }
            std::ostringstream retained, retainedTimes;
            retainedTimes << std::setprecision(17);
            const double restart = output::restartTime(parameters);
            double previous = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < times.size(); ++i) {
                if (!std::isfinite(times[i]) || times[i] <= previous)
                    throw std::runtime_error("Text timestamps must increase: " + dataFile);
                previous = times[i];
                if (times[i] <= restart + 1.e-12) {
                    retained << rows[i] << '\n';
                    retainedTimes << times[i] << '\n';
                    lastTime_ = times[i];
                }
            }
            output::replace(dataFile, retained.str());
            output::replace(timeFile, retainedTimes.str());
            append = true;
        } else if (resumeOutput_ && std::filesystem::exists(timeFile)) {
            throw std::runtime_error("Missing text data for timestamps: " + dataFile);
        }
        txt_out_.open(dataFile, append ? ios_base::app : ios_base::out);
        if (!txt_out_) throw std::runtime_error("Cannot open text output " + dataFile);
        if (trackTime_) {
            txt_out_ << std::setprecision(17);
            time_out_.open(timeFile, append ? ios_base::app : ios_base::out);
            time_out_ << std::setprecision(17);
            if (!time_out_) throw std::runtime_error("Cannot open text timestamps " + timeFile);
        }
    });
}

bool ExporterTxt::acceptTime(double time)
{
    if (!std::isfinite(time)) throw std::runtime_error("Nonfinite output timestamp");
    return !resumeOutput_ || time > lastTime_ + 1.e-12;
}

void ExporterTxt::recordTime(double time)
{
    if (trackTime_) {
        time_out_ << time << '\n';
        time_out_.flush();
        if (!txt_out_ || !time_out_) throw std::runtime_error("Cannot write text output history");
    }
    lastTime_ = time;
}

void ExporterTxt::exportDataAtTime(double time, double data)
{
    if (verbose_ && acceptTime(time)) {
        writeTxt(data);
        recordTime(time);
    }
}
    
void ExporterTxt::exportData(double data){
   
    writeTxt(data);
    
}
void ExporterTxt::exportData(std::string data1, double data2){
   
    writeTxt(data1,data2);
    
}

void ExporterTxt::exportData(std::string data1, std::string data2){
   
    writeTxt(data1,data2);
    
}
void ExporterTxt::exportData(double data1, double data2){
   
    writeTxt(data1,data2);
    
}

void ExporterTxt::writeTxt(double data){
    if (verbose_) {
        txt_out_ << data << "\n";
        txt_out_.flush();
    }

}


void ExporterTxt::writeTxt(double data1, double data2){
    if (verbose_) {
        if (timeColumn_ && !acceptTime(data1)) return;
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
        if (timeColumn_) recordTime(data1);
    }

}

void ExporterTxt::writeTxt(std::string data1, double data2){
    if (verbose_) {
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::writeTxt(std::string data1, std::string data2){
    if (verbose_) {
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::closeExporter(){
    if (verbose_) {
        txt_out_.close();
        if (time_out_.is_open()) time_out_.close();
    }

}
}

