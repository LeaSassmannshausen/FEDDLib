#pragma once
#include <mpi.h>
#include <sys/resource.h>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <string>
   
namespace FEDD {

class MemoryLogger {
public:
    MemoryLogger(MPI_Comm comm, const std::string &basename = "memlog")
        : comm_(comm)
    {
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);

    }

    ~MemoryLogger() {
        if (file_.is_open())
            file_.close();
    }

    // Call this once per time step
    void log(int step,double timeStep=0.0) {

        // Create per-rank log file, e.g. memlog_rank000.csv
        if(file_.is_open()==false)
        {
            char fname[256];
            sprintf(fname, "%s_rank%03d.csv", basename.c_str(), rank_);
            file_.open(fname, std::ios::out);
            file_ << "Time step, Newton step,current_MB,peak_MB" << std::endl;
        }

        double current = getCurrentRSS();
        double peak = getPeakRSS();

        file_ << timeStep << ","
              << step << ","
              << std::fixed << std::setprecision(2)
              << current << "," << peak << std::endl;
        file_.flush();  // make sure it’s written to disk
    }

     void print_mem_usage() {
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);
        std::cout << "Memory: " << (r.ru_maxrss / 1024.0) << " MB" << std::endl;
    }
    double get_global_memory_usage() {
        double local_mem = get_mem_usage_mb();
        double total_mem;
        MPI_Reduce(&local_mem, &total_mem, 1, MPI_DOUBLE, MPI_SUM, 0, comm_);
        
        return total_mem;
    }

    double get_mem_usage_mb() {
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);
        return r.ru_maxrss / 1024.0; // MB
    }


private:
    MPI_Comm comm_;
    int rank_, size_;
    std::ofstream file_;

    // Peak resident set size (ru_maxrss) – kilobytes → MB
    static double getPeakRSS() {
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);
        return r.ru_maxrss / 1024.0;
    }

    // Current resident set size (Linux only, via /proc/self/statm)
    static double getCurrentRSS() {
        std::ifstream statm("/proc/self/statm");
        long size = 0, resident = 0;
        statm >> size >> resident;
        long page_kb = sysconf(_SC_PAGESIZE) / 1024;
        return (resident * page_kb) / 1024.0;
    }
};
};