#pragma once
#include <ostream>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>

#define MAX_DATE 128

class log_stream {

    std::ostream& out1_;
    std::fstream log_file;

public:

    log_stream() : out1_(std::cout) {
        time_t now = time(NULL);
        char the_date[MAX_DATE];
        the_date[0] = '\0';
        struct tm timeinfo;
        gmtime_s(&timeinfo, &now);

        if (now != -1) strftime(the_date, MAX_DATE, "%S_%M_%H_%d_%m_%Y", &timeinfo);

        std::string file_name = "log_" + std::string(the_date) + ".txt";
        log_file = std::fstream(file_name, std::ios_base::app);

        if (log_file.is_open()) {
            log_file << "Starting log file at:" << the_date << std::endl;
        }
        
    }
    
    ~log_stream() {
        log_file.close();
    }

    template <typename T>
    log_stream& operator<<(T const& t) {
        out1_ << t;
        log_file << t;
        return *this;
    }

    typedef log_stream& (*MyStreamManipulator)(log_stream&);

    // take in a function with the custom signature
    log_stream& operator<<(MyStreamManipulator manip) {
        // call the function, and return it's value
        return manip(*this);
    }

    // this is the type of std::cout
    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

    // this is the function signature of std::endl
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    log_stream& operator<<(StandardEndLine manip) {
        // call the function, but we cannot return it's value
        manip(std::cout);
        log_file << manip;

        return *this;
    }
};