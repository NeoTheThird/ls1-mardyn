/*
 * CRusage.hpp
 *
 *  Created on: Jun 23, 2011
 *      Author: Martin Schreiber, Wolfgang Eckhardt
 */

#ifndef CRUSAGE_HPP_
#define CRUSAGE_HPP_

#include "utils/Logger.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>

using namespace Log;

class CProcessMemoryInformation {

public:
	typedef unsigned long t_size;

	class CMemoryInformation {
	public:
		// program size
		t_size size;

        // resident
		t_size resident;

        // shared pages
		t_size share;

        // code
		t_size text;

        // libraries
		t_size lib;

        // data + stack
		t_size data;

	};

private:

	static CMemoryInformation memInfo;

	static std::string getStringFromSize(t_size size) {
	    char buf[32];

		if (size <= 1024) {
		    sprintf(buf, "%lu", size);

			return std::string(buf) + "B";
		}

		if (size <= 1024*1024) {
			double ns = size;
			ns /= 1024.0;

			sprintf(buf, "%f", ns);

			return std::string(buf) + "KB";
		}

		if (size <= 1024*1024*1024) {
			double ns = size;
			ns /= 1024.0*1024.0;

			sprintf(buf, "%f", ns);

			return std::string(buf) + "MB";
		}


		double ns = size;
		ns /= 1024.0*1024.0*1024.0;

		sprintf(buf, "%f", ns);

		return std::string(buf) + "GB";
	}


	/**
	 * return 0 on success.
	 */
	static int getUsageInformation(CMemoryInformation& memoryInformation) {
        FILE* file = fopen("/proc/self/statm", "r");
        if (file == NULL) {
        	global_log->error() << "ERROR: cannot open statm" << std::endl;
        	return 1;
        }

        int retval = fscanf(	file,
        		"%lu %lu %lu %lu %lu %lu",
        		&memoryInformation.size, &memoryInformation.resident, &memoryInformation.share, &memoryInformation.text, &memoryInformation.lib, &memoryInformation.data
        	);

        fclose(file);
        if (retval != 6) {
        	global_log->error() << "Error during reading values from /proc/self/statm" << std::endl;
        	return 1;
        }
        return 0;
	}

	static void outputUsageInformation(CMemoryInformation& memoryInformation) {
        // determine page size
        long page_size = sysconf(_SC_PAGESIZE);

        global_log->info() << " assuming page size " << page_size << std::endl;
        global_log->info() << " + Program size: " << getStringFromSize(memoryInformation.size*page_size) << std::endl;
        global_log->info() << " +     Resident: " << getStringFromSize(memoryInformation.resident*page_size) << std::endl;
        global_log->info() << " +       Shared: " << getStringFromSize(memoryInformation.share*page_size) << std::endl;
        global_log->info() << " +         Text: " << getStringFromSize(memoryInformation.text*page_size) << std::endl;
        global_log->info() << " +          Lib: " << getStringFromSize(memoryInformation.lib*page_size) << std::endl;
        global_log->info() << " +         Data: " << getStringFromSize(memoryInformation.data*page_size) << std::endl;
	}


public:

	static int collectMaxUsage() {
		CMemoryInformation tmpInfo;
		if (getUsageInformation(tmpInfo)) {
			return 1;
		}

		if (tmpInfo.size > memInfo.size) {
			memInfo.size = tmpInfo.size;
		}

		if (tmpInfo.resident > memInfo.resident) {
			memInfo.resident = tmpInfo.resident;
		}

		if (tmpInfo.share > memInfo.share) {
			memInfo.share = tmpInfo.share;
		}

		if (tmpInfo.text > memInfo.text) {
			memInfo.text = tmpInfo.text;
		}

		if (tmpInfo.lib > memInfo.lib) {
			memInfo.lib = tmpInfo.lib;
		}

		if (tmpInfo.data > memInfo.data) {
			memInfo.data = tmpInfo.data;
		}

		return 0;
	}



	static void outputUsageInformation() {
		CMemoryInformation memoryInformation;
		if (getUsageInformation(memoryInformation)) {
			return;
		}
		outputUsageInformation(memoryInformation);
	}


	static void outputMaxUsageInformation() {
		outputUsageInformation(memInfo);
	}
};
#endif /* CRUSAGE_HPP_ */
