#pragma once

#include <string>
#include <time.h>
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

using namespace std;

class CBenchmark {
public:
	CBenchmark(void);
	~CBenchmark(void);
	// Sets the start timepoint for the benchmark object
	void Start(void);
	// Sets the stop timepoint for the benchmark object
	void Stop(void);
	// Retrieves the elapsed cpu time
	double GetElapsedCpuTime(void);
	// Retrieves the elapsed wall time
	double GetElapsedWallTime(void);
	// Writes the elapsed CPU and wall seconds into the specified string.
	void DisplayTime(const string& programName);
private:
	// Set to true if the start timepoint has been set
	bool mIsStartSet;
	// Set to true if the stop timepoint has been set
	bool mIsStopSet;
	// Stores the CPU start timepoint
	clock_t mCpuStart;
	// Stores the CPU stop timepoint
	clock_t mCpuStop;
	// Stores the total CPU time
	clock_t mCpuTotal;
#ifdef WIN32
	// Stores the wall start timepoint
	FILETIME mWallStart;
	// Stores the wall stop timepoint
	FILETIME mWallStop;
	// Stores the total wall time
	unsigned long long mWallTotal;
#else
	// Stores the wall start timepoint
	struct timeval mWallStart;
	// Stores the wall stop timepoint
	struct timeval mWallStop;
	// Stores the total wall time
	struct timeval mWallTotal;
#endif
};
