#include "Benchmark.h"

CBenchmark::CBenchmark(void)
: mIsStartSet(false)
, mIsStopSet(false)
, mCpuTotal(0)
{
	// reset our total walltime counters
#ifdef WIN32
	mWallTotal         = 0;
#else
	mWallTotal.tv_sec  = 0;
	mWallTotal.tv_usec = 0;
#endif
}

CBenchmark::~CBenchmark(void) {}

// Sets the start timepoint for the benchmark object
void CBenchmark::Start(void) {
	mIsStartSet = true;
	mCpuStart = clock();

#ifdef WIN32
	GetSystemTimeAsFileTime(&mWallStart);
#else
	gettimeofday(&mWallStart, (struct timezone *)0);
#endif
}

// Sets the stop timepoint for the benchmark object
void CBenchmark::Stop(void) {
	mIsStopSet = true;
	mCpuStop = clock();

	// calculate the elapsed CPU time
	mCpuTotal += (mCpuStop - mCpuStart);

#ifdef WIN32
	GetSystemTimeAsFileTime(&mWallStop);

	// calculate the elapsed wall time
	unsigned long long start, stop;

	memcpy((char*)&start, (char*)&mWallStart, 8);
	memcpy((char*)&stop, (char*)&mWallStop, 8);

	mWallTotal += (stop - start);
#else
	gettimeofday(&mWallStop, (struct timezone *)0);

	// calculate the elapsed wall time
	mWallTotal.tv_sec += (mWallStop.tv_sec - mWallStart.tv_sec);

	if(mWallStop.tv_usec < mWallStart.tv_usec) {
		mWallStop.tv_usec += 1000000;
		--mWallStop.tv_sec;
	}

	mWallTotal.tv_usec += (mWallStop.tv_usec - mWallStart.tv_usec);

	if(mWallTotal.tv_usec > 1000000) {
		mWallTotal.tv_usec -= 1000000;
		++mWallStop.tv_sec;
	}
#endif
}

// Retrieves the elapsed cpu time
double CBenchmark::GetElapsedCpuTime(void) {

	// return if both timepoints haven't been set
	if(!(mIsStartSet && mIsStopSet)) return -1.0;

	// return the number of elapsed seconds
	return mCpuTotal / (double)CLOCKS_PER_SEC;
}

// Retrieves the elapsed wall time
double CBenchmark::GetElapsedWallTime(void) {

	double seconds = -1.0;

	// return if both timepoints haven't been set
	if(!(mIsStartSet && mIsStopSet)) return seconds;

#ifdef WIN32
	// calculate the difference and return the number of elapsed seconds
	seconds = mWallTotal / 10000000.0;
#else
	seconds = mWallTotal.tv_sec + mWallTotal.tv_usec / 1000000.0;
#endif

	return seconds;
}

// Writes the elapsed CPU and wall seconds into the specified string.
void CBenchmark::DisplayTime(const string& programName) {

	// retrieve the elapsed time
	double cpuSeconds  = GetElapsedCpuTime();
	double wallSeconds = GetElapsedWallTime();

	printf("%s ", programName.c_str());

	if(cpuSeconds >= 0.0) {
#ifdef WIN32
		printf("wall time: ");
		printf("%0.3f s\n", wallSeconds);
#else
		printf("CPU time: ");
		printf("%0.3f s", cpuSeconds);
		printf(", wall time: ");
		printf("%0.3f s\n", wallSeconds);
#endif
	} else printf("CPU time: error measuring cpu time.\n");
}
