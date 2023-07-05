#pragma once
#include <Windows.h>


class HighResolutionTimer {
public:
				  HighResolutionTimer();
	double		  GetFrequency() const;
	void		  Start();
	void		  Stop();
	double		  GetTime() const;
	double		  GetElapsedMin() const;
	double		  GetElapsedSec() const;

private:
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
	double		  elapsed;
	double		  frequency;
};

HighResolutionTimer::HighResolutionTimer() {
	frequency = GetFrequency();
}

double HighResolutionTimer::GetFrequency() const {
	LARGE_INTEGER proc_freq;
	if (!::QueryPerformanceFrequency(&proc_freq)) return -1;
	return static_cast< double >(proc_freq.QuadPart);
}

void HighResolutionTimer::Start() {
	elapsed = 0;
	DWORD_PTR oldmask = ::SetThreadAffinityMask(::GetCurrentThread(), 0);
	::QueryPerformanceCounter(&start);
	::SetThreadAffinityMask(::GetCurrentThread(), oldmask);
}

void HighResolutionTimer::Stop() {
	DWORD_PTR oldmask = ::SetThreadAffinityMask(::GetCurrentThread(), 0);
	::QueryPerformanceCounter(&stop);
	::SetThreadAffinityMask(::GetCurrentThread(), oldmask);
	elapsed = ((stop.QuadPart - start.QuadPart) / frequency);
}

double HighResolutionTimer::GetTime() const {
	LARGE_INTEGER time;
	::QueryPerformanceCounter(&time);
	return time.QuadPart / frequency;
}

double HighResolutionTimer::GetElapsedMin() const {
	return std::floor((elapsed / 60.0));
}

double HighResolutionTimer::GetElapsedSec() const {
	return elapsed;
}