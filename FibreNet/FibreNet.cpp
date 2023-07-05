#pragma once
#include "stdafx.h"
#include "SimulationManager.h"
#include "log_stream.h"

int _tmain(int argc, _TCHAR* argv[]) {
	std::shared_ptr<log_stream> p = std::shared_ptr<log_stream>(new log_stream());
	SimulationManager simulation_manager;	
	if (simulation_manager.Initialize(p)) {
		simulation_manager.Run();
		simulation_manager.Close();
	}
	system("PAUSE");
	return 0;
}

