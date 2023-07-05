#pragma once
#include <typeindex> 
#include <iostream>
#include <fstream>
#include <memory>
#include <array>
#include <vector>
#include "openCLManager.h"
#include "Bead.h"
#include "SpringStructure.h"
#include "AngleBond.h"
#include "BufferChunk.h"
#include "BeadInfo.h"
#include "TextParser.h"
#include "Type.h"
#include "Pair.h"
#include "SimulationSettings.h"
#include "GridParameters.h"
#include  <map>
#include  <set>
#include <algorithm>
#include <random>
#include "log_stream.h"

/*
 *
 * Main class containing host code resposible for setting 
 * up the simulation environment and running the device code.
 *
 */
class SimulationManager {
	bool									init_status;
	bool									simulation_restart;
	bool									export_position;
	bool									export_velocity;
	bool									export_force;
	bool									minimalization_reset;
	int										simulation_steps;
	int										minimalization_sub_steps;
	int										minimalization_steps;
	int										grid_update_freq;
	int										beads_num;
	int										springs_num;
	int										angles_num;
	double									skin_radius_ratio;
	int										dump_frequency;
	std::string								dump_file_name;
	std::string								restart_file_name;
	std::string								vel_dump_file_name;
	std::string								force_dump_file_name;
	SimulationSettings						settings;
	std::map<std::string,int>				type_id;
	std::vector<Type>						type_list;
	std::vector<int>						export_bead_index;
	std::vector<std::string>				export_bead_name;
	std::vector<int> 						export_type_list;
	std::shared_ptr<log_stream>				log;
	OpenclManager							openClManager;

	bool									SetDumpFile(const std::string& value);
	bool									InitializeDataDump();
	bool									DumpData(int proegress);
	bool									DumpRestartData(int proegress);
	bool									ReadInputFile(std::string inputFile);
	bool									LoadModelData(std::vector<std::string>& input_file_content);
	bool									ReadRestartFile(std::vector<std::string>& input_file_content);
	bool									ReadKernelData(std::string);
	bool									ReadBeadsData(std::vector<std::string>&	fileContents);
	bool									ReadSpringData(std::vector<std::string>& fileContents);
	bool									ReadAngleData(std::vector<std::string>&	fileContents);
	bool									RunSimulationLoop(int loop_size);
	bool									MinimizeSystemEnergy();
	bool									LoadBeadModifications(std::vector<std::string>& file_contents);
	bool									LoadSpringModifications(std::vector<std::string>& file_contents);
	bool									LoadTypeAndPairInfo(std::vector<std::string>&	file_contents);
	bool									LoadExportInfo(std::vector<std::string>&	file_contents);
	bool									LoadSimulationBoxInfo(std::vector<std::string>& file_contents);
	bool									LoadSolverInfo(std::vector<std::string>& file_contents);
	bool									LoadGlobalForceInfo(std::vector<std::string>& file_contents);

public:

	bool									Initialize(std::shared_ptr<log_stream> p);
	bool									Run();
	void									Close();
};