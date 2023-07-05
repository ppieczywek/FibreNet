#pragma once
#include "stdafx.h"
#include "SimulationManager.h"
#include "HighResolutionTimer.h"
#include <map>
#include <iomanip>

bool SimulationManager::Initialize(std::shared_ptr<log_stream> p) {
	if (p == nullptr) return false;
	log = p;
	
	*log << std::endl << std::endl;
	*log << "\t#####################################" << std::endl;
	*log << "\t#                                   #" << std::endl;
	*log << "\t#  SOLUTION MANAGER INITIALIZATION  #" << std::endl;
	*log << "\t#                                   #" << std::endl;
	*log << "\t#####################################" << std::endl << std::endl;
	
	simulation_steps		= -1;
	dump_frequency			= -1;
	beads_num				=  0;
	springs_num				=  0;
	grid_update_freq		=  8;
	minimalization_steps	=  0;
	minimalization_sub_steps = 0;
	skin_radius_ratio		= 0.5;
	simulation_restart		= false;

	export_position = false;
	export_velocity = false;
	export_force = false;

	minimalization_reset = false;

	if (!openClManager.Initialize(log)) {
		*log << "\tError occured during initialization of openClManager!" << std::endl;
		*log << "\tClosing program!" << std::endl;
		return false; 
	}

	std::string input_file;

	*log << "\tPlease specify the name of the input file: ";
	std::cin >> input_file;

	if (!ReadInputFile(input_file)) return false;
	if (!ReadKernelData("kernel.cl")) return false;

	init_status = true;
	return true;
};


bool SimulationManager::Run()
{
	if (init_status) {
		HighResolutionTimer timer;
			
		*log << std::endl << std::endl << "%-------  Running Simulation  -------%" << std::endl << std::endl;
		*log << "\tBox size:\t" << settings.size.x << " " << settings.size.y << " " << settings.size.z << std::endl;
		*log << "\tGrid cell size:\t" << settings.step.x << " x " << settings.step.y << " x " << settings.step.z << std::endl;
		*log << "\tGrid cells:\t" << settings.x_cells_num << " x " << settings.y_cells_num << " x " << settings.z_cells_num << std::endl;
		*log << "\tSimulation steps:\t" << simulation_steps << std::endl;
		*log << "\tData dump frequency:\t" << dump_frequency << std::endl;
		*log << "\tTime step:\t" << settings.dt << std::endl << std::endl << std::endl;

		if (!InitializeDataDump()) return false;

		timer.Start();

		if(minimalization_steps > 0) {
			if (!MinimizeSystemEnergy()) return false;}

		int progress = 0;
		if (!DumpData(progress)) return false;

		while(progress < simulation_steps) {
			if ((simulation_steps - progress) < dump_frequency) {
				if (!RunSimulationLoop(simulation_steps - progress)) return false;
				progress += simulation_steps - progress;
			} else {
				if (!RunSimulationLoop(dump_frequency)) return false;
				progress += dump_frequency;
			}
			
			if (!DumpData(progress)) return false;
			*log << "\tSimulation progress: finished " << progress << " out of " << simulation_steps << " simulation steps" << std::endl;
		}

		timer.Stop();

		if (!DumpRestartData(progress)) {
			*log << "\tError: unable to dump final restart data frame." << std::endl;
			return false;
		}

		*log << "\tSimulation finished after " << timer.GetElapsedMin() << " minutes" << std::endl;
		*log << "\tClosing program." << std::endl;

		return true;
	}
	return false;
};

/*
 * Runs simulation energy minimization loop by sequentaily executing pre-defined kernel functions.
 *
 * @param loop_size - number of simulation steps to execute
 * @return bool, true if manage to all simulation steps, otherwise false
 */
bool SimulationManager::MinimizeSystemEnergy()
{
	*log << "\tMinimizeSystemEnergy: initializing system energy minimalization." << std::endl; 
	if (minimalization_steps < 0) {
		*log << "\tMinimizeSystemEnergy: ERROR - loop size should be a positive number." << std::endl;
		return false;
	}
		
	if (minimalization_reset == true) {
		if (!openClManager.ExecuteKernel("ResetBeadEnergy")) return false;
	}
	
	if (minimalization_steps == 0) {
		*log << "\tMinimizeSystemEnergy: number of minimalization steps set to zero, skipping minimalization." << std::endl;
		return true;
	}

	int loops_number = minimalization_steps / minimalization_sub_steps;
	for (int pp = 0; pp < loops_number; pp++) {
		if (!RunSimulationLoop(minimalization_sub_steps)) {
			return false;
		} else {
			*log << "\tMinimizeSystemEnergy: minimalization progress - finished " << pp * minimalization_sub_steps
					  << " out of " << minimalization_sub_steps * loops_number << " minimalization steps" << std::endl;
		}

		if (minimalization_reset == true) {
			if (!openClManager.ExecuteKernel("ResetBeadEnergy")) return false;
		}
	}
	*log << "\tMinimizeSystemEnergy: system energy minimalization finished." << std::endl << std::endl;
	return true;
};

/*
 * Runs simulation main loop by sequentaily executing pre-defined kernel functions.
 *
 * @param loop_size - number of simulation steps to execute
 * @return bool, true if manage to all simulation steps, otherwise false
 */
bool SimulationManager::RunSimulationLoop(int loop_size) {
	int		   grid_update_cnt			= 0;
	if (loop_size <= 0) {
		*log << "\tRunSimulationLoop: ERROR - loop size should be a positive number" << std::endl;
		return false; }

	for (auto step = 0; step < loop_size; step++) {
		if (grid_update_cnt == grid_update_freq) grid_update_cnt = 0;
		if (!openClManager.ExecuteKernel("AdvancePosition")) return false;
		if (springs_num > 0) {
			if (!openClManager.ExecuteKernel("ResolveSprings")) return false;
		}

		if (angles_num > 0) {
			if (!openClManager.ExecuteKernel("ResolveAngles")) return false;
		}

		if (grid_update_cnt == 0) {
			if (!openClManager.ExecuteKernel("ResetGrid")) return false;
			if (!openClManager.ExecuteKernel("ResetContactList")) return false;
			if (!openClManager.ExecuteKernel("FillGrid")) return false;
			if (!openClManager.ExecuteKernel("BuildContactList")) return false;
		}

		if (!openClManager.ExecuteKernel("ResolveContacts")) return false;
		if (!openClManager.ExecuteKernel("AdvanceVelocity")) return false;
		if (!openClManager.RunKernels()) return false;
		grid_update_cnt++; 
	}
	return true;
};

/*
 * Reads new settings of selected bead structures from input file.
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to set new data of beads, otherwise false
 */
bool SimulationManager::LoadBeadModifications(std::vector<std::string>& file_contents)
{
	if (!openClManager.GetStatus()) {
		*log << "\tLoadBeadModifications: ERROR - incorrect number of uploaded beads." << std::endl;
		return false;
	}

	if (BufferIndex::bead_buffer < 0) {
		*log << "\tLoadBeadModifications: ERROR - bead buffer not initialized." << std::endl;
		return false;
	}

	if (BufferIndex::bead_info_buffer < 0) {
		*log << "\tLoadBeadModifications: ERROR - bead info buffer not initialized." << std::endl;
		return false;
	}

	if (openClManager.GetBufferLength(BufferIndex::bead_buffer) == -1) {
		*log << "\tLoadBeadModifications: ERROR - bead buffer is empty." << std::endl;
		return false;
	}

	if (openClManager.GetBufferLength(BufferIndex::bead_info_buffer) == -1)
	{
		*log << "\tLoadBeadModifications: ERROR - bead info buffer is empty." << std::endl;
		return false;
	}

	std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
	std::vector<BeadInfo> bead_info(openClManager.GetBufferLength(BufferIndex::bead_info_buffer), BeadInfo());

	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) {
		*log << "\tLoadBeadModifications: ERROR - unagle to read bead buffer data." << std::endl;
		return false;
	}

	if (!openClManager.ReadBuffer(BufferIndex::bead_info_buffer, &bead_info.front())) {
		*log << "\tLoadBeadModifications: ERROR - unagle to read bead info buffer data." << std::endl;
		return false;
	}

	if (type_id.empty()) {
		*log << "\tLoadBeadModifications: ERROR - type list not initialized." << std::endl;
		return false;
	}

	std::vector<std::string> tokens;
	bool modified = false;
	int my_start = 0;

	while ((my_start = TextParser::ScanText(file_contents, "modif_bead", tokens, my_start)) != -1) {
		if (tokens.size() > 4) {
			int type = -100;
			selection_type selection = none;

			if (tokens[1] == "id") {
				selection = by_id;
			} else if (tokens[1] == "tag") {
				selection = by_tag;
			} else if (tokens[1] == "type") {
				selection = by_type;
			} else {
				return false;
			}

			if (selection == by_type) {
				std::string type_name = tokens[2];
				if (type_id.find(type_name) == type_id.end()) {
					*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - first interction type id does not match any of the defined types." << std::endl;
					return false;
				}
				type = type_id.find(type_name)->second;
			}

			if (selection == by_id || selection == by_tag) {
				type = std::stoi(tokens[2]);
			}

			if (tokens[3] == "mass") {
				*log << "\tLoadBeadModifications: applying bead mass modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_mass = std::stod(tokens[4]);
					if (new_mass >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::bead_info_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_type) {
									if (bead_info[ii].type == type)	match_found = true;
								}
								
								if (selection == by_tag) {
									if (bead_info[ii].tag == type) match_found = true;
								}

								if (selection == by_id)	{
									if (ii == type) match_found = true;
								}

								if(match_found)	{
									bead_info[ii].mass = new_mass;
									if (bead_info[ii].mass > 0.0f) {
										bead_info[ii].mass_inv = 1.0f / bead_info[ii].mass;
										bead_info[ii].lambda_dt_mass_inv = settings.lambda_dt * bead_info[ii].mass_inv;
										bead_info[ii].half_dt_mass_inv = settings.half_dt * bead_info[ii].mass_inv;
										bead_info[ii].half_dt2_mass_inv = settings.half_dt2 * bead_info[ii].mass_inv;
									} else {
										bead_info[ii].mass_inv = 0.0f;
										bead_info[ii].lambda_dt_mass_inv = 0.0f;
										bead_info[ii].half_dt_mass_inv = 0.0f;
										bead_info[ii].half_dt2_mass_inv = 0.0f;
									}
								}
							}
							modified = true;
						} else {
							*log << "\tLoadBeadModifications: ERROR - bead_info_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - incorrect value for 'mass', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_bead	  mass'." << std::endl;
					return false;
				}
			} else if (tokens[3] == "velocity") {
				*log << "\tLoadBeadModifications: applying bead velocity modifications" << std::endl;
				if (tokens.size() == 7) {
					Vector3 new_velocity;
					new_velocity.x = std::stod(tokens[4]);
					new_velocity.y = std::stod(tokens[5]);
					new_velocity.z = std::stod(tokens[6]);

					int length = openClManager.GetBufferElements(BufferIndex::bead_buffer);
					if (length > 0) {
						bool match_found = false;
						for (auto ii = 0; ii < length; ii++) {
							match_found = false;
							if (selection == by_type) {
								if (bead_info[ii].type == type)	match_found = true;
							}

							if (selection == by_tag) {
								if (bead_info[ii].tag == type) match_found = true;
							}

							if (selection == by_id) {
								if (ii == type)	match_found = true;
							}

							if (match_found) {
								bead_data[ii].velocity = new_velocity;
								bead_data[ii].velocity_old = new_velocity;
							}
						}
						modified = true;
					} else {
						*log << "\tLoadBeadModifications: ERROR - bead_buffer is empty." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_bead	  velocity'." << std::endl;
					return false;
				}
			} else if (tokens[3] == "skin_radius") {
				*log << "\tLoadBeadModifications: applying bead skin_radius modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_skin_radius;
					new_skin_radius = std::stod(tokens[4]);
					int length = openClManager.GetBufferElements(BufferIndex::bead_info_buffer);
					if (length > 0) {
						bool match_found = false;
						for (auto ii = 0; ii < length; ii++) {
							match_found = false;
							if (selection == by_type) {
								if (bead_info[ii].type == type) match_found = true;
							}

							if (selection == by_tag) {
								if (bead_info[ii].tag == type) match_found = true;
							}

							if (selection == by_id) {
								if (ii == type) match_found = true;
							}

							if (match_found) {
								bead_info[ii].skin_radius = new_skin_radius;
								bead_info[ii].skin_radius_sq = new_skin_radius * new_skin_radius;
							}
						}
						modified = true;
					} else {
						*log << "\tLoadBeadModifications: ERROR - bead_buffer is empty." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_bead	  skin_radius'." << std::endl;
					return false;
				}
			} else {
				*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - unrecognized second parameter of 'modif_bead' command." << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadBeadModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_bead' command." << std::endl;
			return false;
		}
	}

	if (modified == true) {
		if (!openClManager.WriteBuffer<Bead>(BufferIndex::bead_buffer, bead_data)) {
			*log << "\tLoadBeadModifications: ERROR - unable to write 'BufferIndex::bead_buffer'" << std::endl;
			return false;
		}

		if (!openClManager.WriteBuffer<BeadInfo>(BufferIndex::bead_info_buffer, bead_info)) {
			*log << "\tLoadBeadModifications: ERROR - unable to write 'BufferIndex::bead_info_buffer'" << std::endl;
			return false;
		}
		*log << "\tLoadBeadModifications: bead_buffer updated with new data" << std::endl;
	}
	return true;
}

/*
 * Reads new settings of selected spring structures from input file.
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to set new data of spring bonds, otherwise false
 */
bool SimulationManager::LoadSpringModifications(std::vector<std::string>& file_contents) {
	if (!openClManager.GetStatus()) {
		*log << "\tLoadSpringModifications: ERROR - OpenClManager not initialized." << std::endl;
		return false;
	}

	if (BufferIndex::spring_buffer < 0) {
		*log << "\tLoadSpringModifications: ERROR - spring buffer not initialized." << std::endl;
		return false;
	}

	if (openClManager.GetBufferLength(BufferIndex::spring_buffer) == -1) {
		*log << "\tLoadSpringModifications: ERROR - spring buffer is empty." << std::endl;
		return false;
	}

	std::vector<SpringStructure> spring_data(openClManager.GetBufferLength(BufferIndex::spring_buffer), SpringStructure());
	if (!openClManager.ReadBuffer(BufferIndex::spring_buffer, &spring_data.front())) {
		*log << "\tLoadSpringModifications: ERROR - unagle to read spring buffer data." << std::endl;
		return false;
	}

	std::vector<std::string> tokens;
	bool modified = false;
	int my_start = 0;

	while ((my_start = TextParser::ScanText(file_contents, "modif_spring", tokens, my_start)) != -1) {
		if (tokens.size() > 4) {
			int type = -100;
			selection_type selection = none;

			if (tokens[1] == "id") {
				selection = by_id;
			} else if (tokens[1] == "tag") {
				selection = by_tag;
			} else {
				//unknown selection method
				return false;
			}

			if (selection == by_id || selection == by_tag) type = std::stoi(tokens[2]);

			if (tokens[3] == "c1") {
				*log << "\tLoadSpringModifications: applying spring c1 modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_c1 = std::stod(tokens[4]);
					if (new_c1 >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::spring_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_tag) {
									if (spring_data[ii].tag == type) match_found = true;
								}
								if (selection == by_id)	{
									if (ii == type)	match_found = true;
								}
								if (match_found) spring_data[ii].c1 = new_c1;
							}
							modified = true;
						} else {
							*log << "\tLoadSpringModifications: ERROR - spring_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect value for 'c1', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring	  c1'." << std::endl;
					return false;
				}
			} else if (tokens[3] == "c2") {
				*log << "\tLoadSpringModifications: applying spring c2 modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_c2 = std::stod(tokens[4]);
					if (new_c2 >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::spring_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_tag) {
									if (spring_data[ii].tag == type) match_found = true;
								}
								if (selection == by_id) {
									if (ii == type)	match_found = true;
								}
								if (match_found) spring_data[ii].c2 = new_c2;
							}
							modified = true;
						} else {
							*log << "\tLoadSpringModifications: ERROR - spring_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect value for 'c2', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring	  c2'." << std::endl;
					return false;
				}
			} else if (tokens[3] == "rest_length") {
				*log << "\tLoadSpringModifications: applying spring rest_lenght modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_rest_length = std::stod(tokens[4]);
					if (new_rest_length >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::spring_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_tag) {
									if (spring_data[ii].tag == type) match_found = true;
								}
								if (selection == by_id)	{
									if (ii == type) match_found = true;
								}
								if (match_found) spring_data[ii].rest_length = new_rest_length;
							}
							modified = true;
						} else {
							*log << "\tLoadSpringModifications: ERROR - spring_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect value for 'rest_length', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring	 rest_length'." << std::endl;
					return false;
				}
			} else if (tokens[3] == "stiffness") {
				*log << "\tLoadSpringModifications: applying spring stiffness modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_stiffness = std::stod(tokens[4]);
					if (new_stiffness >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::spring_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_tag) {
									if (spring_data[ii].tag == type) match_found = true;
								}
								if (selection == by_id) {
									if (ii == type)	match_found = true;
								}
								if (match_found) spring_data[ii].stiffness = new_stiffness;
							}
							modified = true;
						} else {
							*log << "\tLoadSpringModifications: ERROR - spring_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect value for 'stiffness', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring	 stiffness'." << std::endl;
				}
			} else if (tokens[3] == "damping") {
				*log << "\tLoadSpringModifications: applying spring damping modifications" << std::endl;
				if (tokens.size() == 5) {
					double new_damping = std::stod(tokens[4]);
					if (new_damping >= 0.0f) {
						int length = openClManager.GetBufferElements(BufferIndex::spring_buffer);
						if (length > 0) {
							bool match_found = false;
							for (auto ii = 0; ii < length; ii++) {
								match_found = false;
								if (selection == by_tag) {
									if (spring_data[ii].tag == type) match_found = true;
								}
								if (selection == by_id) {
									if (ii == type)	match_found = true;
								}
								if (match_found) spring_data[ii].damping = new_damping;
							}
							modified = true;
						} else {
							*log << "\tLoadSpringModifications: ERROR - spring_buffer is empty." << std::endl;
							return false;
						}
					} else {
						*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect value for 'damping', should be zero or positive real number." << std::endl;
						return false;
					}
				} else {
					*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring	 damping'." << std::endl;
				}
			} else {
				*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - unrecognized second parameter of 'modif_spring' command." << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSpringModifications: ERROR - line " << my_start << " - incorrect number of arguments for 'modif_spring' command." << std::endl;
			return false;
		}
	}

	if (modified == true) {
		if (!openClManager.WriteBuffer<SpringStructure>(BufferIndex::spring_buffer, spring_data)) {
			*log << "\tLoadSpringModifications: ERROR - unable to write 'BufferIndex::spring_buffer'" << std::endl;
			return false;
		}
		*log << "\tLoadSpringModifications: spring_buffer updated with new data" << std::endl;
	}
	return true;
}

/*
 * Reads general simulation solver settings from input file.
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to read solver data, otherwise false
 */
bool SimulationManager::LoadGlobalForceInfo(std::vector<std::string>& file_contents) {
	std::vector<std::string>	tokens;
	int	input_file_line = 0;
	if ((input_file_line = TextParser::ScanText(file_contents, "gravity", tokens, 0)) > 0) {
		if (tokens.size() == 4) {
			settings.gravity.x = std::stod(tokens[1]);
			settings.gravity.y = std::stod(tokens[2]);
			settings.gravity.z = std::stod(tokens[3]);
		} else {
			*log << "\tLoadGlobalForceInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'gravity'" << std::endl;
			return false;
		}
	} else {
		settings.gravity.x = 0.0f;
		settings.gravity.y = 0.0f;
		settings.gravity.z = 0.0f;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "damping", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			settings.damping = std::stod(tokens[1]) * -1.0f;
		} else {
			*log << "\tLoadGlobalForceInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'damping'" << std::endl;
			return false;
		}
	} else {
		settings.damping = 0.0f;
	}
	return true;
}

/*
 * Reads general simulation solver settings from input file.
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to read solver data, otherwise false
 */
bool SimulationManager::LoadSolverInfo(std::vector<std::string>& file_contents) {
	int	input_file_line = 0;
	std::vector<std::string>	tokens;
	if ((input_file_line = TextParser::ScanText(file_contents, "step", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			settings.dt = std::stod(tokens[1]); 
			if (settings.dt > 0.0f)	{
				settings.half_dt = settings.dt * 0.5f;
				settings.half_dt2 = settings.dt * settings.dt * 0.5f;
			} else {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'step'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'step'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadSolverInfo: ERROR - required command value not specified in input file: 'step'" << std::endl;
		return false;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "lambda", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			settings.lambda_dt = std::stod(tokens[1]);
			if (settings.lambda_dt > 0.0f) {
				settings.lambda_dt *= settings.dt;
			} else {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'lambda'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'lambda'" << std::endl;
			return false;
		}
	} else {
		settings.lambda_dt = 0.5f * settings.dt;
		*log << "\tLoadSolverInfo: WARNING - command value not specified in input file, using default value: 'lambda = 0.5'" << std::endl;
		return false;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "run", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			simulation_steps = std::stoi(tokens[1]);
			if (simulation_steps < 0) {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'run'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'run'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadSolverInfo: ERROR - required command value not specified in input file: 'run'" << std::endl;
		return false;
	}


	if ((input_file_line = TextParser::ScanText(file_contents, "minimize_steps", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			minimalization_steps = std::stoi(tokens[1]);
			if (minimalization_steps < 0) {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'minimize'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'minimize'" << std::endl;
			return false;
		}
	} else {
		minimalization_steps = 0;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "minimize_substeps", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			minimalization_sub_steps = std::stoi(tokens[1]);
			if (minimalization_sub_steps < 0) {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'minimize_step'" << std::endl;
				return false;
			}

			if (minimalization_sub_steps > minimalization_steps) {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'minimize_step'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'minimize_step'" << std::endl;
			return false;
		}
	} else {
		if (minimalization_steps > 0) {
			*log << "\tLoadSolverInfo: WARNING - command value not specified in input file, using default value: 'minimize_step = minimize'" << std::endl;
			minimalization_sub_steps = minimalization_steps;
		} else {
			minimalization_sub_steps = 0;
		}
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "minimize_reset", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			if (tokens[1] == "yes") {
				minimalization_reset = true;
			} else if (tokens[1] == "no") {
				minimalization_reset = false;
			} else {
				*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect command value: 'minimize_reset'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'minimize_reset'" << std::endl;
			return false;
		}
	} else {
		if (minimalization_steps > 0) {
			*log << "\tLoadSolverInfo: WARNING - command value not specified in input file, using default value: 'minimize_reset = false'" << std::endl;
			minimalization_reset = false;
		} else {
			minimalization_reset = false;
		}
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "cutoff", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			settings.cutoff = std::stod(tokens[1]);
			settings.cutoff_sq = settings.cutoff * settings.cutoff;
			if (settings.cutoff <= 0.0f) {
				*log << "\tLoadSolverInfo: WARNING - line " << input_file_line << " - incorrect value of command parameter: 'cutoff'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSolverInfo: WARNING - line " << input_file_line << " - incorrect number of command parameters: 'cutoff'" << std::endl;
			return false;
		}
	} else {
		std::cout << "\tLoadSolverInfo: WARNING -  line " << input_file_line << " - cutoff radius was not specified in the input file." << std::endl;
		return false;
	}
	return true;
};

/*
 * Reads the prameters of non-bonded bead-to-bead interactions
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to read intearaction pairs data, otherwise false
 */
bool SimulationManager::LoadTypeAndPairInfo(std::vector<std::string>& file_contents)
{
	std::vector<std::string>	tokens;
	int input_file_line = 0;
	int types_num = 0;
	while ((input_file_line = TextParser::ScanText(file_contents, "type", tokens, input_file_line)) != -1) {
		if (tokens.size() == 2) {
			auto id = tokens[1];
			auto it = type_id.find(id);

			if (it == type_id.end()) {
				Type type;
				type.id = id;
				type_list.push_back(type);
				type_id[id] = types_num;
				types_num++;
			} else {
				*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line << " - bead/paticle type redefinition." << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line << " - incorrect bead/paticle type definition." << std::endl;
			return false;
		}
	}

	if (types_num <= 0) {
		*log << "\tLoadTypeAndPairInfo: ERROR - unable to find bead/paticle type definitions." << std::endl;
		return false;
	}

	settings.types_num = types_num;
	int pairs_num = ((types_num*(types_num - 1)) / 2) + types_num;
	std::vector<Pair> pair_list(pairs_num, Pair());
	for (int ii = 0; ii < pair_list.size(); ii++) {
		pair_list[ii].p1 = -1;
		pair_list[ii].p2 = -1;
	}
	
	input_file_line = 0;
	while ((input_file_line = TextParser::ScanText(file_contents, "pair", tokens, input_file_line)) != -1) {
		if (tokens.size() > 4) {
			std::string type_name_1 = tokens[2];
			std::string type_name_2 = tokens[3];

			if (type_id.find(type_name_1) == type_id.end()) {
				*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
					" - first intearction type id does not match any of the defined types." << std::endl;
				return false;
			}
			if (type_id.find(type_name_2) == type_id.end()) {
				*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
					" - second intearction type id does not match any of the defined types." << std::endl;
				return false;
			}

			int type_1 = type_id.find(type_name_1)->second;
			int type_2 = type_id.find(type_name_2)->second;
			if (type_1 > type_2) std::swap(type_1, type_2);
			int pair_list_index = (types_num * type_1) + type_2 - ((type_1*(type_1 + 1)) / 2);

			if (pair_list_index >= 0 && pair_list_index < pair_list.size()) {
				Pair interaction_pair;
				
				if (tokens[1] == "lj") {
					if (tokens.size() == 7) {
						interaction_pair.type = 1;
						interaction_pair.p1 = type_1;
						interaction_pair.p2 = type_2;
						interaction_pair.c1 = std::stod(tokens[4]); 
						interaction_pair.c2 = std::stod(tokens[5]); 
						interaction_pair.c3 = std::stod(tokens[6]); 
						double r2 = interaction_pair.c3 * interaction_pair.c3;
						double sigma2 = interaction_pair.c1;
						double epsilon = interaction_pair.c2;
						sigma2 *= sigma2;

						double fr2 = sigma2 / r2;
						double fr6 = fr2 * fr2 * fr2;
						double fpr = 48.0 * epsilon * fr6 * (fr6 - 0.5) / r2;
						interaction_pair.c4 = fpr;
					} else {
						*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
							" - lj intearction definition contains incorrect number of parameters." << std::endl;
					}
				} else if (tokens[1] == "hs") {
					if (tokens.size() == 5) {
						interaction_pair.type = 2;
						interaction_pair.p1 = type_1;
						interaction_pair.p2 = type_2;
						interaction_pair.c1 = std::stod(tokens[4]);
					} else {
						*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
							" - hs intearction definition contains incorrect number of parameters." << std::endl;
					}
				} else if (tokens[1] == "st") {
					if (tokens.size() == 8) {
						interaction_pair.type = 3;
						interaction_pair.p1 = type_1;
						interaction_pair.p2 = type_2;
						interaction_pair.c1 = std::stod(tokens[5]); 
						interaction_pair.c2 = std::stod(tokens[6]); 
						interaction_pair.c3 = std::stod(tokens[7]) - interaction_pair.c2; 
						interaction_pair.c4 = std::stod(tokens[7]) * std::stod(tokens[7]); 
						interaction_pair.c5 = std::stod(tokens[4]); 
						interaction_pair.c6 = 2 * std::stod(tokens[6]) - std::stod(tokens[7]);
					} else {
						*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
							" - st intearction definition contains incorrect number of parameters." << std::endl;
					}
				} else {
					*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
						" - unknown intearction typ definition." << std::endl;
					return false;
				}
				pair_list[pair_list_index] = interaction_pair;
			} else {
				*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
					" - intearction type index exceedes interaction type table bounds." << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadTypeAndPairInfo: ERROR - line " << input_file_line <<
				" - incomplete intearction type definition." << std::endl;
			return false;
		}
	}

	for (int ii = 0; ii < pair_list.size(); ii++) {
		if (pair_list[ii].p1 == -1 && pair_list[ii].p2 == -1) {
			*log << "\tLoadTypeAndPairInfo: ERROR - pair list contains undefined interactions. Check if all type - type intearctions are defined in input file." << std::endl;
			return false;
		}
	}

	if (pairs_num > 0) {
		BufferIndex::bead_pair_table_buffer = openClManager.CreateBuffer<Pair>(pair_list, CL_MEM_READ_ONLY);
		if (BufferIndex::bead_pair_table_buffer < 0) {
			*log << "\tLoadTypeAndPairInfo: ERROR - unable to create 'bead_pair_table_buffer'." << std::endl;
			return false;
		}
	}
	return true;
};

/*
 * Reads the prameters of the data dump procedure based on
 * the information from the input file
 *
 * @param file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to read data damp info, otherwise false
 */
bool SimulationManager::LoadExportInfo(std::vector<std::string>& file_contents)
{
	std::vector<std::string>	tokens;
	int	input_file_line = 0;

	if (type_id.empty()) {
		*log << "\tLoadExportInfo: ERROR - Before reading export data load the particle type data from the input file." << std::endl;
		return false;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "output_file", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			dump_file_name = tokens[1] + "_coord.xyz";
			vel_dump_file_name = tokens[1] + "_vel.xyz";
			force_dump_file_name = tokens[1] + "_force.xyz";
			restart_file_name = tokens[1] + "_restart.res";
		} else {
			*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'output'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadExportInfo: ERROR - output file name not specified. Use 'output' command to set the output file name." << std::endl;
		return false;
	}

	if ((input_file_line = TextParser::ScanText(file_contents, "output_freq", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			dump_frequency = std::stoi(tokens[1]);
			if (dump_frequency <= 0) {
				*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect value for data dump freqeancy. Positive integer expected." << std::endl; //linia
				return false;
			}

			if (dump_frequency > simulation_steps) {
				*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect value for data dump freqeancy. Positive integer expected smaller than number of total simulation steps." << std::endl; //linia
				return false;
			}
		} else {
			*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'output_freq'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadExportInfo: ERROR - data dump frequency was not specified. Use 'output_freq' command to set the data dump frequency." << std::endl;
		return false;
	}

	input_file_line = 0;
	while ((input_file_line = TextParser::ScanText(file_contents, "output_type", tokens, input_file_line)) != -1) {
		if (tokens.size() == 2) {
			std::string type = tokens[1];
			if (type == "all") {
				if (export_type_list.size() > 0) {
					export_type_list.clear();
					*log << "\tLoadExportInfo: WARNING - line " << input_file_line << "\t" << " - exporting all possible types. Current type export list will be overwritten." << std::endl;
				}
					
				for (auto it = type_id.begin(); it != type_id.end(); it++) 	{
					export_type_list.push_back(it->second);
				}
			} else {
				if (type_id.find(type) == type_id.end()) {
					*log << "\tLoadExportInfo: ERROR - line " << input_file_line << "\t" << " - export bead type deos not match any of the defined types." << std::endl;
					return false;
				}

				int export_type = type_id.find(type)->second;
				if (find(export_type_list.begin(), export_type_list.end(), export_type) == export_type_list.end()) {
					export_type_list.push_back(export_type);
				} else {
					*log << "\tLoadExportInfo: WARNING - line " << input_file_line << "\t" << " - export bead type already on list." << std::endl;
				}
			}
		} else {
			*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'export'" << std::endl;
			return false;
		}
	}

	if (export_type_list.empty()) {
		*log << "\tLoadExportInfo: ERROR - list of exported bead types is empty. Please define types of objects to export." << std::endl;
		return false;
	}

	input_file_line = 0;
	while ((input_file_line = TextParser::ScanText(file_contents, "output_data", tokens, input_file_line)) != -1) {
		if (tokens.size() == 2) {

			if (tokens[1] == "position") {
				export_position = true;
			} else if (tokens[1] == "velocity") {
				export_velocity = true;
			} else if (tokens[1] == "force") {
				export_force = true;
			} else if (tokens[1] == "all") {
				export_velocity = true;
				export_position = true;
				export_force = true;
			} else {
				*log << "\tLoadExportInfo: ERROR input file line - " << input_file_line << "\t" << " - export data type label deos not match any of the defined types.!" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadExportInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'export'" << std::endl;
			return false;
		}
	}

	if (export_position == false && export_velocity == false && export_force == false) {
		*log << "\tLoadExportInfo: ERROR - export data type not specified. Pleas specify at least one type to export (position velocity or force)." << std::endl;
		return false;
	}
	return true;
};

/*
 * Pre-defines the prameters of the data dump procedure based on 
 * the information gathered from the input file
 *
 * @return bool, true if manage to initialize data damp procedure, otherwise false
 */
bool SimulationManager::InitializeDataDump()
{
	if (!openClManager.GetStatus()) {
		*log << "\tInitializeDataDump: ERROR - openClManager not initialized." << std::endl;
		return false;
	}

	if (BufferIndex::bead_info_buffer < 0) {
		*log << "\tInitializeDataDump: ERROR - ufferIndex::bead_info_buffer not initialized." << std::endl;
		return false;
	}
	std::vector<BeadInfo> bead_info_buffer(openClManager.GetBufferLength(BufferIndex::bead_info_buffer), BeadInfo());

	if (bead_info_buffer.size() == 0) {
		*log << "\tInitializeDataDump: ERROR - unable to initialize bead_info_buffer ." << std::endl;
		return false;
	}

	if (!openClManager.ReadBuffer(BufferIndex::bead_info_buffer, &bead_info_buffer.front())) {
		*log << "\tInitializeDataDump: ERROR - unable to load bead_info_buffer ." << std::endl;
		return false;
	}
	
	if (export_type_list.empty()) {
		*log << "\tInitializeDataDump: ERROR - export_type_list is empty." << std::endl;
		return false;
	}

	auto number_of_beads = openClManager.GetBufferElements(BufferIndex::bead_info_buffer);
	for (auto jj = 0; jj < number_of_beads; jj++) {
		if (std::find(export_type_list.begin(), export_type_list.end(), bead_info_buffer[jj].type) != export_type_list.end()) {
			auto type = bead_info_buffer[jj].type;
			if (type < 0 || type > type_list.size()) {
				*log << "\tInitializeDataDump: ERROR - bead type does not match any of types from declared type list." << std::endl;
				return false;
			}
			export_bead_index.push_back(jj);
			export_bead_name.push_back(type_list[type].id);
		}
	}
	
	if (export_bead_index.empty() && export_bead_name.empty()) {
		*log << "\tInitializeDataDump: ERROR - declared export types not found in model - there is no data to export." << std::endl;
		return false;
	}
	return true;
};

/*
 * Reads data from input file and initializes simulation
 *
 * @param input_file_path - path to input file
 * @return bool, true if manage to read dat from input file, otherwise false
 */
bool SimulationManager::ReadInputFile(std::string input_file_path) {	
	*log << std::endl << std::endl << "\tReadInputFile: reading simulation input file -" << input_file_path << std::endl << std::endl;
	
	std::vector<std::string>	file_contents;
	if (!TextParser::ReadFileContent(input_file_path, file_contents)) {
		*log << "\tReadInputFile: ERROR - unable to open input file:" << input_file_path << std::endl;
		return false;
	}
	//reads blocks of data from simulation input file
	if (!LoadGlobalForceInfo(file_contents)) return false;
	if (!LoadSolverInfo(file_contents)) return false;
	if (!LoadSimulationBoxInfo(file_contents)) return false;
	if (!LoadTypeAndPairInfo(file_contents)) return false;
	if (!LoadExportInfo(file_contents)) return false;
	if (!LoadModelData(file_contents)) return false;
	if (!LoadBeadModifications(file_contents)) return false;
	if (!LoadSpringModifications(file_contents)) return false;
	return true;
};

/*
 * Reads kernel code from external kernel text file, 
 * sets input variables and registers on OpenCl platform device.
 *
 * @param input_file - path to OpenCl kernel file
 * @return bool, true if manage to perform register kernels, otherwise false
 */
bool SimulationManager::ReadKernelData(std::string input_file) {
	//read OpenCl kernels code
	std::ifstream file(input_file);
	if (!file.is_open())	return false;
	if (!file.good())		return false;
	std::string	 kernel_code((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

	// name all kernel functions
	// prvided names must match names of kernel 
	// functions from kernel code text file
	std::vector<std::string> kernel_names;
	kernel_names.push_back("AdvancePosition");
	kernel_names.push_back("ResetGrid");
	kernel_names.push_back("ResetContactList");
	kernel_names.push_back("FillGrid");
	kernel_names.push_back("BuildContactList");
	kernel_names.push_back("ResolveContacts");
	kernel_names.push_back("AdvanceVelocity");
	kernel_names.push_back("ResetBeadEnergy");
	kernel_names.push_back("ResolveSprings");
	kernel_names.push_back("ResolveAngles");

	//register all kernel functions
	if (!openClManager.RegisterKernel(kernel_names, kernel_code)) {
		return false;
	}

	//get instances of all data buffers
	auto bead_buffer = openClManager.GetBuffer(BufferIndex::bead_buffer);
	auto angle_buffer = openClManager.GetBuffer(BufferIndex::angle_buffer);
	auto spring_buffer = openClManager.GetBuffer(BufferIndex::spring_buffer);
	auto bead_info_buffer = openClManager.GetBuffer(BufferIndex::bead_info_buffer);
	auto grid_head_buffer = openClManager.GetBuffer(BufferIndex::grid_head_buffer);
	auto grid_tail_buffer = openClManager.GetBuffer(BufferIndex::grid_tail_buffer);
	auto grid_cell_nh_buffer = openClManager.GetBuffer(BufferIndex::grid_cell_nh_buffer);
	auto bead_nh_head_buffer = openClManager.GetBuffer(BufferIndex::bead_nh_head_buffer);
	auto bead_nh_tail_buffer = openClManager.GetBuffer(BufferIndex::bead_nh_tail_buffer);
	auto bead_pair_coeff_buffer = openClManager.GetBuffer(BufferIndex::bead_pair_coeff_buffer);
	auto bead_pair_table_buffer = openClManager.GetBuffer(BufferIndex::bead_pair_table_buffer);

	if (bead_buffer == nullptr) return false; 
	if (grid_head_buffer == nullptr) return false;
	if (grid_tail_buffer == nullptr) return false;
	if (grid_cell_nh_buffer == nullptr) return false;
	if (bead_nh_head_buffer == nullptr) return false;
	if (bead_nh_tail_buffer == nullptr) return false;
	if (bead_pair_coeff_buffer == nullptr) return false;
	if (bead_pair_table_buffer == nullptr) return false;

	int beads_num = openClManager.GetBufferElements(BufferIndex::bead_buffer);
	int grid_size = openClManager.GetBufferElements(BufferIndex::grid_head_buffer);
	int list_size = openClManager.GetBufferElements(BufferIndex::bead_nh_head_buffer);

	// Assigns input paramters of all kernel functions.
	// Order and type of input must be the same as defined in 
	// kernel code text file.
	if (!openClManager.SetKernelArg("AdvancePosition", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 1, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 2, settings)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 3, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("AdvancePosition", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	if (spring_buffer != nullptr) {
		int springs_num = openClManager.GetBufferElements(BufferIndex::spring_buffer);
		if (!openClManager.SetKernelArg("ResolveSprings", 0, *bead_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveSprings", 1, *spring_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveSprings", 2, springs_num)) return false;
		openClManager.SetKernelGlobalWorkRange("ResolveSprings", openClManager.GetBufferLength(BufferIndex::spring_buffer));
	}

	if (angle_buffer != nullptr) {
		if (!openClManager.SetKernelArg("ResolveAngles", 0, *angle_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveAngles", 1, *spring_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveAngles", 2, *bead_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveAngles", 3, angles_num)) return false;
		openClManager.SetKernelGlobalWorkRange("ResolveAngles", openClManager.GetBufferLength(BufferIndex::angle_buffer));
	}

	if (!openClManager.SetKernelArg("ResetGrid", 0, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetGrid", 1, grid_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetGrid", openClManager.GetBufferLength(BufferIndex::grid_head_buffer));

	if (!openClManager.SetKernelArg("ResetContactList", 0, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetContactList", 1, list_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetContactList", openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));

	if (!openClManager.SetKernelArg("FillGrid", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 1, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 2, *grid_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 3, settings)) return false;   //tu tez zmiana
	if (!openClManager.SetKernelArg("FillGrid", 4, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("FillGrid", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	if (!openClManager.SetKernelArg("BuildContactList", 0, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 1, *grid_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 2, *grid_cell_nh_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 3, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 4, *bead_nh_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 5, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 6, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 7, *bead_pair_coeff_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 8, *bead_pair_table_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 9, settings)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 10, grid_size)) return false;
	openClManager.SetKernelGlobalWorkRange("BuildContactList", openClManager.GetBufferLength(BufferIndex::grid_head_buffer));

	if (!openClManager.SetKernelArg("ResolveContacts", 0, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 1, *bead_nh_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 2, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 3, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 4, *bead_pair_coeff_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 5, *bead_pair_table_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 6, settings)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 7, list_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResolveContacts", openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));
	
	if (!openClManager.SetKernelArg("AdvanceVelocity", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 1, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 2, beads_num)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 3, openClManager.GetBufferLength(BufferIndex::bead_buffer))) return false;
	openClManager.SetKernelGlobalWorkRange("AdvanceVelocity", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	if (!openClManager.SetKernelArg("ResetBeadEnergy", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetBeadEnergy", 1, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetBeadEnergy", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	return true;
};

/*
 * Uploads data from model file.
 *
 * @param input_file_contents - container holding text lines of simulation input file
 * @return bool, true if manage to perform upload, otherwise false
 */
bool SimulationManager::LoadModelData(std::vector<std::string>& input_file_content)
{
	std::string					model_file_name;
	std::vector<std::string>	model_file_content;
	std::vector<std::string>	tokens;
	int	input_file_line = 0;

	if (!openClManager.GetStatus()) {
		*log << "\tLoadModelData: ERROR - openClManager not initialized." << std::endl;
		return false;
	}

	// finds name of the model file in simulation input file
	if ((input_file_line = TextParser::ScanText(input_file_content, "model_file", tokens, 0)) > 0) {
		if (tokens.size() == 2)	{
			model_file_name = tokens[1];
		} else {
			*log << "\tLoadModelData: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'model'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadModelData: ERROR - model data file was not specified in the input file!" << std::endl;
		return false;
	}

	// reads model file content
	*log << "\tLoadModelData: reading model data from file:  " << model_file_name << std::endl;
	if (!TextParser::ReadFileContent(model_file_name, model_file_content))
	{
		*log << "\tLoadModelData: ERROR - unable to open model file: " << model_file_name << std::endl;
		system("PAUSE");
		return false;
	}

	if (model_file_content.empty())	{
		*log << "\tLoadModelData: ERROR - model file is empty." << std::endl;
		system("PAUSE");
		return false;
	}

	//uploads all structures from model file
	if (!ReadBeadsData(model_file_content)) return false;
	if (!ReadSpringData(model_file_content)) return false;
	if (!ReadAngleData(model_file_content)) return false;
	//gets initial simulation data from restart file
	if (!ReadRestartFile(input_file_content)) return false;

	*log << "\tLoadModelData: model data succesfully uploaded." << std::endl;
	return true;
};

/*
 * Uploads bond data from model file.
 *
 * @param file_contents - container holding text lines of model file
 * @return bool, true if manage to perform upload, otherwise false
 */
bool  SimulationManager::ReadBeadsData(std::vector<std::string>& file_contents) {	
	std::vector<Bead>					bead_buffer;
	std::vector<BeadInfo>				bead_info_buffer;
	std::vector<std::string>		    tokens;
	int									file_line = 0;
	
	bead_buffer.reserve(20000);
	bead_info_buffer.reserve(20000);

	//processing each line of model file containing "BEAD" token
	while ((file_line = TextParser::ScanText(file_contents, "BEAD", tokens, file_line)) != -1) {
		if (tokens[0] == "BEAD") {
			if (tokens.size() == 18) {
				BeadInfo		bead_info;
				Bead			bead;
				
				bead.position.x = std::stod(tokens[6]);
				bead.position.y = std::stod(tokens[7]);
				bead.position.z = std::stod(tokens[8]);
				
				if (bead.position.x > settings.size.x ||
					bead.position.y > settings.size.y ||
					bead.position.z > settings.size.z ||
					bead.position.x < 0.0f ||
					bead.position.y < 0.0f ||
					bead.position.z < 0.0f) {
					*log << "\tReadBeadsData: ERROR - model geometry does not fit the simulation box size." << std::endl;
					return false;
				}

				bead.velocity.x = std::stod(tokens[9]);
				bead.velocity.y = std::stod(tokens[10]);
				bead.velocity.z = std::stod(tokens[11]);
				bead.velocity_old.x = bead.velocity.x;
				bead.velocity_old.y = bead.velocity.x;
				bead.velocity_old.z = bead.velocity.x;
				bead.force.x = std::stod(tokens[12]);
				bead.force.y = std::stod(tokens[13]);
				bead.force.z = std::stod(tokens[14]);
				bead.force_old.x = bead.force.x;
				bead.force_old.y = bead.force.y;
				bead.force_old.z = bead.force.z;
				bead.image_position.x = std::stod(tokens[15]);
				bead.image_position.y = std::stod(tokens[16]);
				bead.image_position.z = std::stod(tokens[17]);

				auto it = type_id.find(tokens[1]);
				if (it == type_id.end()) {
					*log << "\tReadBeadsData: ERROR - model file line - " << file_line << "\t"
						      << " - bead type does not match any of defined types." << std::endl;
					return false; 
				}
				
				auto index = it->second;
				bead_info.type = index;
				bead_info.tag = std::stoi(tokens[2]);
				bead_info.loc_id = std::stoi(tokens[3]);
				bead_info.radius = std::stod(tokens[4]);
				bead_info.radius_sq = bead_info.radius * bead_info.radius;
				bead_info.skin_radius = bead_info.radius + skin_radius_ratio;
				bead_info.skin_radius_sq = bead_info.skin_radius * bead_info.skin_radius;
				bead_info.mass = std::stod(tokens[5]);

				if (bead_info.mass > 0.0f) {
					bead_info.mass_inv = 1.0f / bead_info.mass;
					bead_info.lambda_dt_mass_inv = settings.lambda_dt * bead_info.mass_inv;
					bead_info.half_dt_mass_inv = settings.half_dt * bead_info.mass_inv;
					bead_info.half_dt2_mass_inv = settings.half_dt2 * bead_info.mass_inv;
				} else {
					bead_info.mass_inv = 0.0f;
					bead_info.lambda_dt_mass_inv = 0.0f;
					bead_info.half_dt_mass_inv = 0.0f;
					bead_info.half_dt2_mass_inv = 0.0f;
				}
				
				bead_buffer.push_back(bead);
				bead_info_buffer.push_back(bead_info);
			} else {
				*log << "\tReadBeadsData: ERROR - invalid Particle data definition at line:  "
					      << file_line << " - solver terminates." << std::endl;
				return false;
			}
		}
	}

	beads_num = static_cast<int>(bead_buffer.size());

	if (beads_num > 0) {	
		// basic buffer of Bead structures 
		if ((BufferIndex::bead_buffer = openClManager.CreateBuffer<Bead>(bead_buffer, CL_MEM_READ_WRITE)) < 0) {
			*log << "\tReadBeadsData: ERROR - unable to load BufferIndex::particles !" << std::endl;
			return false;
		}
	
		// buffer of BeadInfo structures, holding bead additional data and precalculated parameteres
		if ((BufferIndex::bead_info_buffer = openClManager.CreateBuffer<BeadInfo>(bead_info_buffer, CL_MEM_READ_ONLY)) < 0) {
			*log << "\tReadBeadsData: ERROR - unable to load BufferIndex::particles_info !" << std::endl;
			return false;
		}

		// buffer of holding nmubers of 
		if ((BufferIndex::bead_nh_head_buffer = openClManager.CreateBuffer<int>(beads_num, CL_MEM_READ_WRITE)) < 0) {
			*log << "\tReadBeadsData: ERROR - unable to load BufferIndex::listHead !" << std::endl;
			return false;
		}
				
		if ((BufferIndex::bead_nh_tail_buffer = openClManager.CreateBuffer<ParticleNh>(beads_num, CL_MEM_READ_WRITE)) < 0) {
			*log << "\tReadBeadsData: ERROR - unable to load BufferIndex::listTail  !" << std::endl;
			return false;
		}

		if ((BufferIndex::bead_pair_coeff_buffer = openClManager.CreateBuffer<ParticleNh>(beads_num, CL_MEM_READ_WRITE)) < 0) {
			*log << "\tReadBeadsData: ERROR - unable to load BufferIndex::intercationCoeffs !" << std::endl;
			return false;
		}
		*log << "\tReadBeadsData: Bead data uploaded succesfully." << std::endl;
		return true;
	} else {
		*log << "\tReadBeadsData: ERROR - model data requires bead definition." << std::endl;
		return false;
	}
}

/*
 * Uploads harmonic bond data from model file.
 *
 * @param file_contents - container holding text lines of model file
 * @return bool, true if manage to perform upload, otherwise false
 */
bool  SimulationManager::ReadSpringData(std::vector<std::string>& file_contents) {
	int								file_line = 0;
	std::vector<SpringStructure>	spring_data;
	std::vector<std::string>		tokens;
	
	spring_data.reserve(2000);

	//processing each line of model file containing "SPRING" token
	while ((file_line = TextParser::ScanText(file_contents, "SPRING", tokens, file_line)) != -1) {
		if (tokens[0] == "SPRING") {
			if (tokens[1] == "HR_FL" && tokens.size() != 9)	{
				*log << "\tReadSpringData: ERROR - invalid spring data definition at line - "
					<< file_line << " - solver terminates." << std::endl;
				return false;
			}

			if (tokens[1] == "HR" && tokens.size() < 7)	{
				*log << "\tReadSpringData: ERROR - invalid spring data definition at line - "
					<< file_line << " - solver terminates." << std::endl;
				return false;
			}

			SpringStructure spring;
			if (tokens[1] == "HR") {
				spring.type = 0;
			} else if (tokens[1] == "HR_FL") {
				spring.type = 1;
				spring.c1 = std::stoi(tokens[7]);
				spring.c2 = std::stoi(tokens[8]);
			} else {
				*log << "Error: invalid spring data definition at line - "
					<< file_line << " - solver terminates." << std::endl;
				return false;
			}

			spring.tag = std::stoi(tokens[2]);

			spring.p1 = std::stoi(tokens[3]);
			if (spring.p1 < 0 && spring.p1 >= beads_num) {
				*log << "\tReadSpringData: ERROR - model file line - " << file_line << "\t"
					<< " - wrong bead index reference." << std::endl;
				return false;
			}

			spring.p2 = std::stoi(tokens[4]);
			if (spring.p2 < 0 && spring.p2 >= beads_num) {
				*log << "\tReadSpringData: ERROR - model file line - " << file_line << "\t"
					<< " - wrong bead index reference." << std::endl;
				return false;
			}

			spring.rest_length = std::stod(tokens[5]);
			if (spring.rest_length < 0.0f) {
				*log << "\tReadSpringData: ERROR - model file line - " << file_line << "\t"
					<< " - wrong rest length value." << std::endl;
				return false;
			}

			spring.stiffness = std::stod(tokens[6]);
			if (spring.stiffness < 0.0f) {
				*log << "\tReadSpringData: ERROR - model file line - " << file_line << "\t"
					<< " - wrong stiffness value." << std::endl;
				return false;
			}
			spring.status = 1;
			spring_data.push_back(spring);
		}		
	}

	if (spring_data.size() > 0) {
		springs_num = static_cast<int>(spring_data.size());
		BufferIndex::spring_buffer = openClManager.CreateBuffer<SpringStructure>(spring_data, CL_MEM_READ_WRITE);
		if (BufferIndex::spring_buffer >= 0) {
			*log << "\tReadSpringData: spring bonds data uploaded succesfully." << std::endl;
			return true;
		} else {
			*log << "\tReadSpringData: ERROR - unable to upload spring bonds data." << std::endl;
			return false;
		}
	} else {
		return true;
	}
};

/*
 * Uploads angle bond data from model file.
 *
 * @param file_contents - container holding text lines of model file
 * @return bool, true if manage to perform upload, otherwise false
 */
bool  SimulationManager::ReadAngleData(std::vector<std::string>& file_contents) {
	int								file_line = 0;
	std::vector<AngleBond>			angle_data;
	std::vector<std::string>		tokens;
	angle_data.reserve(2000);

	//processing each line of model file containing "ANGLE" token
	while ((file_line = TextParser::ScanText(file_contents, "ANGLE", tokens, file_line)) != -1) {
		if (tokens[0] == "ANGLE") {
			if (tokens.size() == 10) {
				
				AngleBond angle_bond;
				if (tokens[1] == "HR") {
					angle_bond.type = 0;
				} else if (tokens[1] == "HR_FL") {
					angle_bond.type = 1;
				} else {
					*log << "\tReadAngleData: ERROR - invalid able bond data definition at line - "
						<< file_line << " - solver terminates." << std::endl;
					return false;
				}

				angle_bond.tag = std::stoi(tokens[2]);
				angle_bond.b1 = std::stoi(tokens[3]);
				if (angle_bond.b1 < 0 && angle_bond.b1 >= beads_num) {
					*log << "\tReadAngleData: ERROR - model file line - " << file_line << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}

				angle_bond.b2 = std::stoi(tokens[4]);
				if (angle_bond.b2 < 0 && angle_bond.b2 >= beads_num) {
					*log << "\tReadAngleData: ERROR - model file line - " << file_line << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}
				
				angle_bond.b3 = std::stoi(tokens[5]);
				if (angle_bond.b3 < 0 && angle_bond.b3 >= beads_num) {
					*log << "\tReadAngleData: ERROR - model file line - " << file_line << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}

				angle_bond.s1 = std::stoi(tokens[6]);
				angle_bond.s2 = std::stoi(tokens[7]);
				angle_bond.angle = std::stof(tokens[8]);
				angle_bond.c1 = std::stod(tokens[9]);
				angle_bond.status = 1;
				angle_data.push_back(angle_bond);
			} else {
				*log << "\tReadAngleData: ERROR - invalid angle data definition at line - "
					<< file_line << " - solver terminates." << std::endl;
				return false;
			}
		}
	}
	if (angle_data.size() > 0) {
		angles_num = static_cast<int>(angle_data.size());
		BufferIndex::angle_buffer = openClManager.CreateBuffer<AngleBond>(angle_data, CL_MEM_READ_ONLY);
		if (BufferIndex::angle_buffer >= 0)	{
			*log << "\tReadAngleData: angle bonds data uploaded succesfully." << std::endl;
			return true;
		} else {
			*log << "\tReadAngleData: ERROR - unable to upload angle bonds data." << std::endl;
			return false;
		}
	} else {
		return true;
	}
};

/*
 * Uploads bead data from restart file based on informations
 * from input file.
 *
 * @param input_file_content - container holding text lines of input file
 * @return bool, true if manage to perform upload, otherwise false
 */
bool SimulationManager::ReadRestartFile(std::vector<std::string>& input_file_content)
{
	int input_file_line = 0;
	std::string					restart_file_name;
	std::vector<std::string>	restart_file_content;
	std::vector<std::string>	tokens;
		
	bool restart_position;
	bool restart_velocity;
	bool restart_force;

	restart_position = false;
	restart_velocity = false;
	restart_force = false;

	if (!openClManager.GetStatus()) {
		*log << "\tReadRestartFile: ERROR - openClManager not initialized." << std::endl;
		return false;
	}
	
	int my_start = 0;
	while ((my_start = TextParser::ScanText(input_file_content, "restart_file", tokens, my_start)) != -1) {
		if (tokens.size() == 2) {
			restart_file_name = tokens[1];
		} else {
			*log << "\tReadRestartFile: ERROR - line " << my_start << " - incorrect number of command parameters - 'restart'" << std::endl;
			return false;
		}
	}
	
	my_start = 0;
	while ((my_start = TextParser::ScanText(input_file_content, "restart_data", tokens, my_start)) != -1) {
		if (tokens.size() == 2)	{
			if (tokens[1] == "position") {
				restart_position = true;
			} else if (tokens[1] == "velocity") {
				restart_velocity = true;
			} else if (tokens[1] == "force") {
				restart_force = true;
			} else if (tokens[1] == "all") {
				restart_position = true;
				restart_velocity = true;
				restart_force = true;
			} else {
				*log << "\tReadRestartFile: ERROR - line " << my_start << " - incorrect definition of command parameter - 'restart data'" << std::endl;
				return false;
			}
		} else {
			*log << "\tReadRestartFile: ERROR - line " << my_start << " - incorrect number of command parameters - 'restart'" << std::endl;
			return false;
		}
	}


	if (restart_file_name.empty()) {
		*log << "\tReadRestartFile: restart data file was not specified." << std::endl;
		return true;
	}


	*log << "\tReadRestartFile: reading restart data from file:  " << restart_file_name << std::endl;
	if (!TextParser::ReadFileContent(restart_file_name, restart_file_content)) {
		*log << "\tReadRestartFile: ERROR - unable to open restart file: " << restart_file_name << std::endl;
		system("PAUSE");
		return false;
	}
	
	if (!restart_file_content.empty()) {
		std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
		if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) {
			*log << "\tReadRestartFile: ERROR - unable to read BufferIndex::bead_buffer data." << std::endl;
			return false;
		}

		std::vector<BeadInfo> bead_info_data(openClManager.GetBufferLength(BufferIndex::bead_info_buffer), BeadInfo());
		if (!openClManager.ReadBuffer(BufferIndex::bead_info_buffer, &bead_info_data.front())) {
			*log << "\tReadRestartFile: ERROR - unable to read BufferIndex::bead_info_buffer data." << std::endl;
			return false;
		}

		int number_of_beads = openClManager.GetBufferElements(BufferIndex::bead_buffer);
		int bead_index = 0;

		for (int ii = 0; ii < static_cast<int>(restart_file_content.size()); ii++) {
			auto str = restart_file_content[ii];
			std::istringstream iss(str);
			std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
			
			if (!tokens.empty()) {
				if (tokens.size() == 13) {
					if (bead_index < number_of_beads) {
						if (type_id.find(tokens[0]) != type_id.end()) {
							int this_type = type_id.find(tokens[0])->second;
							if (bead_info_data[bead_index].type != this_type) {
								*log << "\tReadRestartFile: ERROR - the type of bead from restart file does not match the type of bead id = " << bead_index << std::endl;
								return false;
							}
						} else {
							*log << "\tReadRestartFile: ERROR - unable to recognize the type of bead from restart file." << std::endl;
							return false;
						}

						if (restart_position == true) {
							bead_data[bead_index].position.x = std::stod(tokens[1]);
							bead_data[bead_index].position.y = std::stod(tokens[2]);
							bead_data[bead_index].position.z = std::stod(tokens[3]);
							bead_data[bead_index].image_position.x = std::stod(tokens[4]);
							bead_data[bead_index].image_position.y = std::stod(tokens[5]);
							bead_data[bead_index].image_position.z = std::stod(tokens[6]);
						}

						if (restart_velocity == true) {
							bead_data[bead_index].velocity.x = std::stod(tokens[7]);
							bead_data[bead_index].velocity.y = std::stod(tokens[8]);
							bead_data[bead_index].velocity.z = std::stod(tokens[9]);
							bead_data[bead_index].velocity_old.x = std::stod(tokens[7]); 
							bead_data[bead_index].velocity_old.y = std::stod(tokens[8]);
							bead_data[bead_index].velocity_old.z = std::stod(tokens[9]);
						}

						if (restart_force == true) {
							bead_data[bead_index].force.x = std::stod(tokens[10]);
							bead_data[bead_index].force.y = std::stod(tokens[11]);
							bead_data[bead_index].force.z = std::stod(tokens[12]);
							bead_data[bead_index].force_old.x = 0.0f;
							bead_data[bead_index].force_old.y = 0.0f;
							bead_data[bead_index].force_old.z = 0.0f;
						}
					} else {
						*log << "\tReadRestartFile: ERROR - the number of bead from restart file is larger than number of bead in model file." << std::endl;
						return false;
					}
					bead_index++;
				}
			}
			tokens.clear();
		}
		
		if (bead_index != number_of_beads) {
			*log << "\tReadRestartFile: ERROR - incorrect number of uploaded beads." << std::endl;
			return false;
		}
		
		if (!openClManager.WriteBuffer<Bead>(BufferIndex::bead_buffer, bead_data)) {
			*log << "\tReadRestartFile: ERROR - unable to write data into BufferIndex::bead_buffer." << std::endl;
			return false;
		}
		simulation_restart = true;
		return true;
	} else {
		*log << "\tReadRestartFile: ERROR - restart file: " << restart_file_name <<  " is empty." <<std::endl;
		return false;
	}
};

/*
 * Creates data dump text file.
 *
 * @param value - strin containing name of the data dump file.
 * @return bool, true if manage to create data dump file, otherwise false
 */
bool SimulationManager::SetDumpFile(const std::string& value) {
	std::cout << "\tOutput file: " << value << std::endl;
	if (std::ifstream(value)) {
		*log << "\tWarning - output file already exists - overwriting existing file" << std::endl;
	}
	
	std::fstream dump_file;	
	dump_file.open(value, std::ios::out);
	if (dump_file.is_open()) {
		if (dump_file.good()) {
			dump_file.close();
			return true;
		}
	}
	*log << "Error: Unable to set output file!" << std::endl;
	return false;
};

/*
 * Performs a full simulation data dump into restart file.
 *
 * @param progress - integer indicating current simulation step.
 * @return bool, true if manage to perform data dump, otherwise false
 */
bool SimulationManager::DumpRestartData(int progress) {
	if (!openClManager.GetStatus())	{
		*log << "\tDumpRestartData: ERROR - incorrect number of uploaded beads." << std::endl;
		return false;
	}
	
	if (BufferIndex::bead_buffer < 0) {
		*log << "\tDumpRestartData: ERROR - bead buffer not initialized." << std::endl;
		return false;
	}
	
	if (BufferIndex::bead_info_buffer < 0) {
		*log << "\tDumpRestartData: ERROR - bead info buffer not initialized." << std::endl;
		return false;
	}

	if (openClManager.GetBufferLength(BufferIndex::bead_buffer) == -1) {
		*log << "\tDumpRestartData: ERROR - bead buffer is empty." << std::endl;
		return false;
	}
	
	if (openClManager.GetBufferLength(BufferIndex::bead_info_buffer) == -1) {
		*log << "\tDumpRestartData: ERROR - bead info buffer is empty." << std::endl;
		return false;
	}

	std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
	std::vector<BeadInfo> bead_info(openClManager.GetBufferLength(BufferIndex::bead_info_buffer), BeadInfo());

	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) {
		*log << "\tDumpRestartData: ERROR - unagle to read bead buffer data." << std::endl;
		return false;
	}

	if (!openClManager.ReadBuffer(BufferIndex::bead_info_buffer, &bead_info.front())) {
		*log << "\tDumpRestartData: ERROR - unagle to read bead info buffer data." << std::endl;
		return false;
	}

	std::fstream dump_file(restart_file_name, std::fstream::out | std::fstream::trunc);
	if (!dump_file.is_open()) {
		*log << "\tDumpRestartData: ERROR - unable to open restart data dump file." << std::endl;
		return false;
	}

	if (!dump_file.good()) {
		*log << "\tDumpRestartData: ERROR - restart data dump file corrupted." << std::endl;
		return false;
	}

	int number_of_beads = openClManager.GetBufferElements(BufferIndex::bead_buffer);
	dump_file << number_of_beads << std::endl;
	dump_file << "FRAME" << "\t" << progress << "\t" << settings.dt << "\t" << settings.size.x << "\t" << settings.size.y << "\t" << settings.size.z << std::endl;
	for (auto ii = 0; ii < number_of_beads; ii++) {
		auto id = bead_info[ii].type;
		auto type = std::find_if(type_id.begin(), type_id.end(), [&id](const std::pair<std::string, int> &p){return p.second == id; });

		dump_file << "\t" << type->first;
		dump_file << "\t" << bead_data[ii].position.x;
		dump_file << "\t" << bead_data[ii].position.y;
		dump_file << "\t" << bead_data[ii].position.z;
		dump_file << "\t" << bead_data[ii].image_position.x;
		dump_file << "\t" << bead_data[ii].image_position.y;
		dump_file << "\t" << bead_data[ii].image_position.z;
		dump_file << "\t" << bead_data[ii].velocity.x;
		dump_file << "\t" << bead_data[ii].velocity.y;
		dump_file << "\t" << bead_data[ii].velocity.z;
		dump_file << "\t" << bead_data[ii].force.x;
		dump_file << "\t" << bead_data[ii].force.y;
		dump_file << "\t" << bead_data[ii].force.z << std::endl;
	}
	*log << "\tDumpRestartData: restart data dump successful." << std::endl;
	return true;
}

/*
 * Uploads bead data from OpenCL buffer and writes selected
 * data into output data files.
 *
 * @param progress - integer indicating current simulation step.
 * @return bool, true if manage to perform data dump, otherwise false
 */
bool SimulationManager::DumpData(int progress)
{
	//static variables for data output
	static int dump_counter = 1;
	static float restart_dump_stage = 0.25f;
	static std::vector<Bead> bead_data; 
	static std::fstream dump_file;
	static std::fstream vel_dump_file;
	static std::fstream force_dump_file;

	if (!openClManager.GetStatus())	{
		*log << "\tDumpData: ERROR - openClManager is not initialized." << std::endl;
		return false;
	}
	if (BufferIndex::bead_buffer < 0) {
		*log << "\tDumpData: ERROR - BufferIndex::bead_buffer is not initialized." << std::endl;
		return false;
	}
	if (export_bead_index.empty()) {
		*log << "\tDumpData: ERROR - export_bead_index is empty." << std::endl;
		return false;
	}
	
	// initialize static variables at first data dump
	if (dump_counter == 1) {
		
		bead_data = std::vector<Bead>(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
		if (bead_data.empty()) {
			*log << "\tDumpData: ERROR - unable to initialize bead_data buffer." << std::endl;
			return false;
		}

		if (export_position == true) {
			dump_file = std::fstream(dump_file_name, std::ios_base::app);
			if(!dump_file.good()) {
				*log << "\tDumpData: ERROR - unable to initialize dump_file." << std::endl;
			}
		}
		
		if (export_velocity == true) {
			vel_dump_file = std::fstream(vel_dump_file_name, std::ios_base::app);
			if (!vel_dump_file.good()) {
				*log << "\tDumpData: ERROR - unable to initialize vel_dump_file." << std::endl;
			}
		}

		if (export_force == true) {
			force_dump_file = std::fstream(force_dump_file_name, std::ios_base::app);
			if (!force_dump_file.good()) {
				*log << "\tDumpData: ERROR - unable to initialize force_dump_file." << std::endl;
			}
		}

		dump_file << std::fixed << std::setprecision(5);
		vel_dump_file << std::fixed << std::setprecision(5);
		force_dump_file << std::fixed << std::setprecision(5);
	}
	
	// read bead data from openCL memory buffer
	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) {
		*log << "\tDumpData: ERROR - unable to read data from BufferIndex::bead_buffer." << std::endl;
		return false;
	}

	if (export_position == true) {
		if (dump_file.is_open() && dump_file.good()) {
			dump_file << export_bead_index.size() << std::endl;
			dump_file << "FRAME" << "\t" << progress << "\t" << settings.dt << "\t" << settings.size.x << "\t" << settings.size.y << "\t" << settings.size.z << std::endl;
			for (auto ii = 0; ii < export_bead_index.size(); ii++) {
				dump_file << export_bead_name[ii];
				dump_file << "\t" << bead_data[export_bead_index[ii]].position.x;
				dump_file << "\t" << bead_data[export_bead_index[ii]].position.y;
				dump_file << "\t" << bead_data[export_bead_index[ii]].position.z << std::endl;
			}
		} else {
			*log << "\tDumpData: ERROR - unable to write data to dump_file" << std::endl;
			return false;
		}
	}

	if (export_velocity == true) {
		if (vel_dump_file.good() && vel_dump_file.is_open()) {
			vel_dump_file << export_bead_index.size() << std::endl;
			vel_dump_file << "FRAME" << "\t" << progress << "\t" << settings.dt << "\t" << settings.size.x << "\t" << settings.size.y << "\t" << settings.size.z << std::endl;
			
			for (auto ii = 0; ii < export_bead_index.size(); ii++) {
				vel_dump_file << export_bead_name[ii];
				vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.x;
				vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.y;
				vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.z << std::endl;
			}
		} else {
			*log << "\tDumpData: ERROR - unable to write data to vel_dump_file" << std::endl;
			return false;
		}
	}

	if (export_force == true) {
		if (force_dump_file.good() && force_dump_file.is_open()) {
			force_dump_file << export_bead_index.size() << std::endl;
			force_dump_file << "FRAME" << "\t" << progress << "\t" << settings.dt << "\t" << settings.size.x << "\t" << settings.size.y << "\t" << settings.size.z << std::endl;

			for (auto ii = 0; ii < export_bead_index.size(); ii++) {
				force_dump_file << export_bead_name[ii];
				force_dump_file << "\t" << bead_data[export_bead_index[ii]].force.x;
				force_dump_file << "\t" << bead_data[export_bead_index[ii]].force.y;
				force_dump_file << "\t" << bead_data[export_bead_index[ii]].force.z << std::endl;
			}
		} else {
			*log << "\tDumpData: ERROR - unable to write data to force_dump_file" << std::endl;
			return false;
		}
	}

	if ((static_cast<float>(progress) / static_cast<float>(simulation_steps)) > restart_dump_stage) {
		restart_dump_stage += 0.25f;
		if (!DumpRestartData(progress)) {
			*log << "\tDumpData: ERROR - unable to dump restart data at step " << progress << std::endl;
			return false;
		}
	}

	dump_counter++;
	return true;
};

/*
 * Loads input parametes for simulation box creates predefined 
 * data structures and loads them into the openCL memory buffer.
 *
 * @param file_contents - container of input file text lines.
 * @return bool, true if manage to create a simulation box, otherwise false
 */
bool SimulationManager::LoadSimulationBoxInfo(std::vector<std::string>& file_contents)
{
	int	input_file_line = 0;
	std::vector<std::string>	tokens;
	
	//reads box_unit_cell size from input file
	if ((input_file_line = TextParser::ScanText(file_contents, "box_unit_cell", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			settings.box_unit_cell = std::stod(tokens[1]);
			if (settings.box_unit_cell <= 0.0f) {
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect value of command parameter: 'box_unit_cell'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'box_unit_cell'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - box_unit_cell was not specified in the input file." << std::endl;
		return false;
	}
	
	//reads box_size size from input file
	if ((input_file_line = TextParser::ScanText(file_contents, "box_size", tokens, 0)) > 0) {
		if (tokens.size() == 4)	{
			settings.size.x = std::stod(tokens[1]);
			settings.size.y = std::stod(tokens[2]);
			settings.size.z = std::stod(tokens[3]);

			if (settings.size.x < 0.0f || settings.size.y < 0.0f || settings.size.z < 0.0f) {
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - box dimensions must have positive values." << std::endl;
				return false;
			}

			if ((settings.size.x / settings.box_unit_cell) < 3.0f ||
				(settings.size.y / settings.box_unit_cell) < 3.0f ||
				(settings.size.z / settings.box_unit_cell) < 3.0f) 	{
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - box dimensions must be at least three times bigger than cutoff radius." << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'box_size'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - box dimensions were not specified in the input file." << std::endl;
		return false;
	}

	//reads box_type size from input file
	if ((input_file_line = TextParser::ScanText(file_contents, "box_type", tokens, 0)) > 0) {
		if (tokens.size() == 2) {
			if (tokens[1] == "pbc") {
				settings.bc_type = 0;
			} else if (tokens[1] == "finite") {
				settings.bc_type = 1;
			} else {
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect parameter value: 'box_type'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'box_type'" << std::endl;
			return false;
		}
	} else {
		*log << "\tLoadSimulationBoxInfo: ERROR - 'box_type' not specified in input file." << std::endl;
		return false;
	}

	//reads grid_update_freq size from input file
	if ((input_file_line = TextParser::ScanText(file_contents, "grid_update_freq", tokens, 0)) > 0) {
		if (tokens.size() == 2)	{
			grid_update_freq = std::stoi(tokens[1]);
			if (grid_update_freq <= 0) {
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect value of command parameter: 'grid_update_freq'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'grid_update_freq'" << std::endl;
			return false;
		}
	} else {
		grid_update_freq = 8;
	}

	//reads box_wall_rep size from input file
	if ((input_file_line = TextParser::ScanText(file_contents, "box_wall_rep", tokens, 0)) > 0) {
		if (tokens.size() == 2)	{
			settings.wall_contact_stiffness = std::stod(tokens[1]);
			if (settings.wall_contact_stiffness < 0.0f)	{
				*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect value of command parameter: 'box_wall_rep'" << std::endl;
				return false;
			}
		} else {
			*log << "\tLoadSimulationBoxInfo: ERROR - line " << input_file_line << " - incorrect number of command parameters: 'box_wall_rep'" << std::endl;
			return false;
		}
	} else {
		settings.bc_type = 1;
		*log << "\tLoadSimulationBoxInfo: specified box type reqires to set value for box_wall_rep" << std::endl;
	}
	
	if (settings.size.x <= 0.0f || settings.size.y <= 0.0f || settings.size.z <= 0.0f) {
		*log << "\tLoadSimulationBoxInfo: ERROR - invalid grid size." << std::endl;
		return false;
	}

	std::vector<GridCellNh>	grid_cell_nhood;
	settings.x_cells_num = static_cast<int>(floor(settings.size.x / settings.box_unit_cell));
	settings.y_cells_num = static_cast<int>(floor(settings.size.y / settings.box_unit_cell));
	settings.z_cells_num = static_cast<int>(floor(settings.size.z / settings.box_unit_cell));
	settings.xy_cells_num = settings.x_cells_num * settings.y_cells_num;

	settings.step.x = settings.size.x / static_cast<float>(settings.x_cells_num);
	settings.step.y = settings.size.y / static_cast<float>(settings.y_cells_num);
	settings.step.z = settings.size.z / static_cast<float>(settings.z_cells_num);

	settings.half_size_squared = std::min<double>(std::min<double>(settings.size.x, settings.size.y), settings.size.z)*0.5f;
	settings.half_size_squared *= settings.half_size_squared;

	int grid_size = settings.x_cells_num * settings.y_cells_num * settings.z_cells_num;

	//initialize collsion grid array with default cells
	GridCellNh test;
	test.n[0] = -1;
	test.n[1] = -1;
	test.n[2] = -1;
	test.n[3] = -1;
	test.n[4] = -1;
	test.n[5] = -1;
	test.n[6] = -1;
	test.n[7] = -1;
	test.n[8] = -1;
	test.n[9] = -1;
	test.n[10] = -1;
	test.n[11] = -1;
	test.n[12] = -1;
	grid_cell_nhood.resize(grid_size, test);

	//foe each grid cell set adresses of neigbouring cells
	for (int zz = 0; zz < settings.z_cells_num; zz++) {
		for (int yy = 0; yy < settings.y_cells_num; yy++) {
			for (int xx = 0; xx < settings.x_cells_num; xx++) {
				int targetCell = xx + yy * settings.x_cells_num + zz * settings.x_cells_num * settings.y_cells_num;
				int cnt = 0;
				for (int dx = -1; dx < 2; dx++) {
					for (int dy = -1; dy < 2; dy++) {
						int nx = xx + dx;
						int ny = yy + dy;
						int nz = zz - 1;

						if (settings.bc_type == 0) {
							if (nx < 0) nx = settings.x_cells_num - 1;
							if (nx >= settings.x_cells_num) nx = 0;
							if (ny < 0) ny = settings.y_cells_num - 1;
							if (ny >= settings.y_cells_num) ny = 0;
							if (nz < 0) nz = settings.z_cells_num - 1;
							if (nz >= settings.z_cells_num) nz = 0;

							int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
							grid_cell_nhood[targetCell].n[cnt] = nh;
						}

						if (settings.bc_type == 1) {
							if (nx < 0 || nx >= settings.x_cells_num || ny < 0 || ny >= settings.y_cells_num || nz < 0 || nz >= settings.z_cells_num) {
								grid_cell_nhood[targetCell].n[cnt] = -1;
							}
							else {
								int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
								grid_cell_nhood[targetCell].n[cnt] = nh;
							}
						}
						cnt++;
					}
				}

				for (int dy = -1; dy < 2; dy++) {
					int nx = xx + 1;
					int ny = yy + dy;
					int nz = zz;

					if (settings.bc_type == 0) {
						if (nx < 0) nx = settings.x_cells_num - 1;
						if (nx >= settings.x_cells_num) nx = 0;
						if (ny < 0) ny = settings.y_cells_num - 1;
						if (ny >= settings.y_cells_num) ny = 0;
						if (nz < 0) nz = settings.z_cells_num - 1;
						if (nz >= settings.z_cells_num) nz = 0;

						int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
						grid_cell_nhood[targetCell].n[cnt] = nh;
					}

					if (settings.bc_type == 1) {
						if (nx < 0 || nx >= settings.x_cells_num || ny < 0 || ny >= settings.y_cells_num || nz < 0 || nz >= settings.z_cells_num) {
							grid_cell_nhood[targetCell].n[cnt] = -1;
						}
						else {
							int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
							grid_cell_nhood[targetCell].n[cnt] = nh;
						}
					}
					cnt++;
				}

				int nx = xx;
				int ny = yy - 1;
				int nz = zz;
				if (settings.bc_type == 0) {
					if (nx < 0) nx = settings.x_cells_num - 1;
					if (nx >= settings.x_cells_num) nx = 0;
					if (ny < 0) ny = settings.y_cells_num - 1;
					if (ny >= settings.y_cells_num) ny = 0;
					if (nz < 0) nz = settings.z_cells_num - 1;
					if (nz >= settings.z_cells_num) nz = 0;

					int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
					grid_cell_nhood[targetCell].n[cnt] = nh;
				}
				if (settings.bc_type == 1) {
					if (nx < 0 || nx >= settings.x_cells_num || ny < 0 || ny >= settings.y_cells_num || nz < 0 || nz >= settings.z_cells_num) {
						grid_cell_nhood[targetCell].n[cnt] = -1;
					}
					else {
						int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
						grid_cell_nhood[targetCell].n[cnt] = nh;
					}
				}
				cnt++;
			}
		}
	}
	// head buffer of collision grid; each element corresponds to one grid cell; 
	// numerical value indicates the number of beads inside cell
	BufferIndex::grid_head_buffer = openClManager.CreateBuffer<int>(grid_size, CL_MEM_READ_WRITE);
	if (BufferIndex::grid_head_buffer < 0) {
		*log << "\tLoadSimulationBoxInfo: ERROR - unable to load BufferIndex::gridHead." << std::endl;
		return false;
	}

	// tail buffer of collision grid cells; each element contains vector of 
	// bead ID's inside specific cell
	BufferIndex::grid_tail_buffer = openClManager.CreateBuffer<GridCell>(grid_size, CL_MEM_READ_WRITE);
	if (BufferIndex::grid_tail_buffer < 0) {
		*log << "\tLoadSimulationBoxInfo: ERROR - unable to load BufferIndex::gridTail." << std::endl;
		return false;
	}
	// buffer of structures containing addresses of adjacent cells in the collision grid
	// for each individial grid cell
	BufferIndex::grid_cell_nh_buffer = openClManager.CreateBuffer<GridCellNh>(grid_cell_nhood, CL_MEM_READ_ONLY);
	if (BufferIndex::grid_cell_nh_buffer < 0) {
		*log << "\tLoadSimulationBoxInfo: ERROR - unable to load BufferIndex::nHoodBuffer." << std::endl;
		return false;
	}
	return true;
};

void SimulationManager::Close() {
};
