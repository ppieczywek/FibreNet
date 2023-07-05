#pragma once
#include "stdafx.h"
#include "OpenclManager.h"

bool OpenclManager::Initialize(std::shared_ptr<log_stream> p) {
	if (p == nullptr) return false;
	log = p;

	localWorkGroupSize = 64;
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);

	if (all_platforms.size() == 0) {
		*log << "\tOpenClManager:" << "\t\t" << "No platforms found. Check OpenCL installation!\n";
		return false;
	}
	*log << "\tSelect available OpenCL platform:" << std::endl << std::endl;

	for (int ii = 0; ii < all_platforms.size(); ii++) {
		*log << "\t[" << ii + 1 << "]" << "\t" << all_platforms[ii].getInfo<CL_PLATFORM_NAME>() << std::endl;
	}
	*log << "\t[" << all_platforms.size()+1 << "]" << "\t" << "terminate program" << std::endl;
		
	int input = 0;
	while (1) {
		*log << std::endl << "\tType platform ID: ";
		std::cin >> input;
		*log << std::endl;

		if ((all_platforms.size() + 1) == input) return false;

		if ((input-1) >= 0 && (input-1) < all_platforms.size()) {
			default_platform = all_platforms[input-1];
			*log << "\tSelected platform:\t" << default_platform.getInfo<CL_PLATFORM_NAME>() << std::endl << std::endl;;
			break;
		}
	}
	
	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
	if (all_devices.size() == 0) {
		*log << "\tOpenClManager:" << "\t\t" << " No GPU devices found !\n";
		default_platform.getDevices(CL_DEVICE_TYPE_CPU, &all_devices);
		if (all_devices.size() == 0) return false;
	}

	*log << std::endl << std::endl << "\tDevices available on selected platform: \n\n";
	for (int ii = 0; ii < all_devices.size(); ii++) {
		bool device_available = false;
		if (all_devices[ii].getInfo<CL_DEVICE_AVAILABLE>() == CL_TRUE) device_available = true;
		
		*log << "\t[" << ii + 1 << "]" << "\t" << all_devices[ii].getInfo<CL_DEVICE_NAME>() << std::endl;
		if (device_available == true) {
			*log << "\t\t\t global memory size: " << "\t\t" <<
						 all_devices[ii].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / static_cast<int64_t>(1024*1024) << " MB" << std::endl;
			
			*log << "\t\t\t local memory size: " << "\t\t" <<
						 all_devices[ii].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() / static_cast<int64_t>(1024) << " Bytes" << std::endl;

			*log << "\t\t\t max work group size: " << "\t\t" <<
						 all_devices[ii].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << std::endl << std::endl;
		} else {
			*log << "\t - device not available" << std::endl;
		}
	}
	*log << "\t[" << all_devices.size() + 1 << "]" << "\t" << "terminate program" << std::endl;
	input = 0;
	while (1) {
		*log << std::endl << "\tType device ID: ";
		std::cin >> input;
		*log << std::endl << std::endl;

		if ((all_devices.size() + 1) == input) return false;
		if ((input - 1) >= 0 && (input - 1) < all_devices.size()) {
			default_device = all_devices[input - 1];
			*log << "\tSelected device:\t" << default_device.getInfo<CL_DEVICE_NAME>() << "\r";
			break;
		}
	}

	context = cl::Context({ default_device });
	queue = cl::CommandQueue(context, default_device);
	managerStatus = true;
	return true;
};

bool OpenclManager::RegisterKernel(const std::vector<std::string> &kernelNames, const std::string &kernelCode) {
	if (!managerStatus) *log << "\tOpenclManager::RegisterKernel: ERROR - manager not initialized" << std::endl;
	cl::Program::Sources sources;
	sources.push_back({ kernelCode.c_str(), kernelCode.length() });
	cl::Program program(context, sources);

	std::string options = "-cl-mad-enable";
	if (program.build({ default_device },options.c_str(), NULL,NULL) != CL_SUCCESS) {
		*log << "\tOpenclManager::RegisterKernel: ERROR - unable to build: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
		return false;
	}

	for (auto ii = 0; ii < kernelNames.size(); ii++) {
		cl_int err;
		cl::Kernel my_kernel(program, kernelNames[ii].c_str(), &err);
		if (err == CL_SUCCESS) {
			kernels.push_back(my_kernel);
			kernelsNames.push_back(kernelNames[ii]);
			kernelGlobalWorkRange.push_back(0);
		} else {
			*log << "\tOpenclManager::RegisterKernel: ERROR - unable to crete kernel:" << kernelNames[ii] << "\n";
			checkErr(err);
			return false;
		}
	}
	return true;
};

cl::Buffer* OpenclManager::GetBuffer(const int& bufferId) {
	if (!clBuffer.empty()) {
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) return &clBuffer[bufferId];
		}
	}
	return nullptr;
};

int OpenclManager::GetLocalWorkGroupSize() {
	return localWorkGroupSize;
}

int OpenclManager::GetBufferSize(const int& bufferId) {
	if (!clBuffer.empty()) {
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) return clBufferSize[bufferId];
		}
	}
	return -1;
};

int OpenclManager::GetBufferElements(const int& bufferId) {
	if (!clBuffer.empty()) {
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) return clBufferElements[bufferId];
		}
	}
	return -1;
};

int OpenclManager::GetBufferLength(const int& bufferId) {
	if (!clBuffer.empty()) {
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) return clBufferLength[bufferId];	
		}
	}
	return -1;
};

int OpenclManager::GetBufferStructureSize(const int& bufferId) {
	if (!clBuffer.empty()) 	{
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) return clBufferStructureSize[bufferId];
		}
	}
	return -1;
};

bool OpenclManager::ReadBuffer(int bufferId, void* ptr) {
	if (!clBuffer.empty()) {
		if (bufferId < clBuffer.size()) {
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0) {
				cl_int err = queue.enqueueReadBuffer(clBuffer[bufferId], 
													 CL_TRUE,
													 0,
													 clBufferLength[bufferId]*clBufferStructureSize[bufferId],
													 ptr);	
				if (err == CL_SUCCESS) {
					return true;
				} else {
					*log << "\tOpenclManager::ReadBuffer: ERROR - unable to read buffer data" << std::endl;
					checkErr(err);
					return false;
				}
			}
			else {
				*log << "\tOpenclManager::ReadBuffer: ERROR - buffer is empty!" << std::endl;
			}
		} else {
			*log << "\tOpenclManager::ReadBuffer: ERROR - unable to find any buffer data" << std::endl;
		}
	}
	return false;
};

bool OpenclManager::RunKernels() {
	if (managerStatus) {
		cl_int err = queue.finish();
		if (err == CL_SUCCESS) {
			return true;
		} else {
			*log << "\tOpenclManager::RunKernels: ERROR - unable to run kernels" << std::endl;
			checkErr(err);
			return false;
		}
	} else {
		*log << "\tOpenclManager::RunKernels: ERROR -manager not initialized" << std::endl;
		return false;
	}
}

bool OpenclManager::ExecuteKernel(std::string kernelName) {
	if (managerStatus) {
		int kernelId = GetKernelId(kernelName);

		if (kernelId >= 0 && kernelId < kernels.size()) {
			if (kernelGlobalWorkRange[kernelId] > 0) {
				cl_int err = queue.enqueueNDRangeKernel(kernels[kernelId],
					cl::NullRange,
					cl::NDRange(kernelGlobalWorkRange[kernelId]),
					cl::NDRange(localWorkGroupSize),
					NULL,
					NULL);
				if (err == CL_SUCCESS) {
					return true;
				}
				else {
					*log << "\tOpenclManager::ExecuteKernel: ERROR - unable to execute kernel: " << kernelsNames[kernelId] << std::endl;
					checkErr(err);
					return false;
				}
			}
			else {
				*log << "\tOpenclManager::ExecuteKernel: ERROR - incorrect kernel work range" << std::endl;
				return false;
			}
		}
		else {
			*log << "\tOpenclManager::ExecuteKernel: ERROR - kernel structure is empty" << std::endl;
			return false;
		}
	}
	else {
		*log << "\tOpenclManager::ExecuteKernel: ERROR - manager not initialized" << std::endl;
		return false;
	}
};

bool OpenclManager::ExecuteKernel(int kernelId) {
	if (managerStatus) {
		if (kernelId >= 0 && kernelId < kernels.size()) {
			if (kernelGlobalWorkRange[kernelId] > 0) {
				cl_int err = queue.enqueueNDRangeKernel(kernels[kernelId],
														cl::NullRange,
														cl::NDRange(kernelGlobalWorkRange[kernelId]),
														cl::NDRange(localWorkGroupSize),
														NULL,
														NULL);
				if (err == CL_SUCCESS) {
					return true;
				} else {
					*log << "\tOpenclManager::ExecuteKernel: ERROR - unable to execute kernel: " << kernelsNames[kernelId] << std::endl;
					checkErr(err);
					return false;
				}
			} else {
				*log << "\tOpenclManager::ExecuteKernel: ERROR - incorrect kernel work range" << std::endl;
				return false;
			}
		} else {
			*log << "\tOpenclManager::ExecuteKernel: ERROR - kernel structure is empty" << std::endl;
			return false;
		}
	} else {
		*log << "\tOpenclManager::ExecuteKernel: ERROR - manager not initialized" << std::endl;
		return false;
	}
};

bool OpenclManager::GetStatus() {
	return managerStatus;
};

int	OpenclManager::GetNumberOfKernels() {
	if (!kernels.empty()) {
		return static_cast<int>(kernels.size());
	} else {
		return -1;
	}
}

void	OpenclManager::SetKernelGlobalWorkRange(std::string name, int workRange) {
	if (managerStatus) {
		int kernelId = -1;
		for (int ii = 0; ii < kernelsNames.size(); ii++){
			if (name == kernelsNames[ii]) kernelId = ii;
		}

		if (!kernels.empty()) {
			if (kernelId != -1 && kernelId < kernels.size()) {
				kernelGlobalWorkRange[kernelId] = workRange;
			} else {
				*log << "\tOpenclManager::SetKernelGlobalWorkRange: ERROR - wrong kernel ID invoked!!" << std::endl;
			}
		} else {
			*log << "\tOpenclManager::SetKernelGlobalWorkRange: ERROR - there are no kernels present in this system!!" << std::endl;
		}
	}
};

void OpenclManager::Close() {
	*log << "\tOpenClManager:" << "\t" << "manager cleanup...." << std::endl;
};

int	OpenclManager::GetKernelId(std::string name) {
	if (managerStatus) {
		for (int ii = 0; ii < kernelsNames.size(); ii++) {
			if (name == kernelsNames[ii]) {
				return ii;
			}
		}
		return -1;
	} else {
		return -1;
	}
}

void OpenclManager::checkErr(cl_int err) {
	switch (err) {
		case CL_INVALID_COMMAND_QUEUE:
			*log << "\t OpenclManager::checkErr: command_queue is not a valid command-queue" << std::endl;
			break;
			
		case CL_INVALID_CONTEXT:
			*log << "\t OpenclManager::checkErr: context associated with command_queue and buffer are not the	same or if the context associated with command_queue and events in	event_wait_list are not the same" << std::endl;
			break;

		case CL_INVALID_MEM_OBJECT:
			*log << "\t OpenclManager::checkErr: buffer is not a valid buffer object" << std::endl;
			break;

		case CL_INVALID_VALUE:
			*log << "\t OpenclManager::checkErr: the region being read specified by (offset, cb) is out of bounds or if ptr is a NULL value" << std::endl;
			break;
		
		case CL_INVALID_EVENT_WAIT_LIST:
			*log << "\t OpenclManager::checkErr: event_wait_list is NULL and num_events_in_wait_list greater than 0, or event_wait_list is not NULL and num_events_in_wait_list is 0, or if event objects in event_wait_list are not valid events" << std::endl;
			break;

		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			*log << "\t OpenclManager::checkErr: there is a failure to allocate memory for data store associated with buffer" << std::endl;
			break;

		case CL_OUT_OF_HOST_MEMORY:
			*log << "\t OpenclManager::checkErr: there is a failure to allocate resources required by the OpenCL implementation on the host" << std::endl;
			break;

		case CL_OUT_OF_RESOURCES:
			*log << "\t OpenclManager::checkErr: out of resources" << std::endl;
			break;

		case CL_INVALID_PROGRAM:
			*log << "\t OpenclManager::checkErr: program is not a valid program object!" << std::endl;
			break;
	
		case CL_INVALID_PROGRAM_EXECUTABLE:
			*log << "\t OpenclManager::checkErr: there is no successfully built executable for program!" << std::endl;
			break;
	
		case CL_INVALID_KERNEL_NAME:
			*log << "\t OpenclManager::checkErr: kernel_name is not found in program!" << std::endl;
			break;
				
		case CL_INVALID_KERNEL:
			*log << "\t OpenclManager::checkErr: kernel is not a valid kernel object" << std::endl;
			break;

		case CL_INVALID_KERNEL_ARGS:
			*log << "\t OpenclManager::checkErr: the kernel argument values have not been specified" << std::endl;
			break;
		
		case CL_INVALID_WORK_DIMENSION:
			*log << "\t OpenclManager::checkErr: work_dim is not a valid value(i.e.a value between 1 and 3)" << std::endl;
			break;
		
		case CL_INVALID_WORK_GROUP_SIZE:
			*log << "\t OpenclManager::checkErr: local_work_size is specified and number of work - items specified by global_work_size is not evenly divisable by size of work - group given by local_work_size or does not match the work - group size specified for kernel using the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier in program source" << std::endl;
			break;

		case CL_INVALID_WORK_ITEM_SIZE:
			*log << "\t OpenclManager::checkErr: the number of work - items specified in any of local_work_size[0], ... local_work_size[work_dim - 1] is greater than the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0], ....CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim - 1]" << std::endl;
			break;

		case CL_INVALID_GLOBAL_OFFSET:
			*log << "\t OpenclManager::checkErr: global_work_offset is not NULL" << std::endl;
			break;

		default:
			*log << "\t OpenclManager::checkErr: unknown error code" << std::endl;
			break;
	}
}