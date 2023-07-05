#pragma once
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <conio.h>
#include "CL/cl.hpp"
#include "BufferIndex.h"
#include "log_stream.h"
#include <memory>

/*
 * Wrapper clss allowing to manage the memory,
 * configure the device and perform calculations suing OpenCl library.
 */
class OpenclManager
{
	bool								managerStatus;
	int									localWorkGroupSize;
	int									global_mem_avalilable;
	int									global_mem_allocated;
	std::shared_ptr<log_stream>			log;


	std::vector<std::string>			kernelsNames;
	std::vector<cl::Kernel>				kernels;
	std::vector<int>					kernelGlobalWorkRange;

	std::vector<cl::Buffer>				clBuffer;
	std::vector<int>					clBufferSize; 
	std::vector<int>					clBufferElements;
	std::vector<int>					clBufferLength;
	std::vector<int>					clBufferStructureSize;

	cl::Device							default_device;
	cl::Platform						default_platform;
	cl::Context							context;
	cl::CommandQueue					queue;

	void								checkErr(cl_int err);

public:

	bool								GetStatus();
	bool								RegisterKernel(const std::vector<std::string> &kernelNames, const std::string &kernelCode);
	bool								ExecuteKernel(int kernelId);
	bool								ExecuteKernel(std::string kernelName);
	bool								RunKernels();
	bool								Initialize(std::shared_ptr<log_stream> p);
	cl::Buffer*							GetBuffer(const int&);
	bool								ReadBuffer(int, void*);
	int									GetBufferSize(const int&);
	int									GetBufferElements(const int&);
	int									GetBufferLength(const int&);
	int									GetBufferStructureSize(const int&);
	int									GetLocalWorkGroupSize();
	int									GetNumberOfKernels();

	template<class systemType> bool		SetKernelArg(std::string name, cl_uint argId, const systemType &argPtr);
	template<class type>       int	    CreateBuffer(std::vector<type>& data, cl_mem_flags flag);
	template<class type>       int		CreateBuffer(int size, cl_mem_flags flag);
	template<class type>	   int	    WriteBuffer(int bufferId, std::vector<type>& data);
	void								SetKernelGlobalWorkRange(std::string name, int workRange);
	int									GetKernelId(std::string name);
	void								Close();
	
};


template<class systemType> bool	OpenclManager::SetKernelArg(std::string name, cl_uint argId, const systemType &argPtr) {
	if (managerStatus) {
		int kernelId = -1;
		for (int ii = 0; ii < kernelsNames.size(); ii++) {
			if (name == kernelsNames[ii]) kernelId = ii;
		}

		if (kernelId < kernels.size() && kernelId >= 0)	{
			cl_int err = kernels[kernelId].setArg<systemType>(argId, argPtr);
			if (err != CL_SUCCESS) {
				checkErr(err);
				*log << "\tOpenclManager::SetKernelArg: ERROR - " << "unable to set kernel argument" << std::endl;
				*log << "\t" << "Kernel id: " << kernelsNames[kernelId] << "\t" << "Argument id: " << argId << std::endl;
				return false;
			} else {
				return true;
			}
		} else {
			*log << "\tOpenclManager::SetKernelArg: ERROR - " << "unable to set kernel argument - system not initialized" << std::endl;
			return false;
		}
	} else {
		return false;
	}
};

template<class type> int OpenclManager::WriteBuffer(int bufferId, std::vector<type>& data) {
	if (data.size() > 0) {
		unsigned structureSize = sizeof(type);
		int      bufferLength = static_cast<int>(data.size());
		if (clBufferLength[bufferId] != bufferLength ) return false;
		if (clBufferStructureSize[bufferId] != structureSize) return false;
		
		auto err = queue.enqueueWriteBuffer(clBuffer[bufferId], CL_TRUE, 0, structureSize * bufferLength, &data.front());
		if (err == CL_SUCCESS) {			
			return static_cast<int>(1);
		} else {
			*log << "\tOpenclManager::WriteBuffer: ERROR - unable to wrtie data into buffer" << std::endl;
			checkErr(err);
			return -1;
		}
	} else {
		*log << "\tOpenclManager::WriteBuffer: ERROR -  source buffer is empty" << std::endl;
		return -2;
	}
};

template<class type> int OpenclManager::CreateBuffer(std::vector<type>& data, cl_mem_flags flag) {
	if (data.size() > 0) {
		unsigned structureSize = sizeof(type);
		int      bufferElements = static_cast<int>(data.size());
		int      adjustedSize = static_cast<int>(localWorkGroupSize * std::ceil((float)bufferElements / localWorkGroupSize));
		data.resize(adjustedSize);

		cl::Buffer buffer(context, flag, structureSize *  adjustedSize);
		auto err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, structureSize * adjustedSize, &data.front());
		if (err == CL_SUCCESS) {
			clBuffer.push_back(buffer);
			clBufferSize.push_back(bufferElements*structureSize);
			clBufferElements.push_back(bufferElements);
			clBufferLength.push_back(adjustedSize);
			clBufferStructureSize.push_back(structureSize);
			return static_cast<int>(clBuffer.size() - 1);
		} else {
			*log << "\tOpenclManager::CreateBuffer: ERROR - unable to create buffer" << std::endl;
			checkErr(err);
			return -1;
		}
	} else {
		*log << "\tOpenclManager::CreateBuffer: ERROR -  buffer is empty" << std::endl;
		return -2;
	}
};

template<class type> int OpenclManager::CreateBuffer(int size, cl_mem_flags flag) {
	if (size > 0) {
		unsigned			structureSize  = sizeof(type);
		int					bufferElements = size;
		int					adjustedSize   = static_cast<int>(localWorkGroupSize * std::ceil((float)bufferElements / localWorkGroupSize));
		std::vector<type>	empty_buffer;
		empty_buffer.resize(adjustedSize);

		cl::Buffer buffer(context, flag, structureSize *  adjustedSize);
		auto err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, structureSize * adjustedSize, &empty_buffer.front());
		if (err == CL_SUCCESS) {
			clBuffer.push_back(buffer);
			clBufferSize.push_back(bufferElements*structureSize);
			clBufferElements.push_back(bufferElements);
			clBufferLength.push_back(adjustedSize);
			clBufferStructureSize.push_back(structureSize);
			return static_cast<int>(clBuffer.size() - 1);
		} else {
			*log << "\tOpenclManager::CreateBuffer: ERROR - unable to create buffer" << std::endl;
			checkErr(err);
			return -1;
		}
	} else {
		*log << "\tOpenclManager::CreateBuffer: ERROR -  buffer is empty" << std::endl;
		return -1;
	}
};