/*
 * DummyObject.cpp
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#include "DummyObjectWithPtr.h"

int utils::DummyObjectWithPtr::_numInstances = 0;

int utils::DummyObjectWithPtr::_numConstructorCalls = 0;

int utils::DummyObjectWithPtr::_numDestructorCalls = 0;

utils::DummyObjectWithPtr::DummyObjectWithPtr(int data) : _id(_numInstances), _data(data){
	_numInstances++;
	_numConstructorCalls++;

	_ptr = new int[3];
	_ptr[0] = 0;
	_ptr[1] = 1;
	_ptr[2] = 2;
}

utils::DummyObjectWithPtr::DummyObjectWithPtr(const DummyObjectWithPtr& other) {
	_numInstances++;
	_id = other._id;
	_data = other._data;
	_numConstructorCalls++;

	_ptr = new int[3];
	_ptr[0] = other._ptr[0];
	_ptr[1] = other._ptr[1];
	_ptr[2] = other._ptr[2];
}


utils::DummyObjectWithPtr::~DummyObjectWithPtr() {
	_numDestructorCalls++;
	delete[] _ptr;
}

int utils::DummyObjectWithPtr::getID() {
	return _id;
}

void utils::DummyObjectWithPtr::setData(int data) {
	_data = data;
}

int utils::DummyObjectWithPtr::getData() {
	return _data;
}
