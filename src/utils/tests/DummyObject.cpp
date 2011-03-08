/*
 * DummyObject.cpp
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#include "DummyObject.h"

int utils::DummyObject::_numInstances = 0;

utils::DummyObject::DummyObject(int data) : _id(_numInstances), _data(data){
	_numInstances++;
}

utils::DummyObject::~DummyObject() {
}

int utils::DummyObject::getID() {
	return _id;
}

void utils::DummyObject::setData(int data) {
	_data = data;
}

int utils::DummyObject::getData() {
	return _data;
}
