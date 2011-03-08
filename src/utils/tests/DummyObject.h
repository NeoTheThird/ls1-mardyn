/*
 * DummyObject.h
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#ifndef DUMMYOBJECT_H_
#define DUMMYOBJECT_H_

namespace utils {
	class DummyObject;
}

/**
 * This is a dummy class to help test the DynamicArray.
 */
class utils::DummyObject {

private:
	static int _numInstances;

	int _id;

	int _data;

public:

	DummyObject(int data);

	virtual ~DummyObject();

	int getID();

	void setData(int data);

	int getData();

};

#endif /* DUMMYOBJECT_H_ */
