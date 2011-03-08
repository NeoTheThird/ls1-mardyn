/*
 * DummyObject.h
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#ifndef DUMMYOBJECTWITHPTR_H_
#define DUMMYOBJECTWITHPTR_H_


namespace utils {
	class DummyObjectWithPtr;
}

/**
 * This is a dummy class to help test the DynamicArray, especially check the
 * performance vs. std::vector, if the DynamicArray just copies the objects instead
 * of copy-constructing. In it's constructor there's an array being allocated.
 */
class utils::DummyObjectWithPtr {

private:
	static int _numInstances;

	int _id;

	int _data;

	int* _ptr;


public:

	static int _numConstructorCalls;

	static int _numDestructorCalls;

	DummyObjectWithPtr(const DummyObjectWithPtr& other);

	DummyObjectWithPtr(int data);

	virtual ~DummyObjectWithPtr();

	int getID();

	void setData(int data);

	int getData();

};

#endif /* DUMMYOBJECTWITHPTR_H_ */
