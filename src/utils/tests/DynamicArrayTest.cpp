/*
 * DynamicArrayTest.cpp
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#include "DynamicArrayTest.h"

#include "utils/DynamicArray.h"
#include "utils/tests/DummyObject.h"


TEST_SUITE_REGISTRATION(utils::DynamicArrayTest);

utils::DynamicArrayTest::DynamicArrayTest() { }

utils::DynamicArrayTest::~DynamicArrayTest() { }

void utils::DynamicArrayTest::testZeroInitialization() {
	DynamicArray<DummyObject, true, true> arr(0);
	ASSERT_TRUE(arr.empty());
	ASSERT_EQUAL(size_t (0), arr.size());
	ASSERT_EQUAL(size_t (0), arr.capacity());

	ASSERT_TRUE( ! (arr.begin() != arr.end())); // assert both are the same
}


void utils::DynamicArrayTest::testInitialization() {
	DynamicArray<DummyObject, true, true> arr(4);

	ASSERT_TRUE(arr.empty());
	ASSERT_EQUAL(size_t (0), arr.size());
	ASSERT_EQUAL(size_t (4), arr.capacity());

	ASSERT_TRUE( ! (arr.begin() != arr.end())); // assert both are the same

	for (int i = 0; i < 3; i++) {
		DummyObject obj(i);
		arr.push_back(obj);
	}

	ASSERT_TRUE(!arr.empty());
	ASSERT_EQUAL(size_t (3), arr.size());
	ASSERT_EQUAL(size_t (4), arr.capacity());

	DynamicArray<DummyObject, true, true>::iterator it = arr.begin();
	int i = 0;
	while (it != arr.end()) {
		ASSERT_EQUAL((*it).getData(), arr[i].getData());
		ASSERT_EQUAL(it->getData(), arr[i].getData());
		i++;
		++it;
	}
}


void utils::DynamicArrayTest::testReallocationForCopyConstruct() {
	DynamicArray<DummyObject, true, true> arr(4, 1.5);
	testReallocation(arr);
}


void utils::DynamicArrayTest::testReallocationForRealloc() {
	DynamicArray<DummyObject, false, true> arr(4, 1.5);
	testReallocation(arr);
}

template <class vector>
void utils::DynamicArrayTest::testReallocation(vector& v) {


	for (int i = 0; i < 5; i++) {
			DummyObject obj(i);
			v.push_back(obj);
	}

	ASSERT_TRUE(!v.empty());
	ASSERT_EQUAL(size_t (5), v.size());
	ASSERT_EQUAL(size_t (6), v.capacity());

	typename vector::iterator it = v.begin();

	while (it != v.end()) {
		if (it->getData() == 1 || it->getData() == 2) {
			v.erase(it);
		} else {
			++it;
		}
	}

	ASSERT_TRUE(!v.empty());
	ASSERT_EQUAL(size_t (3), v.size());
	ASSERT_EQUAL(size_t (4), v.capacity());

	// the last element has to be copied to the 3rd position;
	for (int i = 0; i < 3; i++) {
		if (i != 1 && i != 2) {
			ASSERT_EQUAL(v[i].getData(), i);
		} else if (i == 1){
			ASSERT_EQUAL(v[i].getData(), 4);
		} else if (i == 2) {
			ASSERT_EQUAL(v[i].getData(), 3);
		}
	}
}

