/*
 * DynamicArrayTest.h
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#ifndef DYNAMICARRAYTEST_H_
#define DYNAMICARRAYTEST_H_

#include "utils/Testing.h"

namespace utils {
	class DynamicArrayTest;
}


class utils::DynamicArrayTest : public utils::Test {

	TEST_SUITE(DynamicArrayTest);
	TEST_METHOD(testZeroInitialization);
	TEST_METHOD(testInitialization);
	TEST_METHOD(testReallocationForCopyConstruct);
	TEST_METHOD(testReallocationForRealloc);
	TEST_SUITE_END();

public:

	DynamicArrayTest();

	virtual ~DynamicArrayTest();

	void testZeroInitialization();

	void testInitialization();

	void testReallocationForCopyConstruct();

	void testReallocationForRealloc();

	template <typename vector>
	void testReallocation(vector& v);

};

#endif /* DYNAMICARRAYTEST_H_ */
