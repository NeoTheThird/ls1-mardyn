/*
 * DynamicArrayPerformanceTest.h
 *
 * @Date: 25.02.2011
 * @Author: eckhardw
 */

#ifndef DYNAMICARRAYPERFORMANCETEST_H_
#define DYNAMICARRAYPERFORMANCETEST_H_

#include "utils/Testing.h"

namespace utils {
	class DynamicArrayPerformanceTest;
}

/**
 * The performance tests consists of adding elements to the vector
 * via push_back() and removing them via pop_back().
 *
 * The runtime has to be a certain fraction of that of std::vector.
 */
class utils::DynamicArrayPerformanceTest : public utils::Test {

	TEST_SUITE(DynamicArrayPerformanceTest);
	TEST_METHOD(testPushPopSimpleObjectRealloc);
	TEST_METHOD(testPushPopSimpleObjectReallocShrink);
	TEST_METHOD(testPushPopObjectWithPtrRealloc);
	TEST_METHOD(testPushPopObjectWithPtrReallocShrink);

	TEST_METHOD(testPushSimpleObjectRealloc);
	TEST_METHOD(testPushObjectWithPtrRealloc);
	TEST_SUITE_END();

public:

	DynamicArrayPerformanceTest();

	virtual ~DynamicArrayPerformanceTest();

	//! test push / pop with simple object with realloc, but without shrink
	void testPushPopSimpleObjectRealloc();

	//! test push / pop with simple object with realloc, with shrink
	void testPushPopSimpleObjectReallocShrink();

	//! test push / pop with object with pointer data with realloc, but without shrink
	void testPushPopObjectWithPtrRealloc();

	//! test push / pop with object with pointer data with realloc, with shrink
	void testPushPopObjectWithPtrReallocShrink();


	//! test push with simple object with realloc
	void testPushSimpleObjectRealloc();

	//! test push with object with pointer data with realloc
	void testPushObjectWithPtrRealloc();


	/**
	 * Test the performance for a huge number of objects pushed back and inserted.
	 *
	 * @param v the vector type which is tested against std::vector
	 * @param factor the factor by which the runtime should at least differ.
	 */
	template <class object, bool copyconstruct, bool shrink>
	void testPushPop(double factor, double capacityIncrement);

	/**
	 * Test the performance for pushing back a huge number of objects.
	 *
	 * @param v the vector type which is tested against std::vector
	 * @param factor the factor by which the runtime should at least differ.
	 */
	template <class object, bool copyconstruct, bool shrink>
	void testPush(double factor, double capacityIncrement);
};

#endif /* DYNAMICARRAYPERFORMANCETEST_H_ */
