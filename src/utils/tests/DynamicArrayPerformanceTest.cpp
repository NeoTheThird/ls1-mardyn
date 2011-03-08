/*
 * DynamicArrayPerformanceTest.cpp
 *
 * @Date: 25.02.2011
 * @Author: eckhardw
 */

#include "DynamicArrayPerformanceTest.h"
#include "utils/tests/DummyObject.h"
#include "utils/tests/DummyObjectWithPtr.h"
#include "utils/DynamicArray.h"
#include "utils/Timer.h"

#include <sstream>
#include <vector>

TEST_SUITE_REGISTRATION(utils::DynamicArrayPerformanceTest);


utils::DynamicArrayPerformanceTest::DynamicArrayPerformanceTest() { }

utils::DynamicArrayPerformanceTest::~DynamicArrayPerformanceTest() { }


void utils::DynamicArrayPerformanceTest::testPushPopSimpleObjectRealloc() {
	// pushing the objects the first time is twice as fast, the other operations should
	// be approximately equal to std::vector, so I'd expect a factor of 0.8
	testPushPop<DummyObject, false, false>(0.8, 2.0);
}


void utils::DynamicArrayPerformanceTest::testPushPopSimpleObjectReallocShrink() {
	// pushing is twice as fast, however shrinking may take a factor of 2
	// -> we try not to lose performance despite shrinking.
	testPushPop<DummyObject, false, true>(1.0, 2.0);
}


void utils::DynamicArrayPerformanceTest::testPushPopObjectWithPtrRealloc() {
	// pushing the objects the first time is twice as fast, the other operations should
	// be approximately equal to std::vector -> try factor 0.7 for repeated push/pop
	testPushPop<DummyObjectWithPtr, false, false>(0.8, 2.0);
}


void utils::DynamicArrayPerformanceTest::testPushPopObjectWithPtrReallocShrink() {
	// pushing is twice as fast, however shrinking may take a factor of 2
	// -> we try not to lose performance despite shrinking.
	testPushPop<DummyObjectWithPtr, false, true>(1.0, 2.0);
}


void utils::DynamicArrayPerformanceTest::testPushSimpleObjectRealloc() {
	// the runtime for the push operation should be approximatelly 0.5
	testPush<DummyObject, false, false>(0.6, 2.0);
}


void utils::DynamicArrayPerformanceTest::testPushObjectWithPtrRealloc() {
	// the runtime for the push operation should be approximatelly 0.5
	testPush<DummyObjectWithPtr, false, false>(0.6, 2.0);
}


template <class object, bool copyconstruct, bool shrink>
void utils::DynamicArrayPerformanceTest::testPush(double factor, double capacityIncrement) {

	int numObjects = 5000000;
	Timer vectorTimer;
	Timer stdvectorTimer;

	for (int i = 0; i < 3; i++) {
		DynamicArray<object, copyconstruct, shrink> v(0, capacityIncrement);
		std::vector<object> stdv;

		vectorTimer.start();
		for (int j = 0; j < numObjects; j++) {
			v.push_back(object(j));
		}
		vectorTimer.stop();

		stdvectorTimer.start();
		for (int j = 0; j < numObjects; j++) {
			stdv.push_back(object(j));
		}
		stdvectorTimer.stop();
	}

	std::cout << "TEST PUSH stdvectorTimer: " << stdvectorTimer.get_etime() << " custom vector " << vectorTimer.get_etime() << std::endl;
	std::cout << "TEST PUSH ratio is " << (vectorTimer.get_etime() / stdvectorTimer.get_etime()) << std::endl;

	std::stringstream msg;
	msg << "time for DynamicArray has to be at most " << factor << " * time for std::Vector! (time DynamicArray="
			<< vectorTimer.get_etime() << " std::vector=" << stdvectorTimer.get_etime() << ")" << std::endl;
	msg << "Note: this is a performance test, so it is not very meaningful if you run it on a system with heavy load." << std::endl;
	ASSERT_TRUE_MSG(msg.str(), stdvectorTimer.get_etime() * factor > vectorTimer.get_etime());
}


template <class object, bool copyconstruct, bool shrink>
void utils::DynamicArrayPerformanceTest::testPushPop(double factor, double capacityIncrement) {

	int numIterations = 2;
	int numObjects = 1000000;
	Timer vectorTimer;
	Timer stdvectorTimer;

	for (int i = 0; i < 3; i++) {
		DynamicArray<object, copyconstruct, shrink> v(0, capacityIncrement);
		std::vector<object> stdv;

		vectorTimer.start();
		for (int i = 0; i < numIterations; i++) {
			for (int j = 0; j < numObjects; j++) {
				v.push_back(object(j));
			}

			for (int j = numObjects-1; j >= 0; j--) {
				v.pop_back();
			}
		}
		vectorTimer.stop();

		stdvectorTimer.start();
		for (int i = 0; i < numIterations; i++) {
			for (int j = 0; j < numObjects; j++) {
				stdv.push_back(object(j));
			}

			for (int j = numObjects-1; j >= 0; j--) {
				stdv.pop_back();
			}
		}
		stdvectorTimer.stop();
	}

	std::cout << "TEST PUSH POP stdvectorTimer: " << stdvectorTimer.get_etime() << " custom vector " << vectorTimer.get_etime() << std::endl;
	std::cout << "TEST PUSH POP ratio is " << (vectorTimer.get_etime() / stdvectorTimer.get_etime()) << std::endl;

	std::stringstream msg;
	msg << "time for DynamicArray has to be at most " << factor << " * time for std::Vector! (time DynamicArray="
			<< vectorTimer.get_etime() << " std::vector=" << stdvectorTimer.get_etime() << ")" << std::endl;
	msg << "Note: this is a performance test, so it is not very meaningful if you run it on a system with heavy load." << std::endl;
	ASSERT_TRUE_MSG(msg.str(), stdvectorTimer.get_etime() * factor > vectorTimer.get_etime());
}
