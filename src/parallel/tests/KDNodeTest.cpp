/*
 * KDNodeTest.cpp
 *
 *  Created on: Feb 24, 2012
 *      Author: eckhardw
 */

#include "KDNodeTest.h"
#include "parallel/KDNode.h"

TEST_SUITE_REGISTRATION(KDNodeTest);

KDNodeTest::KDNodeTest() {
}

KDNodeTest::~KDNodeTest() {
}


void KDNodeTest::testEqual() {
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {3, 3, 3};
	bool coversAll[] = {true, true, true};

	KDNode node1(5, lowerEnd, upperEnd, 41, 0, coversAll);
	KDNode node2(5, lowerEnd, upperEnd, 41, 0, coversAll);
	ASSERT_TRUE(node1.equals(node2));

	KDNode node3(7, lowerEnd, upperEnd, 41, 0, coversAll);
	ASSERT_TRUE(! node1.equals(node3));

	KDNode* node4 = new KDNode(7, lowerEnd, upperEnd, 41, 0, coversAll);
	node1._child1 = node4;
	KDNode* node5 = new KDNode(7, lowerEnd, upperEnd, 41, 0, coversAll);
	node2._child1 = node5;
	// node4 and node5 are deleted by their root nodes...
	ASSERT_TRUE(node1.equals(node2));
}

void KDNodeTest::testBuildKDTree() {

	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 3, 3};
	bool coversAll[] = {true, true, true};

	KDNode node1(1, lowerEnd, upperEnd, 1, 0, coversAll);
	node1.buildKDTree();
	KDNode result(1, lowerEnd, upperEnd, 1, 0, coversAll);
	ASSERT_TRUE(node1.equals(result));


	KDNode root(2, lowerEnd, upperEnd, 1, 0, coversAll);
	root.buildKDTree();

	KDNode resultRoot(2, lowerEnd, upperEnd, 1, 0, coversAll);
	int lower1[] = {0, 0, 0};
	int upper1[] = {3, 3, 3};
	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild1 = new KDNode(1, lower1, upper1, 1, 0, childCoversAll);

	int lower2[] = {4, 0, 0};
	int upper2[] = {7, 3, 3};
	KDNode* resultChild2 = new KDNode(1, lower2, upper2, 1, 1, childCoversAll);
	resultRoot._child1 = resultChild1;
	resultRoot._child2 = resultChild2;

	ASSERT_TRUE(root.equals(resultRoot));
}

void KDNodeTest::testFindAreaForProcess() {
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 3, 3};
	bool coversAll[] = {true, true, true};
	KDNode resultRoot(3, lowerEnd, upperEnd, 1, 0, coversAll);

	int lower1[] = {0, 0, 0};
	int upper1[] = {3, 3, 3};
	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild1 = new KDNode(2, lower1, upper1, 1, 0, childCoversAll);
	int lower2[] = {4, 0, 0};
	int upper2[] = {7, 3, 3};
	KDNode* resultChild2 = new KDNode(1, lower2, upper2, 1, 2, childCoversAll);

	resultRoot._child1 = resultChild1;
	resultRoot._child2 = resultChild2;

	int lower3[] = {0, 0, 0};
	int upper3[] = {1, 3, 3};
//	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild3 = new KDNode(1, lower3, upper3, 1, 0, childCoversAll);

	int lower4[] = {2, 0, 0};
	int upper4[] = {3, 3, 3};
//	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild4 = new KDNode(1, lower4, upper4, 1, 1, childCoversAll);

	resultRoot._child1->_child1 = resultChild3;
	resultRoot._child1->_child2 = resultChild4;

	ASSERT_EQUAL(resultRoot.findAreaForProcess(0), resultChild3);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(1), resultChild4);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(2), resultChild2);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(41), (KDNode*) NULL);
}
