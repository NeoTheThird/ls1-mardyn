/*
 * xyVal.h
 *
 *  Created on: 28.08.2013
 *      Author: mheinen
 */

#ifndef XYVAL_H_
#define XYVAL_H_

class xyVal
{
public:
	xyVal(double dXval, double dYval);
	~xyVal();

	double GetXvalue() {return _dXval;}
	double GetYvalue() {return _dYval;}

private:
	double _dXval;
	double _dYval;

};  // class xyVal


#endif /* XYVAL_H_ */
