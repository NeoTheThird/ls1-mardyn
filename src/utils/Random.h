#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>

class Random
{
 public:
   Random();
   void init(int seed);

   float rnd();

   int getIX() { return this->ix; }
   
   /** returns a gaussian distributed deviate with zero mean and a standard deviation of stdDeviation
    the returned value is in the range +/- infinity (better: smallest, largest double number)*/
   double gaussDeviate(double stdDeviation);

 private:
   int ix, iy;
   float am;
};

#endif

