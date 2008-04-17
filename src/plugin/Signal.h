#include <string>
#include <iostream>

using namespace std;

class Signal
{
 public:
   Signal() { }
   virtual ~Signal() {}

   //! @brief output signal via a stream
   //!
   virtual void occur(ostream* out) { *out << "SIGNAL"; }
};
