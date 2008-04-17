#include "Signal.h"

#include <iostream>
#include <string>
#include <queue>

using namespace std;

//! @brief abstract class / interface for plugins
//!
//! This is an abstract class which provides only virtual methods.
//! As we discussed in depth, virtual methods rely on an additional
//! indirection, but this does not matter at all for plugins because
//! they do not significantly affect the overall performance of ls1r1.
//!
//! Actual plugins extend this class by implementing its methods.
//! An exception is the signal queue, which is already implemented by
//! the generic Plugin class.
//!
class Plugin
{
 public:
   Plugin(unsigned long n) { this->_id = n; this->_name = "generic"; }
   virtual ~Plugin() {};

   //! @brief (re)sets the plugin to its initial state
   //!
   virtual int init(int argc, char** argv, ostream* out);
   //
   //! @brief (re)sets the plugin to its initial state
   //!
   virtual int init() { this->init(0, NULL, &cerr); return 0; }

   //! @brief the actual "main" function of the plugin
   //!
   virtual int run(int argc, char** argv, ostream* out);
   //
   //! @brief the actual "main" function of the plugin
   //!
   virtual int run() { this->run(0, NULL, &cerr); return 0; }

   //! @brief output information on the current state of the plugin
   //!
   virtual int status(int argc, char** argv, ostream* out);
   //
   //! @brief output information on the current state of the plugin
   //!
   virtual int status() { this->run(0, NULL, &cerr); return 0; }

   //! @brief are there any outgoing signals?
   //!
   bool hasSignal() { return !this->_outqueue.empty(); }
   //
   //! @brief should be called periodically by the plugin monitor for all active plugins
   //!
   Signal getSignal() { Signal s = _outqueue.front(); _outqueue.pop(); return s; }
   //
   //! @brief process an incoming signal
   //!
   virtual int signalhandler(Signal s, ostream* out);

   unsigned long getID() { return this->_id; }
   string getName() { return this->_name; }

 private:
   //! plugin identifier
   unsigned long _id;
   //! type/action descriptor
   string _name;
   //! stores all outgoing signals
   queue<Signal> _outqueue;

   //! @brief appends a signal to the outqueue, but particular plugins could e.g. implement priority queues
   virtual void signalize(Signal s) { this->_outqueue.push(s); }
};
