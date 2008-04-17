#include "Plugin.h"

int Plugin::init(int argc, char** argv, ostream* out)
{
   *out << "Warning: initializing " << _id << ", a generic plugin. This should never occur.\n";
   return 1;
}

int Plugin::run(int argc, char** argv, ostream* out)
{
   *out << "Warning: executing " << _id << ", a generic plugin. This should never occur.\n";
   return 1;
}

int Plugin::status(int argc, char** argv, ostream* out)
{
   *out << "Warning: monitoring " << _id << ", a generic plugin. This should never occur.\n";
   return 1;
}

int Plugin::signalhandler(Signal s, ostream* out)
{
   *out << "Warning: generic plugin " << _id << " received a signal. This should never occur.\n";
   return 1;
}
