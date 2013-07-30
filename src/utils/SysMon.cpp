/** \file SysMon.cpp
  * \brief System monitoring class
  * \author Martin Bernreuther <bernreuther@hlrs.de>
  * license: GPL (http://www.gnu.de/documents/gpl.de.html)
*/

#include "SysMon.h"
#include "utils/Logger.h"
#include <stack>

#include <fstream>
#include <sstream>


using namespace std;


#ifdef MPI_VERSION
	const MPI_Datatype SysMon::mpiTvalue=MPI_DOUBLE;
	const int SysMon::mpiRootRank=0;
#endif


void SysMon::clear()
{
	for (list<Expression>::iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
		exprit->clear();
	_expressions.clear();
	_values.clear();
	_initMinMax.clear();
	_valuesMaxMinPeak.clear();
#ifdef MPI_VERSION
	_valuesMaxMin.clear();
#endif
}

int SysMon::addExpression(const std::string& exprstr)
{
	_expressions.push_back(Expression(string(),_variableset));
	Expression& expr=_expressions.back();
	expr.initializeRPN(exprstr,true);	// generate a Label
	//expr.genLabel();	// done at the end of initializeRPN
	
	_values.push_back(0.);
	//_values.resize(_expressions.size(),0.);
	_initMinMax.push_back(true);
	//_initMinMax.resize(_expressions.size(),true);
	_valuesMaxMinPeak.push_back(0.); _valuesMaxMinPeak.push_back(0.);
	//_valuesMaxMinPeak.resize(_values.size()*2,0.);
#ifdef MPI_VERSION
	_valuesMaxMin.push_back(0.); _valuesMaxMin.push_back(0.);
	//_valuesMaxMin.resize(_valuesMaxMinPeak.size(),0.);
#endif
	return _expressions.size()-1;
}


void SysMon::updateExpressionValues(bool resetMinMax)
{
	if(resetMinMax) _initMinMax.assign(numExpressions(),true);
	
	//sync();
	if(_variableset->existVariableGroup("sysconf")) updateVariables_sysconf();
	if(_variableset->existVariableGroup("sysinfo")) updateVariables_sysinfo();
	if(_variableset->existVariableGroup("procmeminfo")) updateVariables_procmeminfo();
	if(_variableset->existVariableGroup("procvmstat")) updateVariables_procvmstat();
	if(_variableset->existVariableGroup("procloadavg")) updateVariables_procloadavg();
	if(_variableset->existVariableGroup("procselfstatm")) updateVariables_procselfstatm();
	if(_variableset->existVariableGroup("procselfschedstat")) updateVariables_procselfschedstat();
	if(_variableset->existVariableGroup("procselfsched")) updateVariables_procselfsched();
	
	size_t i=0;
	for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
	{
		_values[i]=exprit->evaluateFloat();
		++i;
	}
	
	vector<Tvalue> valuesMaxMin(2*_values.size());
	for(i=0;i<_values.size();++i)
	{
		valuesMaxMin[2*i]=_values[i];
		valuesMaxMin[2*i+1]=-_values[i];
	}
	
	for(i=0;i<_valuesMaxMinPeak.size();++i)
	{
		if(_initMinMax[i/2] || valuesMaxMin[i]>_valuesMaxMinPeak[i])
			_valuesMaxMinPeak[i]=valuesMaxMin[i];
	}
	
#ifdef MPI_VERSION
	int myrank;
	MPI_CHECK( MPI_Comm_rank(_mpicomm,&myrank) );
	MPI_CHECK( MPI_Reduce(valuesMaxMin.data(),_valuesMaxMin.data(),_valuesMaxMin.size(),mpiTvalue,MPI_MAX,0,MPI_COMM_WORLD) );
	
	if(myrank==mpiRootRank)
	{
		for(i=0;i<_valuesMaxMinPeak.size();++i)
		{
			if(_initMinMax[i/2] || _valuesMaxMin[i]>_valuesMaxMinPeak[i])
				_valuesMaxMinPeak[i]=_valuesMaxMin[i];
		}
	}
#endif
	_initMinMax.assign(_initMinMax.size(),false);
}

int SysMon::getExpressionIndex(const std::string& label) const
{
	int idx=0;
	for(list<Expression>::const_iterator it=_expressions.begin();it!=_expressions.end();++it)
	{
		if(it->getLabel()==label) return idx;
		++idx;
	}
	return -1;
}

#ifdef MPI_VERSION
pair<SysMon::Tvalue,SysMon::Tvalue> SysMon::getExpressionMinMaxValues(unsigned int index) const
{
	int myrank;
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD,&myrank) );
	if(myrank==mpiRootRank&&index*2+1<_valuesMaxMin.size())
		return make_pair(_valuesMaxMin[index*2],_valuesMaxMin[index*2+1]);
	else
		return make_pair(0,0);
}
#endif


void SysMon::writeExpressionValues(ostream& ostrm) const
{
	size_t numvalues=_values.size();
	size_t i=0;
#ifdef MPI_VERSION
	// MPI_Comm_get_attr(_mpicomm,MPI_HOST,&mpihost,&flag);
	int myrank;
	MPI_CHECK( MPI_Comm_rank(_mpicomm,&myrank) );
	if(myrank==mpiRootRank) {
		for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
		{
			ostrm << exprit->getLabel();
			if(i>=numvalues || _initMinMax[i])
			{
				ostrm << "\tundefined";
			} else {
				ostrm << "\t" << -_valuesMaxMin[2*i+1] << " .. " << _valuesMaxMin[2*i];
				ostrm << "\t" << -_valuesMaxMinPeak[2*i+1] << " .. " << _valuesMaxMinPeak[2*i];
				++i;
			}
			ostrm << endl;
		}
	}
#else
	for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
	{
		ostrm << exprit->getLabel();
		if(i>=numvalues || _initMinMax[i])
		{
			ostrm << "\tundefined";
		} else {
			ostrm << "\t" << _values[i];
			ostrm << "\t" << -_valuesMaxMinPeak[2*i+1] << " .. " << _valuesMaxMinPeak[2*i];
			++i;
		}
		ostrm << endl;
	}
#endif
}



unsigned int SysMon::updateVariables_sysconf()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_SYSCONF
	long val=0;
	
	val=sysconf(_SC_PHYS_PAGES);
	_variableset->setVariable("sysconf:PHYS_PAGES",val);
	++numvalues;
	val=sysconf(_SC_AVPHYS_PAGES);
	_variableset->setVariable("sysconf:AVPHYS_PAGES",val);
	++numvalues;
	val=sysconf(_SC_PAGESIZE);
	_variableset->setVariable("sysconf:PAGESIZE",val);
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_sysinfo()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_SYSINFO
	long val=0;
	
	struct sysinfo sysinfodata;
	int rc=sysinfo(&sysinfodata);
	
	if(!rc)
	{
		val=long(sysinfodata.uptime);
		_variableset->setVariable("sysinfo","uptime",val);
		++numvalues;
		val=long(sysinfodata.loads[0]);
		_variableset->setVariable("sysinfo","loads1",val);
		++numvalues;
		val=long(sysinfodata.loads[1]);
		_variableset->setVariable("sysinfo","loads5",val);
		++numvalues;
		val=long(sysinfodata.loads[2]);
		_variableset->setVariable("sysinfo","loads15",val);
		++numvalues;
		val=long(sysinfodata.totalram);
		_variableset->setVariable("sysinfo","totalram",val);
		++numvalues;
		val=long(sysinfodata.freeram);
		_variableset->setVariable("sysinfo","freeram",val);
		++numvalues;
		val=long(sysinfodata.sharedram);
		_variableset->setVariable("sysinfo","sharedram",val);
		++numvalues;
		val=long(sysinfodata.bufferram);
		_variableset->setVariable("sysinfo","bufferram",val);
		++numvalues;
		val=long(sysinfodata.totalswap);
		_variableset->setVariable("sysinfo","totalswap",val);
		++numvalues;
		val=long(sysinfodata.freeswap);
		_variableset->setVariable("sysinfo","freeswap",val);
		++numvalues;
		val=long(sysinfodata.procs);
		_variableset->setVariable("sysinfo","procs",val);
		++numvalues;
		val=long(sysinfodata.totalhigh);
		_variableset->setVariable("sysinfo","totalhigh",val);
		++numvalues;
		val=long(sysinfodata.freehigh);
		_variableset->setVariable("sysinfo","freehigh",val);
		++numvalues;
		val=long(sysinfodata.mem_unit);
		_variableset->setVariable("sysinfo","mem_unit",val);
		++numvalues;
	}
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procmeminfo()
{	// 
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCMEMINFO
	ifstream ifstrm("/proc/meminfo");
	if(!ifstrm) return 0;
	long val=0;
	string line,label;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{
		++linenr;
		istringstream iss(line);
		if (!(iss >> label >> val)) { break; }
		//iss.str(string());
		size_t i=label.find_first_of(" :");
		while(i!=string::npos)
		{
			label.erase(i,1);
			i=label.find_first_of(" :");
		}
		_variableset->setVariable("procmeminfo",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procvmstat()
{	//
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCVMSTAT
	ifstream ifstrm("/proc/vmstat");
	if(!ifstrm) return 0;
	long val=0;
	string line,label;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{
		++linenr;
		istringstream iss(line);
		if (!(iss >> label >> val)) { break; }
		//iss.str(string());
		_variableset->setVariable("procvmstat",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procloadavg()
{	// see also `man 5 proc | grep -m 1 -A 12 /proc/loadavg`
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCLOADAVG
	ifstream ifstrm("/proc/loadavg");
	if(!ifstrm) return 0;
	string line,label;
	if (!getline(ifstrm, line)) return 0;
	ifstrm.close();
	size_t i=line.find("/");
	if(i!=string::npos)
	{
		//line[i]=' ';
		line.replace(i,1," ");
	}
	float loadavg1,loadavg5,loadavg15;
	unsigned int numschedentexec,numschedentexist,pidrecent;
	istringstream iss(line);
	if (!(iss >> loadavg1 >> loadavg5 >> loadavg15
	          >> numschedentexec >> numschedentexist >> pidrecent)) { return 0; }
	//iss.str(string());
	_variableset->setVariable("procloadavg","loadavg1",double(loadavg1));
	++numvalues;
	_variableset->setVariable("procloadavg","loadavg5",double(loadavg5));
	++numvalues;
	_variableset->setVariable("procloadavg","loadavg15",double(loadavg15));
	++numvalues;
	_variableset->setVariable("procloadavg","numschedentexec",long(numschedentexec));
	++numvalues;
	_variableset->setVariable("procloadavg","numschedentexist",long(numschedentexist));
	++numvalues;
	//_variableset->setVariable("procloadavg","pidrecent",long(pidrecent)); ++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfstatm()
{	// see also `man 5 proc | grep -m 1 -A 12 /statm`
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSTATM
	ifstream ifstrm("/proc/self/statm");
	if(!ifstrm) return 0;
	unsigned long size,resident,share,text,lib,data;
	ifstrm >> size >> resident >> share >> text >> lib >> data;
	ifstrm.close();
	_variableset->setVariable("procselfstatm","size",long(size));
	++numvalues;
	_variableset->setVariable("procselfstatm","resident",long(resident));
	++numvalues;
	_variableset->setVariable("procselfstatm","share",long(share));
	++numvalues;
	_variableset->setVariable("procselfstatm","text",long(text));
	++numvalues;
	_variableset->setVariable("procselfstatm","lib",long(lib));
	++numvalues;
	_variableset->setVariable("procselfstatm","data",long(data));
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfsched()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSCHED
	ifstream ifstrm("/proc/self/sched");
	if(!ifstrm) return 0;
	double val=0;
	string line,label,colon;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{
		++linenr;
		istringstream iss(line);
		if (!(iss >> label >> colon >> val)) { continue; }
		//iss.str(string());
		if(colon!=":") continue;
		_variableset->setVariable("procselfsched",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfschedstat()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSCHEDSTAT
	ifstream ifstrm("/proc/self/schedstat");
	if(!ifstrm) return 0;
	unsigned long runningtime, waitingtime;
	unsigned int numtasks;
	if (!(ifstrm >> runningtime >> waitingtime >> numtasks)) { return 0; }
	ifstrm.close();
	_variableset->setVariable("procselfschedstat:runningtime",long(runningtime));
	++numvalues;
	_variableset->setVariable("procselfschedstat:waitingtime",long(waitingtime));
	++numvalues;
	_variableset->setVariable("procselfschedstat:numtasks",long(numtasks));
	++numvalues;
#endif
	return numvalues;
}

