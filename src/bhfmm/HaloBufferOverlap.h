/*
 * HaloBufferOverlap.h
 *
 *  Created on: Apr 26, 2016
 *      Author: obi
 */

#ifndef HALOBUFFEROVERLAP_H_
#define HALOBUFFEROVERLAP_H_
#include "bhfmm/utils/Vector3.h"


#ifdef ENABLE_MPI


namespace bhfmm{
template <class T>
class HaloBufferOverlap {
public:
	HaloBufferOverlap(Vector3<int> areaHaloSize, Vector3<int> edgeHaloSize,
		int cornerHaloSize, MPI_Comm comm, std::vector<int>& areaNeighbours,std::vector<int>& edgeNeighbours,std::vector<int>& cornerNeighbours, int isSend);
	virtual ~HaloBufferOverlap();
	void startCommunication();
	//communicate without persistent sends and receives
	void communicate();

	void wait();
	int testIfFinished();
	std::vector<T *>& getAreaBuffers(){
		return _areaBuffers;
	}
	std::vector<T *>& getEdgeBuffers(){
		return _edgeBuffers;
	}
	std::vector<T *>& getCornerBuffers(){
		return _cornerBuffers;
	}

	int  getAreaSize(){
		return _areaHaloSize;
	}
	int  getEdgeSize(){
		return _edgeHaloSize;
	}
	int  getCornerSize(){
		return _cornerHaloSize;
	}

	void clear();
private:
	void initCommunicationDouble();
	void fillArraySizes(Vector3<int> areaSizes, Vector3<int> edgeSizes);
	std::vector<T *>  _areaBuffers,  _edgeBuffers,  _cornerBuffers; //arrays for MPI halo transfer (send)
	Vector3<int> _areaHaloSize,_edgeHaloSize;
	int _cornerHaloSize;
	//these arrays save for every halo element the specific size -> if neighbour rank order is changed these arrays have to be changed too but nothing else
	int* _areaHaloSizes, *_edgeHaloSizes; //corner Arrays are always of the same size!
	std::vector<int> _areaNeighbours, _edgeNeighbours,_cornerNeighbours;
	MPI_Request * _areaRequests, *_edgeRequests, *_cornerRequests;
	MPI_Comm _comm;
	int _isSend;
};

#include <algorithm>


template <class T>
HaloBufferOverlap<T>::HaloBufferOverlap(Vector3<int> areaHaloSize, Vector3<int> edgeHaloSize,
		int cornerHaloSize, MPI_Comm comm, std::vector<int>& areaNeighbours,std::vector<int>& edgeNeighbours,std::vector<int>& cornerNeighbours, int isSend):
_areaBuffers(6), _edgeBuffers(12), _cornerBuffers(8), _areaNeighbours(areaNeighbours), _edgeNeighbours(edgeNeighbours), _cornerNeighbours(cornerNeighbours) {
//	for(int i=0; i<areaNeighbours.size();i++){
//		_areaNeighbours[i] = areaNeighbours[i];
//	}
//	for(int i=0; i<edgeNeighbours.size();i++){
//		_edgeNeighbours[i] = edgeNeighbours[i];
//	}
//	for(int i=0; i<cornerNeighbours.size();i++){
//		_cornerNeighbours[i] = cornerNeighbours[i];
//	}

	_areaHaloSize = areaHaloSize;
	_edgeHaloSize = edgeHaloSize;
	_cornerHaloSize = cornerHaloSize;
	_areaRequests = new MPI_Request[_areaBuffers.size()];
	_edgeRequests = new MPI_Request[_edgeBuffers.size()];
	_cornerRequests = new MPI_Request[_cornerBuffers.size()];

	fillArraySizes(areaHaloSize,edgeHaloSize);

//	_cornerHaloSizes = new int[_cornerBuffers.size()];
	_comm = comm;
	for(unsigned int i=0; i<_areaBuffers.size();i++){
		_areaBuffers[i] = new T[_areaHaloSizes[i]];
	}
	for(unsigned int i=0; i<_edgeBuffers.size();i++){
		_edgeBuffers[i] = new T[_edgeHaloSizes[i]];
	}
	for(unsigned int i=0; i<_cornerBuffers.size();i++){
		_cornerBuffers[i] = new T[cornerHaloSize];
	}
	_isSend = isSend;
	//initCommunicationDouble();

}

template <class T>
void HaloBufferOverlap<T>::fillArraySizes(Vector3<int> areaSizes, Vector3<int> edgeSizes){
	_areaHaloSizes = new int[_areaBuffers.size()];
	_edgeHaloSizes = new int[_edgeBuffers.size()];
	for(unsigned int i = 0; i < _areaBuffers.size(); i++){
		_areaHaloSizes[i] = areaSizes[i/2];
	}
	for(unsigned int i = 0; i < _edgeBuffers.size(); i++){
		_edgeHaloSizes[i] = edgeSizes[2-i/4];
	}
}
template <class T>
HaloBufferOverlap<T>::~HaloBufferOverlap() {
	for(unsigned int i=0; i<_areaBuffers.size();i++){
//		MPI_Request_free(&_areaRequests[i]);
		delete[] _areaBuffers[i];
	}
	for(unsigned int i=0; i<_edgeBuffers.size();i++){
//		MPI_Request_free(&_edgeRequests[i]);
		delete[] _edgeBuffers[i];
	}
	for(unsigned int i=0; i<_cornerBuffers.size();i++){
//		MPI_Request_free(&_cornerRequests[i]);
		delete[] _cornerBuffers[i];
	}
}

template <class T>
void HaloBufferOverlap<T>::clear(){
	for(unsigned int i = 0; i < _areaBuffers.size(); i++){
		std::fill(_areaBuffers[i], _areaBuffers[i] + _areaHaloSizes[i] , 0.0);
	}
	for(unsigned int i = 0; i < _edgeBuffers.size(); i++){
		std::fill(_edgeBuffers[i], _edgeBuffers[i] + _edgeHaloSizes[i] , 0.0);
	}
	for(unsigned int i = 0; i < _cornerBuffers.size(); i++){
		std::fill(_cornerBuffers[i], _cornerBuffers[i] + _cornerHaloSize , 0.0);
	}
}

//we assume here that the neighbour arrays are sorted in the way that sites and their opposite sites are always alternating in the array
template <class T>
void HaloBufferOverlap<T>::initCommunicationDouble(){
	for (unsigned int i = 0; i < _areaBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);
			//MPI_Rsend_init(_areaBuffers[i], _areaHaloSize, MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);

		}
		else{
			//adjusts that the tag of receive corresponds to send
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42 + indexShift, _comm, &_areaRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _edgeBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);
			//MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSize, MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);

		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42 + indexShift, _comm, &_edgeRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _cornerBuffers.size(); i++){
		if(_isSend){
			MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		//	MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Recv_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42 + indexShift, _comm, &_cornerRequests[i]);
		}
	}
}

template <class T>
void HaloBufferOverlap<T>::communicate(){
	for (unsigned int i = 0; i < _areaBuffers.size(); i++){
		if(_isSend){
			MPI_Irsend(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);
			//MPI_Rsend_init(_areaBuffers[i], _areaHaloSize, MPI_DOUBLE, _areaNeighbours[i], i + 42, _comm, &_areaRequests[i]);

		}
		else{
			//adjusts that the tag of receive corresponds to send
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Irecv(_areaBuffers[i], _areaHaloSizes[i], MPI_DOUBLE, _areaNeighbours[i], i + 42 + indexShift, _comm, &_areaRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _edgeBuffers.size(); i++){
		if(_isSend){
			MPI_Irsend(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);
			//MPI_Rsend_init(_edgeBuffers[i], _edgeHaloSize, MPI_DOUBLE, _edgeNeighbours[i], i + 42, _comm, &_edgeRequests[i]);

		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Irecv(_edgeBuffers[i], _edgeHaloSizes[i], MPI_DOUBLE, _edgeNeighbours[i], i + 42 + indexShift, _comm, &_edgeRequests[i]);
		}
	}
	for (unsigned int i = 0; i < _cornerBuffers.size(); i++){
		if(_isSend){
			MPI_Irsend(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		//	MPI_Rsend_init(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42, _comm, &_cornerRequests[i]);
		}
		else{
			int indexShift = (i%2 == 0)? +1: -1;
			MPI_Irecv(_cornerBuffers[i], _cornerHaloSize, MPI_DOUBLE, _cornerNeighbours[i], i + 42 + indexShift, _comm, &_cornerRequests[i]);
		}
	}
}

template <class T>
void HaloBufferOverlap<T>::startCommunication(){
	 MPI_Startall(_areaBuffers.size(), _areaRequests);
	 MPI_Startall(_edgeBuffers.size(), _edgeRequests);
	 MPI_Startall(_cornerBuffers.size(), _cornerRequests);
//	 std::cout << _areaBuffers.size() << _edgeBuffers.size() << _cornerBuffers.size() <<"\n";
}

template <class T>
void HaloBufferOverlap<T>::wait(){

	MPI_Status * areaStatusArray = new MPI_Status[_areaBuffers.size()];
	MPI_Waitall(_areaBuffers.size(),_areaRequests, areaStatusArray);

	MPI_Status * edgeStatusArray = new MPI_Status[_edgeBuffers.size()];
	MPI_Waitall(_edgeBuffers.size(),_edgeRequests, edgeStatusArray);

	MPI_Status * cornerStatusArray = new MPI_Status[_cornerBuffers.size()];
	MPI_Waitall(_cornerBuffers.size(),_cornerRequests, cornerStatusArray);

}

template <class T>
int HaloBufferOverlap<T>::testIfFinished(){
	int areaFlag, edgeFlag, cornerFlag;
	MPI_Status * areaStatusArray = new MPI_Status[_areaBuffers.size()];
	MPI_Testall(_areaBuffers.size(),_areaRequests, &areaFlag, areaStatusArray);

	MPI_Status * edgeStatusArray = new MPI_Status[_edgeBuffers.size()];
	MPI_Testall(_edgeBuffers.size(),_edgeRequests, &edgeFlag, edgeStatusArray);

	MPI_Status * cornerStatusArray = new MPI_Status[_cornerBuffers.size()];
	MPI_Testall(_cornerBuffers.size(),_cornerRequests, &cornerFlag, cornerStatusArray);
	return areaFlag * edgeFlag * cornerFlag;
}
}
#endif
#endif /* HALOBUFFEROVERLAP_H_ */
