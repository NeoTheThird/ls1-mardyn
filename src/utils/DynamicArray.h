/*
 * DynamicArray.h
 *
 * @Date: 22.02.2011
 * @Author: Wolfgang Eckhardt
 */

#ifndef DYNAMICARRAY_H_
#define DYNAMICARRAY_H_

#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <iostream>

namespace utils {
	template<typename T, bool copyconstruct, bool shrink>
	class DynamicArray;
}

// defined in MemoryManager.cpp
extern unsigned long dynamicArraySize;

/**
 * This is the implementation of a dynamic array (cfg. std::vector). In contrast
 * to std::vector this implementation may free memory (if shrink == true) if
 * elements are removed (i.e. it can shrink). Additionally, it may use realloc
 * (if copyconstruct==false) when the internal memory is relocated instead of
 * copy-constructing the elements.
 *
 * The public interface should be a subset of std::vector, and the methods implemented
 * here should comply with the specification of the corresponding methods of std::vector.
 * The rationale is that it should be possible to exchange both implementations.
 *
 * Elements stored in this dynamic array have to meet the following conditions:
 * - assignment operator fully implemented
 * - copy constructor fully implemented
 *
 * Performance:
 * - is best for a _capacityIncrement = 2.0, copyconstruct=false, shrink=false
 *   (runtime approx. 70 Perc. of std::vector)
 * - with a capacityIncrement = 2.0, copyconstruct=false, shrink=true
 *   (runtime approx. 80 Perc. of std::vector)
 * - with a capacityIncrement != 2.0, copyconstruct=false, shrink=true, the runtime
 *  equals approximately that of std::vector
 * As measure of performance I take the addition of 10e7 elements via push_back,
 * followed by their deletion via pop_back();
 *
 * @note the option copyconstruct=false will cause objects to be moved in memory,
 *       something which is not defined by the standard. Thus be sure this won't
 *       cause troubles for your objects!
 *
 * @todo Try to implement the resize not with allocate / assign / delete, but just
 *       using a realloc(). That should be fine if we use a placement-new().
 *
 * @todo Try to implement the vector so that it instantiates all objects which may
 *       be contained, then realize deletion just by decreasing the _finish-pointer
 *       and insertion by assignment?
 *
 */
template<typename T, bool copyconstruct, bool shrink>
class utils::DynamicArray {

public:

	/**
	 * Iterator class for the dynamic array. Use it to delete elements from
	 * the dynamic array. I strongly recommend to use it only in a while-loop:
	 *
	 * \verbatim
	 * iterator it = array.begin();
	 * while (it != array.end()) {
	 *   ....
	 *   if (delete) {
	 *   	it = array.erase(it);
	 *   } else {
	 *     ++it;
	 *   }
	 * }
	 * \endverbatim
	 *
	 */
	class iterator {

		friend class DynamicArray;

	private:

		T* _current;

		iterator(T* position) : _current(position) { }

	public:

		iterator() : _current(NULL) {}
		iterator(const DynamicArray<T,copyconstruct, shrink>::iterator& other) : _current(other._current) {}

		iterator& operator=(const iterator rhs) {
			this->_current = rhs._current;
			return *this;
		}

		//! support only prefix operator as it should be more efficient
		iterator& operator++() {
			_current++;
			return *this;
		}

		T& operator*() const {
			return *_current;
		}

		T* operator->() const {
			return _current;
		}

		bool operator!=(const iterator& rhs) const {
			return this->_current != rhs._current;
		}
	};


	DynamicArray(size_t capacity = 1, double capacityIncrement = 2.0)
		: _start(NULL), _finish(NULL), _endOfStorage(NULL), _capacityIncrement(capacityIncrement) {
		reallocate(capacity);
		//std::cout << "Created DynamicArray with size=" << capacity << std::endl;
	}

	//! copy construct a DynamicArray from another. The content is copied.
	DynamicArray(const DynamicArray& other)
		: _start(NULL), _finish(NULL), _endOfStorage(NULL), _capacityIncrement(other._capacityIncrement) {
		reallocate(other.capacity());
		insert(end(), other.begin(), other.end());
		//std::cout << "Created (by Copy) DynamicArray with size=" << other.capacity() << std::endl;
	}


	virtual ~DynamicArray() {

		dynamicArraySize -= capacity() * sizeof(T);

		// call the destructors for elements stored
		for (unsigned int i = 0; i < size(); i++) {
			_start[i].~T();
		}

		// free buffer
		if (copyconstruct) {
			operator delete(_start);
		} else {
			free(_start);
		}
	}


	//! @return the number of elements contained.
	size_t size() const {
		return size_t (_finish - _start);
	}


	//! @return the number of elements which could be stored.
	size_t capacity() const {
		return size_t (_endOfStorage - _start);
	}

	//! @return if the array is empty.
	bool empty() const {
		return _start == _finish;
	}


	//! @return a reference to the first element
	T& front() {
		return *_start;
	}


	//! @return a reference to the last element
	T& back() {
		return *(_finish-1);
	}

	void push_back(const T& data) {
		if (_finish == _endOfStorage) {
			// for constellations of small capacity() and capacityIncrement their product might equal the
			// old capacity, so we have to make sure the old capacity is at least increased by 1.
			size_t newCapacity = (capacity() * _capacityIncrement) > capacity()+1 ? (capacity() * _capacityIncrement) : capacity()+1;
			reallocate(newCapacity);
		}
		new (_finish) T(data);
		++_finish;
	}


	/**
	 * remove the last element from the array and shrink the array, if neccessary.
	 */
	void pop_back() {
		--_finish;
		_finish->~T();

		if (shrink && size() < capacity() / _capacityIncrement) {
			reallocate(capacity() / _capacityIncrement);
		}
	}


	/**
	 * Remove the element at position from this array and reallocate, if neccessary.
	 *
	 * @return iterator pointing to the next element
	 *
	 * @note To prevent defragmentation of memory, the last element is copy-constructed
	 * at the position of the element deleted.
	 */
	iterator& erase(iterator& position) {
		if (_finish > _start) {
			--_finish;
			*position._current = *_finish;
		} else {
			assert(_finish == _start);
			assert(_finish == position._current);
		}
		_finish->~T();

		if (shrink && size() < capacity() / _capacityIncrement) {
			size_t index = position._current - _start;
			reallocate(capacity() / _capacityIncrement);
			position._current = &_start[index];
		}

		return position;
	}

	/**
	 * @note This method should acutally be called append. For compatibility reasons
	 * 		 however I call it insert. Elements are always appended at the end.
	 *
	 * @param first / last: Iterators specifying a range of elements. Inserts the elements
	 *                      in the range [first,last). Notice that the range includes all
	 *                      the elements between first and last, including the element pointed
	 *                      by first but not the one pointed by last.
	 */
	void insert(const iterator& position, const iterator& first, const iterator& last) {
		if (copyconstruct) {
			T* elementToCopy = first._current;
			while (elementToCopy < last._current) {
				push_back(*elementToCopy);
				elementToCopy++;
			}
		} else {
			size_t sizeDelta = last._current - first._current;
			if (size() + sizeDelta > capacity()) {
				reallocate(size() + sizeDelta);
			}
			memcpy(_finish, first._current, sizeDelta * sizeof (T));
		}
	}

	//! @return an iterator pointing to the first element.
	iterator begin() const {
		return iterator(_start);
	}

	//! @return an iterator pointing behind the last element.
	iterator end() const {
		return iterator (_finish);
	}

	/**
	 * Remove all the elements, but don't free the memory. Consequently, use this
	 * method if you want to clear all data and refill the vector subsequently again.
	 */
	void clear() {
		// call the destructors for elements in old location
		for (size_t i = 0; i < size(); i++) {
			_start[i].~T();
		}

		_finish = _start;
	}


	//! access element at index
	T& operator[](size_t index) {
		return _start[index];
	}


	//! access element at index
	const T& operator[](size_t index) const {
		return _start[index];
	}

private:

	//! perform the reallocation if the capacity should be in- or decreased.
	void reallocate(size_t newSize) {

		dynamicArraySize += newSize * sizeof (T);
		dynamicArraySize -= capacity() * sizeof(T);

		T* _newStart;
		size_t numElements = size();

		if (copyconstruct) {
			_newStart = static_cast<T*> ( operator new (newSize * sizeof (T)));
			// copy-construct elements in new location
			for (unsigned int i = 0; i < numElements; i++) {
				new (&_newStart[i]) T(_start[i]);
			}

			// call the destructors for elements in old location
			for (unsigned int i = 0; i < numElements; i++) {
				_start[i].~T();
			}

			// free old buffer
			operator delete(_start);
		} else {
			_newStart = static_cast<T*> (realloc(_start, newSize * sizeof (T)));
		}

		// we need _start to determine size(), so assign _finish before changing _start
		_finish = &_newStart[numElements];
		_start = _newStart;
		_endOfStorage = &_start[newSize];

//		std::cout << "XXXXXXXXXXXXX" << std::endl;
//		std::cout << "DynamicArray reallocated to newsize=" << newSize << std::endl;
//		std::cout << "size="<< size() << " Capacity="<< capacity() << " _start=" << _start << " finish=" << _finish << " endOfStorage=" << _endOfStorage << std::endl;
//		std::cout << "XXXXXXXXXXXXX" << std::endl;

	}


	//! pointer to the beginning of the buffer where the data is contained.
	T* _start;

	//! point to the first byte behind the memory currently used for storing elements.
	T* _finish;

	//! point to the first byte behind the memory allocated as buffer.
	T* _endOfStorage;

	//! The factor by which the capacity is increased or decreased if the
	//! threshold is encountered.
	double _capacityIncrement;

};

#endif /* DYNAMICARRAY_H_ */
